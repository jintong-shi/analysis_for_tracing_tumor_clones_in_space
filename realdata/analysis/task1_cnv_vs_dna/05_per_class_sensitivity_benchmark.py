#!/usr/bin/env python
"""Per-class gain/loss sensitivity benchmark for real SpatialCNV Task 1.

Input is the Step 05 metrics table
`<DATA_ROOT>/realdata/results/task1_cnv_vs_dna/pseudobulk_vs_wes_metrics_all_datasets.csv`,
which supplies each caller's pseudobulk prediction path (`pseudo_file`). Ding
ground truth is loaded from the sample-local `wes_gene_call_matrix.csv`;
metastatic/skin ground truth falls back to `truth_links/wes_cnv` segments
projected to genes with the same gene reference used by Step 05. The current
benchmark is gene-level, so "bin" in the output means one aligned gene/bin unit
from these files.

Figure inputs produced here (all under `--results-dir`):
  ding_gain_loss_per_sample_sensitivity.csv  -> figS_task1C_real_detail.R
  ding_gain_loss_tool_boxplot_stats.csv      -> figS_task1C_real_detail.R
The remaining outputs (`per_sample_sensitivity.csv`, `per_tool_wilcoxon.csv`,
`ding_gain_loss_overall_boxplot_stats.csv`, PNGs) are the working/analysis side
and do not feed the release figures.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = Path("<DATA_ROOT>/realdata/results/task1_cnv_vs_dna")

CONFIG = {
    "baseline": 0.0,
    "grid_quantiles_per_side": 50,
    "fdr_alpha": 0.05,
    "fdr_borderline": 0.10,
    "metrics_csv": RESULTS_DIR / "pseudobulk_vs_wes_metrics_all_datasets.csv",
    "sample_overview_csv": REPO_ROOT / "figures" / "utils" / "sample_overview.csv",
    "gene_pos_file": Path("<REF_ROOT>/gencode_v38_gene_pos.txt"),
    "results_dir": RESULTS_DIR / "per_class_sensitivity",
    "run_priority": (
        "HistoTumor_Benign-ref",
        "Decon_StromaImm-ref",
        "AllSpots_no-ref",
        "Decon_no-ref",
    ),
    "cnv_negative_samples": ("skin_P4_V1", "skin_P4_V2"),
    "wes_track_selector_metric": "Spearman",
    "plot_dpi": 300,
    "stop_if_all_gt_zero_loss": True,
}

GT_LOSS = -1
GT_BASE = 0
GT_GAIN = 1
CLASSES = (GT_GAIN, GT_BASE, GT_LOSS)
MERGE_COORDS = ("chr", "start", "end")


def load_pair(gt_path: Path, pred_path: Path) -> pd.DataFrame:
    """Load and inner-join one ground-truth / prediction pair."""
    if not gt_path.exists():
        raise FileNotFoundError(f"Missing ground-truth file: {gt_path}")
    if not pred_path.exists():
        raise FileNotFoundError(f"Missing prediction file: {pred_path}")

    gt = pd.read_csv(gt_path)
    pred = pd.read_csv(pred_path)
    if "gt" not in gt.columns:
        raise ValueError(f"{gt_path} must contain a 'gt' column")
    if "score" not in pred.columns:
        raise ValueError(f"{pred_path} must contain a 'score' column")

    if "bin_id" in gt.columns and "bin_id" in pred.columns:
        join_cols = ["bin_id"]
    elif all(col in gt.columns for col in MERGE_COORDS) and all(
        col in pred.columns for col in MERGE_COORDS
    ):
        join_cols = list(MERGE_COORDS)
    else:
        raise ValueError(
            "Inputs must share 'bin_id' or the coordinate columns "
            f"{', '.join(MERGE_COORDS)}: {gt_path}, {pred_path}"
        )

    merged = gt[join_cols + ["gt"]].merge(
        pred[join_cols + ["score"]], on=join_cols, how="inner"
    )
    if merged.empty:
        raise ValueError(f"No overlapping bins after inner join: {gt_path}, {pred_path}")

    merged["gt"] = pd.to_numeric(merged["gt"], errors="raise").astype(int)
    merged["score"] = pd.to_numeric(merged["score"], errors="raise").astype(float)
    bad_gt = sorted(set(merged["gt"]) - {GT_LOSS, GT_BASE, GT_GAIN})
    if bad_gt:
        raise ValueError(f"{gt_path} has gt values outside {{-1,0,1}}: {bad_gt}")
    return merged


def classify(scores, lo: float, hi: float) -> np.ndarray:
    """Classify continuous scores with loss <= lo, gain >= hi, base otherwise."""
    scores = np.asarray(scores, dtype=float)
    pred = np.full(scores.shape, GT_BASE, dtype=int)
    pred[scores <= lo] = GT_LOSS
    pred[scores >= hi] = GT_GAIN
    return pred


def multiclass_f1(y_true, y_pred) -> float:
    """Macro mean of one-vs-rest F1 for gain/base/loss."""
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    f1_values = []
    for cls in CLASSES:
        true_pos = y_true == cls
        pred_pos = y_pred == cls
        tp = int(np.sum(true_pos & pred_pos))
        fp = int(np.sum(~true_pos & pred_pos))
        fn = int(np.sum(true_pos & ~pred_pos))
        denom = 2 * tp + fp + fn
        f1_values.append(0.0 if denom == 0 else (2 * tp) / denom)
    return float(np.mean(f1_values))


def _side_thresholds(scores: np.ndarray, baseline: float, side: str, n_grid: int) -> np.ndarray:
    if side == "loss":
        side_scores = scores[scores < baseline]
        fallback = np.nextafter(baseline, -np.inf)
    elif side == "gain":
        side_scores = scores[scores > baseline]
        fallback = np.nextafter(baseline, np.inf)
    else:
        raise ValueError(f"Unknown threshold side: {side}")

    if side_scores.size == 0:
        return np.array([fallback], dtype=float)

    probs = np.linspace(0, 1, max(2, n_grid))
    thresholds = np.unique(np.quantile(side_scores, probs))
    if side == "loss":
        thresholds = thresholds[thresholds < baseline]
    else:
        thresholds = thresholds[thresholds > baseline]
    if thresholds.size == 0:
        return np.array([fallback], dtype=float)
    return thresholds.astype(float)


def pick_thresholds(
    scores,
    y_true,
    baseline: float = CONFIG["baseline"],
    grid_quantiles_per_side: int = CONFIG["grid_quantiles_per_side"],
) -> tuple[float, float, float]:
    """Grid-search loss/gain thresholds maximizing macro gain/base/loss F1."""
    scores = np.asarray(scores, dtype=float)
    y_true = np.asarray(y_true, dtype=int)
    loss_thresholds = _side_thresholds(
        scores, baseline, "loss", grid_quantiles_per_side
    )
    gain_thresholds = _side_thresholds(
        scores, baseline, "gain", grid_quantiles_per_side
    )

    true_loss = y_true == GT_LOSS
    true_base = y_true == GT_BASE
    true_gain = y_true == GT_GAIN
    n_loss = int(np.sum(true_loss))
    n_base = int(np.sum(true_base))
    n_gain = int(np.sum(true_gain))

    pred_loss = scores[None, :] <= loss_thresholds[:, None]
    pred_gain = scores[None, :] >= gain_thresholds[:, None]

    # Counts by true class for each one-sided threshold. Loss/gain prediction
    # regions cannot overlap because thresholds are restricted around baseline.
    loss_on_loss = pred_loss[:, true_loss].sum(axis=1)
    loss_on_base = pred_loss[:, true_base].sum(axis=1)
    loss_on_gain = pred_loss[:, true_gain].sum(axis=1)
    gain_on_loss = pred_gain[:, true_loss].sum(axis=1)
    gain_on_base = pred_gain[:, true_base].sum(axis=1)
    gain_on_gain = pred_gain[:, true_gain].sum(axis=1)

    def f1_array(tp: np.ndarray, fp: np.ndarray, fn: np.ndarray) -> np.ndarray:
        denom = 2 * tp + fp + fn
        out = np.zeros_like(denom, dtype=float)
        np.divide(2 * tp, denom, out=out, where=denom != 0)
        return out

    f1_loss = f1_array(
        loss_on_loss,
        loss_on_base + loss_on_gain,
        n_loss - loss_on_loss,
    )
    f1_gain = f1_array(
        gain_on_gain,
        gain_on_loss + gain_on_base,
        n_gain - gain_on_gain,
    )

    tp_base = n_base - loss_on_base[:, None] - gain_on_base[None, :]
    fp_base = (
        (n_loss - loss_on_loss[:, None] - gain_on_loss[None, :])
        + (n_gain - loss_on_gain[:, None] - gain_on_gain[None, :])
    )
    fn_base = loss_on_base[:, None] + gain_on_base[None, :]
    f1_base = f1_array(tp_base, fp_base, fn_base)

    f1_grid = (f1_loss[:, None] + f1_base + f1_gain[None, :]) / 3
    best_f1 = float(np.nanmax(f1_grid))
    distance_grid = (
        (baseline - loss_thresholds)[:, None]
        + (gain_thresholds - baseline)[None, :]
    )
    tied = np.isclose(f1_grid, best_f1, rtol=0, atol=1e-12)
    distance_grid = np.where(tied, distance_grid, np.inf)
    best_i, best_j = np.unravel_index(np.argmin(distance_grid), distance_grid.shape)
    return float(loss_thresholds[best_i]), float(gain_thresholds[best_j]), best_f1


def per_class_sensitivity(y_true, y_pred) -> dict[str, float | int]:
    """Compute gain/loss sensitivity and true-class counts."""
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    true_gain = y_true == GT_GAIN
    true_loss = y_true == GT_LOSS
    n_gain = int(np.sum(true_gain))
    n_loss = int(np.sum(true_loss))
    sens_gain = (
        float(np.sum(true_gain & (y_pred == GT_GAIN)) / n_gain) if n_gain > 0 else np.nan
    )
    sens_loss = (
        float(np.sum(true_loss & (y_pred == GT_LOSS)) / n_loss) if n_loss > 0 else np.nan
    )
    return {
        "sensitivity_gain": sens_gain,
        "sensitivity_loss": sens_loss,
        "n_true_gain_bins": n_gain,
        "n_true_loss_bins": n_loss,
    }


def run_tool_test(tool_df: pd.DataFrame) -> dict[str, float | int | str]:
    """Run the paired Wilcoxon gain-vs-loss sensitivity test for one tool."""
    tool = str(tool_df["tool"].iloc[0])
    paired = tool_df.dropna(subset=["sensitivity_gain", "sensitivity_loss"]).copy()
    n_paired = int(len(paired))
    if n_paired == 0:
        return {
            "tool": tool,
            "n_paired_samples": 0,
            "median_gain_sensitivity": np.nan,
            "median_loss_sensitivity": np.nan,
            "median_diff": np.nan,
            "wilcoxon_stat": np.nan,
            "p_raw": np.nan,
        }

    gain = paired["sensitivity_gain"].to_numpy(dtype=float)
    loss = paired["sensitivity_loss"].to_numpy(dtype=float)
    diff = gain - loss
    if n_paired < 2:
        stat = np.nan
        p_raw = np.nan
    elif np.allclose(diff, 0.0, equal_nan=False):
        stat = 0.0
        p_raw = 1.0
    else:
        test = stats.wilcoxon(gain, loss, alternative="two-sided")
        stat = float(test.statistic)
        p_raw = float(test.pvalue)

    return {
        "tool": tool,
        "n_paired_samples": n_paired,
        "median_gain_sensitivity": float(np.median(gain)),
        "median_loss_sensitivity": float(np.median(loss)),
        "median_diff": float(np.median(diff)),
        "wilcoxon_stat": stat,
        "p_raw": p_raw,
    }


def _csv_files(path: Path) -> list[Path]:
    return sorted(path.glob("*.csv")) if path.exists() else []


def _read_gt_balance(gt_path: Path) -> dict[str, int | str]:
    gt = pd.read_csv(gt_path)
    if "gt" not in gt.columns:
        raise ValueError(f"{gt_path} must contain a 'gt' column")
    values = pd.to_numeric(gt["gt"], errors="raise").astype(int)
    bad_gt = sorted(set(values) - {GT_LOSS, GT_BASE, GT_GAIN})
    if bad_gt:
        raise ValueError(f"{gt_path} has gt values outside {{-1,0,1}}: {bad_gt}")
    counts = values.value_counts()
    return {
        "sample": gt_path.stem,
        "n_loss": int(counts.get(GT_LOSS, 0)),
        "n_base": int(counts.get(GT_BASE, 0)),
        "n_gain": int(counts.get(GT_GAIN, 0)),
        "n_bins": int(len(values)),
    }


def discover_inventory(gt_dir: Path, pred_dir: Path) -> tuple[list[dict[str, Path | str]], pd.DataFrame]:
    """Discover runnable (tool, sample) pairs and print inventory diagnostics."""
    print(f"[inventory] ground truth dir: {gt_dir}")
    print(f"[inventory] predictions dir: {pred_dir}")
    if not gt_dir.exists():
        raise FileNotFoundError(f"Missing ground-truth directory: {gt_dir}")
    if not pred_dir.exists():
        raise FileNotFoundError(f"Missing predictions directory: {pred_dir}")

    gt_files = {path.stem: path for path in _csv_files(gt_dir)}
    if not gt_files:
        raise FileNotFoundError(f"No ground-truth CSV files found in {gt_dir}")

    balances = pd.DataFrame([_read_gt_balance(path) for path in gt_files.values()])
    balances = balances.sort_values("sample").reset_index(drop=True)
    print("\n[class balance]")
    print(balances.to_string(index=False))
    zero_loss_samples = balances.loc[balances["n_loss"] == 0, "sample"].tolist()
    if zero_loss_samples:
        print(
            "\n[class balance] samples with zero true-loss bins "
            f"(loss sensitivity will be NaN): {', '.join(zero_loss_samples)}"
        )
    if CONFIG["stop_if_all_gt_zero_loss"] and int(balances["n_loss"].sum()) == 0:
        raise RuntimeError(
            "All discovered ground-truth files have zero loss bins; "
            "loss-sensitivity analysis is meaningless."
        )

    tool_dirs = sorted(path for path in pred_dir.iterdir() if path.is_dir())
    if not tool_dirs:
        raise FileNotFoundError(f"No tool prediction directories found in {pred_dir}")

    inventory = []
    unpaired_predictions = []
    for tool_dir in tool_dirs:
        pred_files = {path.stem: path for path in _csv_files(tool_dir)}
        for sample, pred_path in sorted(pred_files.items()):
            gt_path = gt_files.get(sample)
            if gt_path is None:
                unpaired_predictions.append((tool_dir.name, sample, pred_path))
                continue
            inventory.append(
                {
                    "tool": tool_dir.name,
                    "sample": sample,
                    "gt_path": gt_path,
                    "pred_path": pred_path,
                }
            )

    print("\n[inventory] discovered runnable pairs:")
    if inventory:
        inv_df = pd.DataFrame(
            [{"tool": row["tool"], "sample": row["sample"]} for row in inventory]
        )
        counts = inv_df.groupby("tool", sort=True)["sample"].nunique().reset_index()
        counts = counts.rename(columns={"sample": "n_samples"})
        print(counts.to_string(index=False))
    else:
        print("(none)")

    if unpaired_predictions:
        print("\n[inventory] prediction files skipped because matching GT is missing:")
        for tool, sample, pred_path in unpaired_predictions:
            print(f"  tool={tool} sample={sample} path={pred_path}")

    missing_by_tool = []
    gt_sample_set = set(gt_files)
    for tool_dir in tool_dirs:
        pred_sample_set = {path.stem for path in _csv_files(tool_dir)}
        missing = sorted(gt_sample_set - pred_sample_set)
        if missing:
            missing_by_tool.append((tool_dir.name, missing))
    if missing_by_tool:
        print("\n[inventory] GT samples without predictions for each tool:")
        for tool, missing in missing_by_tool:
            print(f"  {tool}: {', '.join(missing)}")

    if not inventory:
        raise RuntimeError("No runnable (tool, sample) pairs were discovered.")
    return inventory, balances


def _read_gene_value_csv(path: Path, value_name: str) -> pd.DataFrame:
    data = pd.read_csv(path)
    if data.empty:
        raise ValueError(f"Empty CSV: {path}")
    gene_col = "gene" if "gene" in data.columns else data.columns[0]
    value_cols = [col for col in data.columns if col != gene_col]
    if not value_cols:
        raise ValueError(f"{path} must contain a value column besides gene")
    out = data[[gene_col, value_cols[0]]].copy()
    out.columns = ["gene", value_name]
    out["gene"] = _normalize_gene_ids(out["gene"])
    out[value_name] = pd.to_numeric(out[value_name], errors="coerce")
    out = out.dropna(subset=["gene", value_name])
    out = out[np.isfinite(out[value_name].to_numpy(dtype=float))]
    return out.drop_duplicates(subset="gene", keep="first")


def _normalize_gene_ids(values) -> pd.Series:
    return pd.Series(values, copy=False).astype(str).str.strip().str.upper()


def _normalize_key(value: str) -> str:
    return "".join(ch for ch in str(value).strip().lower() if ch.isalnum())


def _normalize_chr(values) -> pd.Series:
    out = pd.Series(values, copy=False).astype(str).str.strip()
    out = out.str.replace(r"^chr", "", case=False, regex=True)
    out = out.replace({"M": "M", "MT": "M", "MtDNA": "M", "mitochondria": "M"})
    return "chr" + out


def _normalize_call_to_num(values) -> pd.Series:
    raw = pd.Series(values, copy=False).astype(str).str.strip().str.upper()
    mapped = raw.replace(
        {
            "GAIN": "1",
            "AMP": "1",
            "AMPLIFICATION": "1",
            "+": "1",
            "+1": "1",
            "LOSS": "-1",
            "DEL": "-1",
            "DELETION": "-1",
            "-": "-1",
            "NEUTRAL": "0",
            "DIPLOID": "0",
            "+0": "0",
        }
    )
    num = pd.to_numeric(mapped, errors="coerce")
    return pd.Series(np.sign(num), index=raw.index).replace({np.nan: np.nan})


def _choose_col(columns, candidates, required: bool, label: str) -> str | None:
    key_to_col = {_normalize_key(col): col for col in columns}
    for candidate in candidates:
        col = key_to_col.get(_normalize_key(candidate))
        if col is not None:
            return col
    if required:
        raise ValueError(f"Cannot find {label} column. Tried: {', '.join(candidates)}")
    return None


_GENE_REF_CACHE: pd.DataFrame | None = None


def _gene_reference() -> pd.DataFrame:
    global _GENE_REF_CACHE
    if _GENE_REF_CACHE is not None:
        return _GENE_REF_CACHE
    path = Path(CONFIG["gene_pos_file"])
    if not path.exists():
        raise FileNotFoundError(f"Missing gene position file: {path}")
    gene_ref = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3])
    gene_ref.columns = ["gene", "chr", "start", "end"]
    gene_ref["gene"] = _normalize_gene_ids(gene_ref["gene"])
    gene_ref["chr"] = _normalize_chr(gene_ref["chr"])
    gene_ref["start"] = pd.to_numeric(gene_ref["start"], errors="coerce").round()
    gene_ref["end"] = pd.to_numeric(gene_ref["end"], errors="coerce").round()
    gene_ref = gene_ref.dropna(subset=["gene", "chr", "start", "end"])
    gene_ref["start"] = gene_ref["start"].astype(int)
    gene_ref["end"] = gene_ref["end"].astype(int)
    gene_ref = gene_ref[gene_ref["end"] >= gene_ref["start"]].copy()
    gene_ref = gene_ref.drop_duplicates(subset="gene", keep="first")
    gene_ref["gene_len"] = gene_ref["end"] - gene_ref["start"] + 1
    _GENE_REF_CACHE = gene_ref.sort_values(["chr", "start", "end", "gene"]).reset_index(drop=True)
    return _GENE_REF_CACHE


def _infer_wes_call_path(pseudo_file: Path) -> Path:
    # pseudo_file = run_links/{sample}/{run}/{tool}/pseudobulk_*.csv
    if len(pseudo_file.parents) < 3:
        raise ValueError(f"Cannot infer sample directory from pseudo_file: {pseudo_file}")
    return pseudo_file.parents[2] / "wes_gene_call_matrix.csv"


def _dataset_dir_from_pred_path(pred_path: Path, dataset: str) -> Path:
    for parent in pred_path.parents:
        if parent.name == dataset:
            return parent
    raise ValueError(f"Cannot infer dataset dir for {dataset} from {pred_path}")


def _truth_segment_candidates(dataset_dir: Path, dataset: str, sample: str) -> list[Path]:
    truth_dir = dataset_dir / "truth_links" / "wes_cnv"
    if not truth_dir.exists():
        return []
    direct = [
        truth_dir / f"{sample}_wes_cnv.csv",
        truth_dir / f"{sample}_wes_segments_results.csv",
        truth_dir / sample,
        truth_dir / f"{sample}.csv",
    ]
    mapped = []
    if dataset == "skin_tumor_ff":
        patient = sample
        if patient.startswith("skin_"):
            patient = patient.split("_")[1]
        patient = patient.lower().replace("p", "pt", 1)
        mapped = [
            truth_dir / f"{patient}_wes_segments_results.csv",
            truth_dir / f"{patient}_wes_cnv.csv",
        ]
    out = []
    for path in direct + mapped:
        if path.exists() and path not in out:
            out.append(path)
    return out


def _project_segment_truth(seg_file: Path, wes_sample_id: str) -> pd.DataFrame:
    seg = pd.read_csv(seg_file)
    if seg.empty:
        raise ValueError(f"Empty segment truth file: {seg_file}")
    col_chr = _choose_col(seg.columns, ["CONTIG", "chr", "chromosome", "chrom"], True, "chromosome")
    col_start = _choose_col(seg.columns, ["START", "start", "POS_START"], True, "start")
    col_end = _choose_col(seg.columns, ["END", "end", "POS_END"], True, "end")
    col_call = _choose_col(seg.columns, ["CALL", "call"], True, "call")
    col_wes = _choose_col(
        seg.columns,
        ["sample", "SAMPLE", "wes_sample", "sample_id", "wes_sample_id"],
        False,
        "wes_sample",
    )

    seg = seg.copy()
    if col_wes is None:
        seg["sample_id"] = seg_file.stem
    else:
        seg["sample_id"] = seg[col_wes].astype(str)
    if wes_sample_id in set(seg["sample_id"]):
        seg = seg[seg["sample_id"].eq(wes_sample_id)].copy()
    elif seg["sample_id"].nunique(dropna=True) == 1:
        # Existing Step 05 accepts one-column truth even when the selected WES id
        # is an alias for the same sample.
        pass
    else:
        raise ValueError(f"{seg_file} has no segment rows for WES sample {wes_sample_id}")

    seg["chr"] = _normalize_chr(seg[col_chr])
    seg["start"] = pd.to_numeric(seg[col_start], errors="coerce").round()
    seg["end"] = pd.to_numeric(seg[col_end], errors="coerce").round()
    seg["call_num"] = _normalize_call_to_num(seg[col_call])
    seg = seg.dropna(subset=["chr", "start", "end", "call_num"]).copy()
    seg["start"] = seg["start"].astype(int)
    seg["end"] = seg["end"].astype(int)
    seg = seg[seg["end"] >= seg["start"]]
    if seg.empty:
        raise ValueError(f"No valid segments after parsing: {seg_file}")

    gene_ref = _gene_reference()
    rows = []
    for chrom, seg_chr in seg.groupby("chr", sort=False):
        genes = gene_ref[gene_ref["chr"].eq(chrom)]
        if genes.empty:
            continue
        for _, srow in seg_chr.iterrows():
            overlap = genes[
                (genes["start"] <= int(srow["end"]))
                & (genes["end"] >= int(srow["start"]))
            ].copy()
            if overlap.empty:
                continue
            ovlen = (
                np.minimum(overlap["end"].to_numpy(), int(srow["end"]))
                - np.maximum(overlap["start"].to_numpy(), int(srow["start"]))
                + 1
            )
            keep = ovlen > 0
            if not np.any(keep):
                continue
            rows.append(
                pd.DataFrame(
                    {
                        "gene": overlap.loc[keep, "gene"].to_numpy(),
                        "weighted_call": float(srow["call_num"]) * ovlen[keep],
                    }
                )
            )
    if not rows:
        raise ValueError(f"No gene overlaps after projecting segments: {seg_file}")
    weighted = pd.concat(rows, ignore_index=True)
    agg = weighted.groupby("gene", as_index=False)["weighted_call"].sum()
    gene_len = gene_ref[["gene", "gene_len"]]
    agg = agg.merge(gene_len, on="gene", how="left")
    agg["gt"] = np.sign(agg["weighted_call"] / agg["gene_len"]).astype(int)
    return agg[["gene", "gt"]]


def _load_truth_calls_for_row(row: dict[str, object]) -> pd.DataFrame:
    truth_path = Path(str(row["truth_file"]))
    wes_sample_id = str(row["wes_sample_id"])
    if truth_path.exists():
        truth = pd.read_csv(truth_path)
        if "gene" not in truth.columns:
            raise ValueError(f"{truth_path} must contain a 'gene' column")
        if wes_sample_id not in truth.columns:
            raise ValueError(
                f"{truth_path} does not contain WES sample column '{wes_sample_id}'"
            )
        truth = truth[["gene", wes_sample_id]].copy()
        truth.columns = ["gene", "gt"]
        truth["gene"] = _normalize_gene_ids(truth["gene"])
        truth["gt"] = pd.to_numeric(truth["gt"], errors="coerce")
        truth = truth.dropna(subset=["gene", "gt"]).drop_duplicates(
            subset="gene", keep="first"
        )
        truth = truth[np.isfinite(truth["gt"].to_numpy(dtype=float))]
        truth["gt"] = np.sign(truth["gt"]).astype(int)
        return truth

    pred_path = Path(str(row["pred_path"]))
    dataset = str(row["dataset"])
    sample = str(row["sample"])
    dataset_dir = _dataset_dir_from_pred_path(pred_path, dataset)
    candidates = _truth_segment_candidates(dataset_dir, dataset, sample)
    for seg_file in candidates:
        try:
            truth = _project_segment_truth(seg_file, wes_sample_id)
            row["truth_file"] = str(seg_file)
            row["truth_source"] = "segment_projection"
            return truth
        except Exception:
            continue
    raise FileNotFoundError(
        f"No usable WES truth for {dataset}/{sample}/{wes_sample_id}; "
        f"matrix path tried: {truth_path}; segment candidates: {candidates}"
    )


def load_real_task1_pair(row: dict[str, object]) -> pd.DataFrame:
    """Load one actual Real Task 1 WES-call / pseudobulk-gene pair."""
    pred_path = Path(str(row["pred_path"]))
    if not pred_path.exists():
        raise FileNotFoundError(f"Missing pseudobulk prediction file: {pred_path}")

    pred = _read_gene_value_csv(pred_path, "score")
    truth = _load_truth_calls_for_row(row)

    merged = truth.merge(pred, on="gene", how="inner")
    if merged.empty:
        raise ValueError(
            "No overlapping genes after inner join: "
            f"{truth_path} [{wes_sample_id}], {pred_path}"
        )
    merged["gt"] = merged["gt"].astype(int)
    bad_gt = sorted(set(merged["gt"]) - {GT_LOSS, GT_BASE, GT_GAIN})
    if bad_gt:
        raise ValueError(
            f"{truth_path} [{wes_sample_id}] has gt values outside {{-1,0,1}}: "
            f"{bad_gt}"
        )
    return merged[["gene", "gt", "score"]]


def _load_real_overview(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing sample overview CSV: {path}")
    overview = pd.read_csv(path)
    need = {"sample", "sample_type", "histology_available"}
    missing = need - set(overview.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
    overview = overview.copy()
    overview["sample"] = overview["sample"].astype(str).str.strip()
    overview["sample_type_norm"] = np.where(
        overview["sample_type"].astype(str).str.strip().str.lower().eq("ffpe"),
        "ffpe",
        "ff",
    )
    overview["histology_available_bool"] = (
        overview["histology_available"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("yes")
    )
    return overview[["sample", "sample_type_norm", "histology_available_bool"]]


def _select_real_wes_tracks(metrics: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Choose one WES track per ST sample from existing Real Task 1 metrics."""
    if metric not in metrics.columns:
        raise ValueError(f"WES selector metric '{metric}' not in metrics table")
    scored = metrics.copy()
    scored["_selector_score"] = pd.to_numeric(scored[metric], errors="coerce")
    if scored["_selector_score"].isna().all() and metric != "max_F1":
        scored["_selector_score"] = pd.to_numeric(scored["max_F1"], errors="coerce")
    tracks = (
        scored.groupby(["dataset", "sample", "wes_sample_id"], dropna=False)
        .agg(selector_score=("_selector_score", "mean"), n_rows=("tool", "size"))
        .reset_index()
    )
    tracks = tracks.sort_values(
        ["dataset", "sample", "selector_score", "n_rows", "wes_sample_id"],
        ascending=[True, True, False, False, True],
        na_position="last",
    )
    return tracks.groupby(["dataset", "sample"], as_index=False).head(1)


def _real_run_priority(row: pd.Series) -> int:
    run = str(row["run"])
    histology = bool(row["histology_available_bool"])
    if histology and run == "HistoTumor_Benign-ref":
        return 1
    if histology and run == "Decon_StromaImm-ref":
        return 2
    if not histology and run == "Decon_StromaImm-ref":
        return 1
    if run == "AllSpots_no-ref":
        return 3
    if run == "Decon_no-ref":
        return 4
    return 9


def discover_real_task1_inventory(
    metrics_csv: Path,
    sample_overview_csv: Path,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    """Discover actual Real Task 1 rows from the combined metrics table."""
    print(f"[inventory] metrics CSV: {metrics_csv}")
    print(f"[inventory] sample overview CSV: {sample_overview_csv}")
    if not metrics_csv.exists():
        raise FileNotFoundError(f"Missing metrics CSV: {metrics_csv}")

    raw = pd.read_csv(metrics_csv)
    need = {
        "dataset",
        "sample",
        "run",
        "tool",
        "wes_sample_id",
        "status",
        "pseudo_file",
        "max_F1",
        "Spearman",
    }
    missing = need - set(raw.columns)
    if missing:
        raise ValueError(f"{metrics_csv} missing columns: {', '.join(sorted(missing))}")

    overview = _load_real_overview(sample_overview_csv)
    dat = raw[raw["status"].astype(str).eq("ok")].copy()
    dat["sample"] = dat["sample"].astype(str).str.strip()
    dat["dataset"] = dat["dataset"].astype(str).str.strip()
    dat["run"] = dat["run"].astype(str).str.strip()
    dat["tool"] = dat["tool"].astype(str).str.strip()
    dat["wes_sample_id"] = dat["wes_sample_id"].astype(str).str.strip()
    dat = dat.merge(overview, on="sample", how="inner")
    dat = dat[~dat["sample"].isin(CONFIG["cnv_negative_samples"])].copy()
    if dat.empty:
        raise RuntimeError("No Real Task 1 ok rows remain after scope filtering")

    chosen_tracks = _select_real_wes_tracks(
        dat, metric=str(CONFIG["wes_track_selector_metric"])
    )
    print("\n[inventory] selected WES track per sample:")
    print(
        chosen_tracks[
            ["dataset", "sample", "wes_sample_id", "selector_score", "n_rows"]
        ].to_string(index=False)
    )
    dat = dat.merge(
        chosen_tracks[["dataset", "sample", "wes_sample_id"]],
        on=["dataset", "sample", "wes_sample_id"],
        how="inner",
    )
    dat["sel_pri"] = dat.apply(_real_run_priority, axis=1)
    dat["max_F1_num"] = pd.to_numeric(dat["max_F1"], errors="coerce")
    dat["Spearman_num"] = pd.to_numeric(dat["Spearman"], errors="coerce")
    dat = dat.sort_values(
        ["dataset", "sample", "tool", "sel_pri", "max_F1_num", "Spearman_num"],
        ascending=[True, True, True, True, False, False],
        na_position="last",
    )
    best = dat.groupby(["dataset", "sample", "tool"], as_index=False).head(1).copy()
    best["pred_path"] = best["pseudo_file"].map(lambda value: str(Path(str(value))))
    best["truth_file"] = best["pred_path"].map(
        lambda value: str(_infer_wes_call_path(Path(value)))
    )

    missing_paths = []
    for _, row in best.iterrows():
        if not Path(row["pred_path"]).exists():
            missing_paths.append(("prediction", row["pred_path"]))
    if missing_paths:
        print("\n[inventory] missing files:")
        for kind, path in missing_paths:
            print(f"  {kind}: {path}")
        raise FileNotFoundError("Missing files in discovered Real Task 1 inventory")

    balances = []
    for _, row in best[["dataset", "sample", "wes_sample_id", "truth_file"]].drop_duplicates().iterrows():
        row_dict = row.to_dict()
        sample_rows = best[
            (best["dataset"].eq(row["dataset"]))
            & (best["sample"].eq(row["sample"]))
            & (best["wes_sample_id"].eq(row["wes_sample_id"]))
        ]
        # Any selected prediction from this sample gives the dataset/sample path
        # needed for segment fallback.
        row_dict["pred_path"] = sample_rows.iloc[0]["pred_path"]
        truth = _load_truth_calls_for_row(row_dict)
        values = pd.to_numeric(truth["gt"], errors="coerce").dropna().astype(int)
        bad_gt = sorted(set(values) - {GT_LOSS, GT_BASE, GT_GAIN})
        if bad_gt:
            raise ValueError(
                f"{row_dict['truth_file']} [{row['wes_sample_id']}] has gt values outside {{-1,0,1}}: "
                f"{bad_gt}"
            )
        counts = values.value_counts()
        balances.append(
            {
                "dataset": row["dataset"],
                "sample": row["sample"],
                "wes_sample_id": row["wes_sample_id"],
                "truth_source": row_dict.get("truth_source", "wes_gene_call_matrix"),
                "truth_file": row_dict["truth_file"],
                "n_loss": int(counts.get(GT_LOSS, 0)),
                "n_base": int(counts.get(GT_BASE, 0)),
                "n_gain": int(counts.get(GT_GAIN, 0)),
                "n_bins": int(len(values)),
            }
        )
    balances_df = pd.DataFrame(balances).sort_values(["dataset", "sample"])
    print("\n[class balance: selected WES gene-call tracks]")
    print(balances_df.to_string(index=False))
    zero_loss_samples = balances_df.loc[balances_df["n_loss"] == 0, "sample"].tolist()
    if zero_loss_samples:
        print(
            "\n[class balance] selected samples with zero true-loss genes "
            f"(loss sensitivity will be NaN): {', '.join(zero_loss_samples)}"
        )
    if CONFIG["stop_if_all_gt_zero_loss"] and int(balances_df["n_loss"].sum()) == 0:
        raise RuntimeError(
            "All selected Real Task 1 WES tracks have zero loss genes; "
            "loss-sensitivity analysis is meaningless."
        )

    print("\n[inventory] selected runnable pairs by tool:")
    counts = best.groupby("tool", sort=True)["sample"].nunique().reset_index()
    counts = counts.rename(columns={"sample": "n_samples"})
    print(counts.to_string(index=False))

    inventory = []
    for _, row in best.iterrows():
        inventory.append(
            {
                "dataset": row["dataset"],
                "sample": row["sample"],
                "sample_type": row["sample_type_norm"],
                "run": row["run"],
                "tool": row["tool"],
                "wes_sample_id": row["wes_sample_id"],
                "pred_path": row["pred_path"],
                "truth_file": row["truth_file"],
            }
        )
    return inventory, balances_df


def run_pair(row: dict[str, Path | str], baseline: float, grid_n: int) -> dict[str, float | int | str]:
    if "gt_path" in row:
        data = load_pair(Path(row["gt_path"]), Path(row["pred_path"]))
    else:
        data = load_real_task1_pair(row)
    lo, hi, max_f1 = pick_thresholds(
        data["score"].to_numpy(),
        data["gt"].to_numpy(),
        baseline=baseline,
        grid_quantiles_per_side=grid_n,
    )
    y_pred = classify(data["score"].to_numpy(), lo, hi)
    sens = per_class_sensitivity(data["gt"].to_numpy(), y_pred)
    out = {
        "tool": str(row["tool"]),
        "sample": str(row["sample"]),
        "chosen_loss_threshold": lo,
        "chosen_gain_threshold": hi,
        "max_multiclass_f1": max_f1,
        "n_overlap_bins": int(len(data)),
        **sens,
    }
    for key in ("dataset", "sample_type", "run", "wes_sample_id", "truth_file", "pred_path"):
        if key in row:
            out[key] = str(row[key])
    return out


def adjust_pvalues(tool_results: pd.DataFrame) -> pd.DataFrame:
    tool_results = tool_results.copy()
    tool_results["p_bh"] = np.nan
    valid = tool_results["p_raw"].notna()
    if valid.any():
        _, p_bh, _, _ = multipletests(
            tool_results.loc[valid, "p_raw"].to_numpy(dtype=float),
            method="fdr_bh",
        )
        tool_results.loc[valid, "p_bh"] = p_bh
    return tool_results


def plot_gain_vs_loss(per_sample: pd.DataFrame, out_path: Path, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paired = per_sample.dropna(subset=["sensitivity_gain", "sensitivity_loss"]).copy()
    tools = sorted(per_sample["tool"].unique())
    n_tools = len(tools)
    n_cols = min(4, max(1, n_tools))
    n_rows = math.ceil(n_tools / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(3.2 * n_cols, 3.0 * n_rows),
        sharey=True,
        squeeze=False,
    )

    for index, tool in enumerate(tools):
        ax = axes[index // n_cols][index % n_cols]
        sub = paired[paired["tool"] == tool].sort_values("sample")
        for _, row in sub.iterrows():
            ax.plot(
                [0, 1],
                [row["sensitivity_loss"], row["sensitivity_gain"]],
                color="0.70",
                linewidth=0.8,
                zorder=1,
            )
            ax.scatter(
                [0, 1],
                [row["sensitivity_loss"], row["sensitivity_gain"]],
                color=["#2b6cb0", "#c53030"],
                s=18,
                zorder=2,
            )
        ax.set_title(tool, fontsize=10)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Loss", "Gain"])
        ax.set_ylim(-0.03, 1.03)
        ax.grid(axis="y", color="0.90", linewidth=0.6)
        if index % n_cols == 0:
            ax.set_ylabel("Sensitivity")
        ax.text(
            0.5,
            0.02,
            f"n={len(sub)}",
            ha="center",
            va="bottom",
            transform=ax.transAxes,
            fontsize=8,
            color="0.35",
        )

    for index in range(n_tools, n_rows * n_cols):
        axes[index // n_cols][index % n_cols].axis("off")

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_conclusion_summary(tool_results: pd.DataFrame, out_path: Path, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    data = tool_results.copy()
    data = data[data["n_paired_samples"].fillna(0).astype(int) > 0].copy()
    data = data.sort_values("median_diff", ascending=True)
    y_pos = np.arange(len(data))

    colors = []
    for _, row in data.iterrows():
        diff = row["median_diff"]
        p_bh = row["p_bh"]
        if pd.notna(p_bh) and p_bh < CONFIG["fdr_alpha"] and diff > 0:
            colors.append("#b2182b")  # significant gain > loss
        elif pd.notna(p_bh) and p_bh < CONFIG["fdr_alpha"] and diff < 0:
            colors.append("#2166ac")  # significant loss > gain
        else:
            colors.append("0.70")

    fig_h = max(4.2, 0.36 * len(data) + 1.6)
    fig, ax = plt.subplots(figsize=(7.8, fig_h))
    ax.barh(y_pos, data["median_diff"], color=colors, edgecolor="0.25", linewidth=0.4)
    ax.axvline(0, color="0.15", linewidth=0.9)
    ax.axvline(0.1, color="0.85", linewidth=0.6, zorder=0)
    ax.axvline(-0.1, color="0.85", linewidth=0.6, zorder=0)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(data["tool"])
    ax.set_xlabel("Median sensitivity difference (gain - loss)")
    ax.set_title("No BH-significant gain bias; copyVAE shows loss-biased sensitivity")
    ax.set_xlim(
        min(-0.45, float(np.nanmin(data["median_diff"])) - 0.12),
        max(0.18, float(np.nanmax(data["median_diff"])) + 0.06),
    )
    ax.grid(axis="x", color="0.90", linewidth=0.6)
    ax.set_axisbelow(True)

    for y, (_, row) in zip(y_pos, data.iterrows()):
        diff = float(row["median_diff"])
        n = int(row["n_paired_samples"])
        p_text = "BH=NA" if pd.isna(row["p_bh"]) else f"BH={row['p_bh']:.3g}"
        label = f"n={n}, {p_text}"
        x_text = diff + (0.012 if diff >= 0 else -0.012)
        ha = "left" if diff >= 0 else "right"
        ax.text(x_text, y, label, va="center", ha=ha, fontsize=8, color="0.20")

    legend_handles = [
        Patch(facecolor="#b2182b", edgecolor="0.25", label="BH<0.05, gain > loss"),
        Patch(facecolor="#2166ac", edgecolor="0.25", label="BH<0.05, loss > gain"),
        Patch(facecolor="0.70", edgecolor="0.25", label="Not BH-significant"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _format_p(value: float, prefix: str = "p") -> str:
    if pd.isna(value):
        return f"{prefix}=NA"
    if value < 0.001:
        return f"{prefix}<0.001"
    return f"{prefix}={value:.3g}"


def _ding_tool_stats(ding: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tool, sub in ding.groupby("tool", sort=True):
        paired = sub.dropna(subset=["sensitivity_gain", "sensitivity_loss"])
        if paired.empty:
            rows.append(
                {
                    "tool": tool,
                    "n_paired_samples": 0,
                    "median_gain_sensitivity": np.nan,
                    "median_loss_sensitivity": np.nan,
                    "median_diff": np.nan,
                    "wilcoxon_stat": np.nan,
                    "p_raw": np.nan,
                }
            )
            continue
        gain = paired["sensitivity_gain"].to_numpy(dtype=float)
        loss = paired["sensitivity_loss"].to_numpy(dtype=float)
        diff = gain - loss
        if len(paired) < 2:
            stat = np.nan
            p_raw = np.nan
        elif np.allclose(diff, 0.0, equal_nan=False):
            stat = 0.0
            p_raw = 1.0
        else:
            test = stats.wilcoxon(gain, loss, alternative="two-sided")
            stat = float(test.statistic)
            p_raw = float(test.pvalue)
        rows.append(
            {
                "tool": tool,
                "n_paired_samples": int(len(paired)),
                "median_gain_sensitivity": float(np.median(gain)),
                "median_loss_sensitivity": float(np.median(loss)),
                "median_diff": float(np.median(diff)),
                "wilcoxon_stat": stat,
                "p_raw": p_raw,
            }
        )
    return adjust_pvalues(pd.DataFrame(rows))


def plot_ding_tool_boxplot(
    per_sample: pd.DataFrame,
    out_path: Path,
    stats_path: Path,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ding = per_sample[per_sample["dataset"].eq("Ding_nature_wes")].copy()
    stats_df = _ding_tool_stats(ding)
    stats_df.to_csv(stats_path, index=False)

    stats_plot = stats_df[stats_df["n_paired_samples"] > 0].copy()
    stats_plot = stats_plot.sort_values("median_diff", ascending=True)
    tools = stats_plot["tool"].tolist()
    positions = np.arange(len(tools))
    gain_data = []
    loss_data = []
    for tool in tools:
        paired = ding[ding["tool"].eq(tool)].dropna(
            subset=["sensitivity_gain", "sensitivity_loss"]
        )
        gain_data.append(paired["sensitivity_gain"].to_numpy(dtype=float))
        loss_data.append(paired["sensitivity_loss"].to_numpy(dtype=float))

    fig, ax = plt.subplots(figsize=(max(10, 0.85 * len(tools) + 3), 5.4))
    width = 0.30
    gain_pos = positions - width / 1.5
    loss_pos = positions + width / 1.5
    common = dict(widths=width, patch_artist=True, showfliers=False)
    bp_gain = ax.boxplot(gain_data, positions=gain_pos, **common)
    bp_loss = ax.boxplot(loss_data, positions=loss_pos, **common)
    for patch in bp_gain["boxes"]:
        patch.set_facecolor("#c44e52")
        patch.set_alpha(0.75)
        patch.set_edgecolor("0.25")
    for patch in bp_loss["boxes"]:
        patch.set_facecolor("#4c72b0")
        patch.set_alpha(0.75)
        patch.set_edgecolor("0.25")
    for element in ("whiskers", "caps", "medians"):
        for line in bp_gain[element] + bp_loss[element]:
            line.set_color("0.25")
            line.set_linewidth(0.8)

    rng = np.random.default_rng(0)
    for x, values in zip(gain_pos, gain_data):
        if len(values):
            ax.scatter(
                x + rng.uniform(-0.045, 0.045, len(values)),
                values,
                s=14,
                color="#7f1d1d",
                alpha=0.35,
                linewidth=0,
                zorder=2,
            )
    for x, values in zip(loss_pos, loss_data):
        if len(values):
            ax.scatter(
                x + rng.uniform(-0.045, 0.045, len(values)),
                values,
                s=14,
                color="#1f4e79",
                alpha=0.35,
                linewidth=0,
                zorder=2,
            )

    y_bracket = 1.04
    y_text = 1.06
    for i, tool in enumerate(tools):
        row = stats_plot[stats_plot["tool"].eq(tool)].iloc[0]
        ax.plot(
            [gain_pos[i], gain_pos[i], loss_pos[i], loss_pos[i]],
            [y_bracket - 0.012, y_bracket, y_bracket, y_bracket - 0.012],
            color="0.25",
            linewidth=0.7,
        )
        ax.text(
            positions[i],
            y_text,
            _format_p(row["p_bh"], prefix="BH"),
            ha="center",
            va="bottom",
            fontsize=7,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(tools, rotation=35, ha="right")
    ax.set_ylabel("Sensitivity")
    ax.set_title("Ding-only per-class sensitivity by tool")
    ax.set_ylim(-0.03, 1.13)
    ax.grid(axis="y", color="0.90", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color="#c44e52", lw=8, alpha=0.75, label="Gain"),
            plt.Line2D([0], [0], color="#4c72b0", lw=8, alpha=0.75, label="Loss"),
        ],
        frameon=False,
        loc="lower left",
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_ding_overall_boxplot(
    per_sample: pd.DataFrame,
    out_path: Path,
    stats_path: Path,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ding = per_sample[per_sample["dataset"].eq("Ding_nature_wes")].copy()
    paired = ding.dropna(subset=["sensitivity_gain", "sensitivity_loss"])
    by_sample = (
        paired.groupby("sample", as_index=False)
        .agg(
            sensitivity_gain=("sensitivity_gain", "median"),
            sensitivity_loss=("sensitivity_loss", "median"),
            n_tools=("tool", "nunique"),
        )
        .dropna(subset=["sensitivity_gain", "sensitivity_loss"])
    )
    gain = by_sample["sensitivity_gain"].to_numpy(dtype=float)
    loss = by_sample["sensitivity_loss"].to_numpy(dtype=float)
    diff = gain - loss
    if len(by_sample) < 2:
        stat = np.nan
        p_raw = np.nan
    elif np.allclose(diff, 0.0, equal_nan=False):
        stat = 0.0
        p_raw = 1.0
    else:
        test = stats.wilcoxon(gain, loss, alternative="two-sided")
        stat = float(test.statistic)
        p_raw = float(test.pvalue)
    stats_df = pd.DataFrame(
        [
            {
                "scope": "Ding sample-level median across tools",
                "n_paired_samples": int(len(by_sample)),
                "median_gain_sensitivity": float(np.median(gain)) if len(gain) else np.nan,
                "median_loss_sensitivity": float(np.median(loss)) if len(loss) else np.nan,
                "median_diff": float(np.median(diff)) if len(diff) else np.nan,
                "wilcoxon_stat": stat,
                "p_raw": p_raw,
            }
        ]
    )
    stats_df.to_csv(stats_path, index=False)

    fig, ax = plt.subplots(figsize=(4.4, 4.8))
    bp = ax.boxplot(
        [gain, loss],
        positions=[0, 1],
        widths=0.42,
        patch_artist=True,
        showfliers=False,
    )
    for patch, color in zip(bp["boxes"], ["#c44e52", "#4c72b0"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor("0.25")
    for element in ("whiskers", "caps", "medians"):
        for line in bp[element]:
            line.set_color("0.25")
            line.set_linewidth(0.9)

    rng = np.random.default_rng(1)
    ax.scatter(
        rng.normal(0, 0.035, len(gain)),
        gain,
        s=20,
        color="#7f1d1d",
        alpha=0.45,
        linewidth=0,
        zorder=2,
    )
    ax.scatter(
        1 + rng.normal(0, 0.035, len(loss)),
        loss,
        s=20,
        color="#1f4e79",
        alpha=0.45,
        linewidth=0,
        zorder=2,
    )
    y = 1.04
    ax.plot([0, 0, 1, 1], [y - 0.015, y, y, y - 0.015], color="0.25", linewidth=0.8)
    ax.text(0.5, y + 0.025, _format_p(p_raw), ha="center", va="bottom", fontsize=10)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Gain", "Loss"])
    ax.set_ylabel("Sensitivity")
    ax.set_title("Ding-only overall sensitivity\nsample median across tools")
    ax.set_ylim(-0.03, 1.15)
    ax.grid(axis="y", color="0.90", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def print_summary(tool_results: pd.DataFrame, alpha: float, borderline: float) -> None:
    gain_better = tool_results["median_diff"] > 0
    significant = tool_results[gain_better & (tool_results["p_bh"] < alpha)]
    borderline_df = tool_results[
        gain_better
        & (tool_results["p_bh"] >= alpha)
        & (tool_results["p_bh"] < borderline)
    ]

    print("\n[summary] Gain sensitivity > loss sensitivity")
    if significant.empty:
        print(f"  Significant (BH < {alpha:g}): none")
    else:
        print(f"  Significant (BH < {alpha:g}):")
        for _, row in significant.sort_values("p_bh").iterrows():
            print(
                f"    {row['tool']}: p_bh={row['p_bh']:.4g}, "
                f"median_diff={row['median_diff']:.4f}, n={int(row['n_paired_samples'])}"
            )

    if borderline_df.empty:
        print(f"  Borderline ({alpha:g} <= BH < {borderline:g}): none")
    else:
        print(f"  Borderline ({alpha:g} <= BH < {borderline:g}):")
        for _, row in borderline_df.sort_values("p_bh").iterrows():
            print(
                f"    {row['tool']}: p_bh={row['p_bh']:.4g}, "
                f"median_diff={row['median_diff']:.4f}, n={int(row['n_paired_samples'])}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-csv", type=Path, default=CONFIG["metrics_csv"])
    parser.add_argument(
        "--sample-overview-csv", type=Path, default=CONFIG["sample_overview_csv"]
    )
    parser.add_argument("--results-dir", type=Path, default=CONFIG["results_dir"])
    parser.add_argument("--baseline", type=float, default=CONFIG["baseline"])
    parser.add_argument(
        "--grid-quantiles-per-side",
        type=int,
        default=CONFIG["grid_quantiles_per_side"],
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        inventory, _balances = discover_real_task1_inventory(
            args.metrics_csv, args.sample_overview_csv
        )
        args.results_dir.mkdir(parents=True, exist_ok=True)

        per_sample_rows = [
            run_pair(row, args.baseline, args.grid_quantiles_per_side)
            for row in inventory
        ]
        per_sample = pd.DataFrame(per_sample_rows).sort_values(["tool", "sample"])
        front_cols = [
            "dataset",
            "sample",
            "sample_type",
            "run",
            "tool",
            "wes_sample_id",
            "chosen_loss_threshold",
            "chosen_gain_threshold",
            "max_multiclass_f1",
            "sensitivity_gain",
            "sensitivity_loss",
            "n_true_gain_bins",
            "n_true_loss_bins",
            "n_overlap_bins",
            "truth_file",
            "pred_path",
        ]
        per_sample = per_sample[
            [col for col in front_cols if col in per_sample.columns]
            + [col for col in per_sample.columns if col not in front_cols]
        ]
        per_sample_path = args.results_dir / "per_sample_sensitivity.csv"
        per_sample.to_csv(per_sample_path, index=False)

        # figS_task1C panel c input: Ding rows only, sensitivity columns only.
        ding_cols = ["dataset", "sample", "tool", "sensitivity_gain", "sensitivity_loss"]
        ding_sens = per_sample.loc[
            per_sample["dataset"].eq("Ding_nature_wes"), ding_cols
        ].sort_values(["tool", "sample"], kind="mergesort")
        ding_sens.to_csv(
            args.results_dir / "ding_gain_loss_per_sample_sensitivity.csv", index=False
        )

        excluded_gain = int(per_sample["sensitivity_gain"].isna().sum())
        excluded_loss = int(per_sample["sensitivity_loss"].isna().sum())
        print(
            "\n[paired-test exclusions] "
            f"undefined gain sensitivity rows={excluded_gain}; "
            f"undefined loss sensitivity rows={excluded_loss}"
        )

        tool_rows = [
            run_tool_test(tool_df)
            for _, tool_df in per_sample.groupby("tool", sort=True)
        ]
        per_tool = adjust_pvalues(pd.DataFrame(tool_rows))
        per_tool_path = args.results_dir / "per_tool_wilcoxon.csv"
        per_tool.to_csv(per_tool_path, index=False)

        plot_path = args.results_dir / "gain_vs_loss_sensitivity.png"
        plot_gain_vs_loss(per_sample, plot_path, dpi=CONFIG["plot_dpi"])
        conclusion_path = args.results_dir / "gain_vs_loss_sensitivity_conclusion.png"
        plot_conclusion_summary(per_tool, conclusion_path, dpi=CONFIG["plot_dpi"])
        ding_tool_plot_path = args.results_dir / "ding_gain_loss_tool_boxplot.png"
        ding_tool_stats_path = args.results_dir / "ding_gain_loss_tool_boxplot_stats.csv"
        plot_ding_tool_boxplot(
            per_sample, ding_tool_plot_path, ding_tool_stats_path, dpi=CONFIG["plot_dpi"]
        )
        ding_overall_plot_path = args.results_dir / "ding_gain_loss_overall_boxplot.png"
        ding_overall_stats_path = args.results_dir / "ding_gain_loss_overall_boxplot_stats.csv"
        plot_ding_overall_boxplot(
            per_sample,
            ding_overall_plot_path,
            ding_overall_stats_path,
            dpi=CONFIG["plot_dpi"],
        )

        print(f"\n[done] wrote {per_sample_path}")
        print(f"[done] wrote {per_tool_path}")
        print(f"[done] wrote {plot_path}")
        print(f"[done] wrote {conclusion_path}")
        print(f"[done] wrote {ding_tool_plot_path}")
        print(f"[done] wrote {ding_tool_stats_path}")
        print(f"[done] wrote {ding_overall_plot_path}")
        print(f"[done] wrote {ding_overall_stats_path}")
        print_summary(
            per_tool,
            alpha=CONFIG["fdr_alpha"],
            borderline=CONFIG["fdr_borderline"],
        )
    except Exception as exc:
        print(f"\n[stop] {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
