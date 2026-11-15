import os
import sys
import pandas as pd


METADATA_COLS = {
    "barcode", "barcodes", "spot", "cell", "cell_id", "cellid",
    "in_tissue", "array_row", "array_col",
    "pxl_row_in_fullres", "pxl_col_in_fullres",
    "sample", "_indices", "_scvi_batch", "_scvi_labels",
}


def _parse_keys(s):
    return [x.strip() for x in str(s or "").split(",") if x.strip()]


def _resolve_paths_and_params():
    if "snakemake" in globals():
        input_csv = str(snakemake.input[0])
        output_csv = str(snakemake.output[0])
        params = getattr(snakemake, "params", {})
        epi = params.get("epi_celltypes", "") if hasattr(params, "get") else ""
        stroma = params.get("stroma_celltypes", "") if hasattr(params, "get") else ""
        imm = params.get("imm_celltypes", "") if hasattr(params, "get") else ""
    else:
        if len(sys.argv) < 3:
            print("Usage: python calc_proportions.py <input_csv> <output_csv>")
            sys.exit(1)
        input_csv = sys.argv[1]
        output_csv = sys.argv[2]
        epi = stroma = imm = ""
    requested = _parse_keys(epi) + _parse_keys(stroma) + _parse_keys(imm)
    return input_csv, output_csv, requested


def _select_score_columns(df, requested_keys):
    # Keep numeric-like columns only.
    num_df = df.apply(pd.to_numeric, errors="coerce")
    candidate = [c for c in num_df.columns if not num_df[c].isna().all()]

    # Remove known metadata columns.
    score_cols = [c for c in candidate if c.lower() not in METADATA_COLS]

    # If keys are provided and any key matches, require key-matched columns to avoid noisy columns.
    if requested_keys:
        keys_l = [k.lower() for k in requested_keys]
        matched = [c for c in score_cols if any(k in c.lower() for k in keys_l)]
        if matched:
            # Build denominator from all detected score columns to satisfy score/sum(all scores),
            # but keep matched columns for diagnostics.
            print(f"[calc_proportions] matched key columns ({len(matched)}): {matched}")
        else:
            print("[calc_proportions] warning: no columns matched decon keys; using all inferred score columns.")

    if not score_cols:
        raise ValueError("No cell2location score columns found after filtering metadata columns.")

    return num_df[score_cols]


def main():
    input_csv, output_csv, requested_keys = _resolve_paths_and_params()
    print(f"[calc_proportions] Input:  {input_csv}")
    print(f"[calc_proportions] Output: {output_csv}")

    df = pd.read_csv(input_csv, index_col=0)
    score_df = _select_score_columns(df, requested_keys)

    # Convert score -> proportion per spot: score / sum(all scores in the spot).
    row_sums = score_df.sum(axis=1)
    prop_df = score_df.div(row_sums.replace(0, pd.NA), axis=0).fillna(0)

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    if prop_df.index.name is None:
        prop_df.index.name = "Barcode"
    prop_df.to_csv(output_csv)

    print(f"[calc_proportions] saved {prop_df.shape[0]} spots x {prop_df.shape[1]} cell-type columns.")


if __name__ == "__main__":
    main()
