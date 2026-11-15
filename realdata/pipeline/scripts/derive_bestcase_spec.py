#!/usr/bin/env python3
"""
Derive the best-case run spec: which (sample, reference combo, inferCNV mode)
the pipeline must produce so that every sample x tool record shown in the final
figures can be rebuilt -- and nothing else.

"Best case" is a deterministic rule, not a metric maximum. The rules are the
ones applied in memory by figures/utils/final_record_select.R:

  task1  final_record_select_task1()  sel_pri priority:
           histology available -> HistoTumor_Benign-ref
           otherwise           -> Decon_StromaImm-ref
           (fallbacks AllSpots_no-ref, Decon_no-ref are only reached when the
           preferred run is absent; this script always emits the preferred one)
  task2  final_record_select_task2()
           STARCH / CalicoST / TUSCAN -> AllSpots_no-ref
           otherwise same as task1
  task3  final_record_select_task3()
           numbat -> Decon_StromaImm-ref ; every other tool -> AllSpots_no-ref

The tool set per sample is read from the three frozen metric tables rather than
assumed, so a tool only appears where it actually produced a figure record.

inferCNV is split into two runs because they land in different output folders
and feed different figure tool ids (see 01_standardize_cnv_output_realdata.R):
  infercnv_expr <- <combo>/infercnv_regular   (HMM off)
  infercnv_cnv  <- <combo>/infercnv           (HMM on, needs pred_cnv_genes.dat)
The shell drivers create that split by renaming the non-HMM output folder to
`infercnv_regular` before re-running with HMM on.

Usage:
  derive_bestcase_spec.py --t1 <pseudobulk_vs_wes_metrics_all_datasets.csv> \
                          --t2 <subclone_spatial_assortativity_summary.csv> \
                          --t3 <metrics_combined_relabelled.csv> \
                          [--overview ../../../figures/utils/sample_overview.csv] \
                          [--out ../config/bestcase_spec.tsv]
"""

import argparse
import csv
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OVERVIEW = os.path.join(HERE, "..", "..", "..", "figures", "utils", "sample_overview.csv")
DEFAULT_OUT = os.path.abspath(os.path.join(HERE, "..", "config", "bestcase_spec.tsv"))

# figure tool id -> (pipeline tool name, needs inferCNV HMM)
TOOL_MAP = {
    "infercnv_expr": ("infercnv", False),
    "infercnv_cnv": ("infercnv", True),
    "InferCNV": ("infercnv", False),
    "numbat_expr": ("numbat", False),
    "numbat_cnv": ("numbat", False),
    "Numbat": ("numbat", False),
    "numbat": ("numbat", False),
    "COPYKAT": ("copykat", False),
    "copykat": ("copykat", False),
    "SCEVAN": ("scevan", False),
    "scevan": ("scevan", False),
    "CONICSmat": ("conicsmat", False),
    "copyVAE": ("copyvae", False),
    "copyvae": ("copyvae", False),
    "infercnvpy": ("infercnvpy", False),
    "STARCH": ("starch", False),
    "TUSCAN": ("tuscan", False),
    "tuscan": ("tuscan", False),
    "CalicoST": ("calicost", False),
    "XCLONE": ("xclone", False),
}
REF_FREE_TOOLS = {"STARCH", "CalicoST", "TUSCAN"}
# FFPE samples never run these (Snakefile `ffpe_tools` gate drops them silently)
NOT_ON_FFPE = {"numbat", "xclone", "calicost"}


def read_rows(path):
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t1", required=True)
    ap.add_argument("--t2", required=True)
    ap.add_argument("--t3", required=True)
    ap.add_argument("--overview", default=DEFAULT_OVERVIEW)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    overview = {}
    for r in read_rows(args.overview):
        s = r["sample"].strip()
        if s:
            overview[s] = (
                r["datasets"].strip(),
                r["histology_available"].strip().lower() == "yes",
                r["sample_type"].strip().lower(),
            )

    need = defaultdict(set)
    unmapped = set()

    def add(sample, fig_tool, combo):
        if fig_tool not in TOOL_MAP:
            unmapped.add(fig_tool)
            return
        tool, hmm = TOOL_MAP[fig_tool]
        need[(overview[sample][0], sample, combo, hmm)].add(tool)

    def preferred(sample):
        return "histo_tumor_benign_ref" if overview[sample][1] else "decon_stromaimm_ref"

    for r in read_rows(args.t1):
        s = r.get("sample", "").strip()
        if s in overview and r.get("status", "ok") == "ok":
            add(s, r["tool"].strip(), preferred(s))

    for r in read_rows(args.t2):
        s = r.get("sample", "").strip()
        if s in overview:
            t = r["tool"].strip()
            add(s, t, "allspots_no_ref" if t in REF_FREE_TOOLS else preferred(s))

    for r in read_rows(args.t3):
        s = r.get("sample", "").strip()
        if s in overview:
            t = r["tool"].strip()
            add(s, t, "decon_stromaimm_ref" if t.lower() == "numbat" else "allspots_no_ref")

    if unmapped:
        raise SystemExit("Unmapped figure tool ids (extend TOOL_MAP): %s" % sorted(unmapped))

    rows = []
    dropped = []
    for (dataset, sample, combo, hmm), tools in sorted(need.items(), key=str):
        is_ffpe = overview[sample][2] == "ffpe"
        keep = sorted(t for t in tools if not (is_ffpe and t in NOT_ON_FFPE))
        for t in sorted(tools):
            if is_ffpe and t in NOT_ON_FFPE:
                dropped.append("%s/%s %s: %s is not run on FFPE" % (dataset, sample, combo, t))
        if not keep:
            continue
        rows.append(dict(dataset=dataset, sample=sample, run_combo=combo,
                         infercnv_hmm="true" if hmm else "false",
                         tools=",".join(keep)))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["dataset", "sample", "run_combo", "infercnv_hmm", "tools"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print("wrote %d spec rows -> %s" % (len(rows), args.out))
    print("samples covered: %d" % len({r["sample"] for r in rows}))
    if dropped:
        print("\nFFPE gate dropped %d tool entries:" % len(dropped))
        for d in dropped:
            print("  -", d)


if __name__ == "__main__":
    main()
