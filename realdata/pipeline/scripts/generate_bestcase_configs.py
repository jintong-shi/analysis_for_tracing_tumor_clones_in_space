#!/usr/bin/env python3
"""
Render one Snakemake config per best-case run.

Inputs (both ship with the repo, no external data needed):
  config/samples.tsv         per-sample paths and deconvolution keys
  config/bestcase_spec.tsv   which (sample, combo, inferCNV mode) to run,
                             produced by derive_bestcase_spec.py

Output:
  config/<dataset>/<sample>__<run_combo>[__infercnv_hmm].yaml

One config == one reference combo, so `run_combos` always holds a single entry.
The `__infercnv_hmm` variant reruns inferCNV only, with HMM on; the shell driver
renames the existing non-HMM output folder to `infercnv_regular` first, which is
how `infercnv_expr` and `infercnv_cnv` end up in separate folders.

All paths use the placeholders documented in docs/PATHS.md and must be
substituted before running.

Usage: generate_bestcase_configs.py [--clean]
"""

import argparse
import csv
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.abspath(os.path.join(HERE, "..", "config"))

REF = {
    "gene_pos": "<REF_ROOT>/gencode_v38_gene_pos.txt",
    "chrom_arm": "<TOOLS_ROOT>/CONICSmat/chromosome_arm_positions_grch38.txt",
    "gene_map": "<REF_ROOT>/hg38_gene_map.csv",
    "cell_cycle": "<TOOLS_ROOT>/copyVAE/data/Macosko_cell_cycle_genes.txt",
    "starch_repo": "<TOOLS_ROOT>",
    "numbat_gmap": "<TOOLS_ROOT>/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz",
    "numbat_eagle": "<TOOLS_ROOT>/Eagle_v2.4.1/eagle",
    "numbat_vcf": "<REF_ROOT>/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf",
    "numbat_panel": "<REF_ROOT>/numbat/1000G_hg38",
    "xclone_genes": "<REF_ROOT>/xclone/annotate_genes_hg38_update_20230126.txt",
}
DISPLAY = {
    "allspots_no_ref": "AllSpots_no-ref",
    "decon_no_ref": "Decon_no-ref",
    "decon_stromaimm_ref": "Decon_StromaImm-ref",
    "histo_tumor_benign_ref": "HistoTumor_Benign-ref",
}
DECON_COMBOS = {"decon_no_ref", "decon_stromaimm_ref"}
HISTO_COMBOS = {"histo_tumor_benign_ref"}


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def yaml_list(items):
    return "[" + ", ".join('"%s"' % i for i in items) + "]"


def render(s, spec):
    combo = spec["run_combo"]
    tools = spec["tools"].split(",")
    hmm = spec["infercnv_hmm"] == "true"
    L = ["## sample input ##",
         'sample_type: "%s"' % s["sample_type"],
         'sample: "%s"' % s["sample"],
         'input_path: "%s"' % s["input_path"]]
    if combo in DECON_COMBOS and s["c2l_run"] != "True" and s["decon_csv"]:
        L.append('decon_csv: "%s"' % s["decon_csv"])
    if combo in HISTO_COMBOS and s["histology_csv"]:
        # prostate_eva_ffpe has no per-sample histology CSV; the pipeline exports
        # it from the Seurat object (Snakefile warns rather than fails).
        L.append('histology_csv: "%s"' % s["histology_csv"])
    L += ["threads: 8", "",
          "tools: " + yaml_list(tools),
          "run_combos: " + yaml_list([combo]), ""]

    if combo in DECON_COMBOS:
        L += ["decon:",
              '  epi_celltypes: "%s"' % s["epi_celltypes"],
              '  stroma_celltypes: "%s"' % s["stroma_celltypes"],
              '  imm_celltypes: "%s"' % s["imm_celltypes"],
              "",
              "cell2location:",
              "  run: %s" % s["c2l_run"]]
        if s["sc_ref_folder"]:
            L.append('  sc_ref_folder: "%s"' % s["sc_ref_folder"])
        if s["c2l_run"] == "True":
            L += ['  cell_type_col: "%s"' % s["c2l_cell_type_col"],
                  '  batch_key: "%s"' % s["c2l_batch_key"],
                  "  gpu: False",
                  "  N_cells_per_location: 10",
                  "  detection_alpha: 20"]
        L += ["", "params:", "  epi_threshold: 0.5", "  ref_top_frac: 0.10", ""]

    if combo in HISTO_COMBOS and s["tumor_labels"]:
        L += ["histology:",
              "  tumor_labels: " + yaml_list(s["tumor_labels"].split(",")),
              "  benign_labels: " + yaml_list(s["benign_labels"].split(",")),
              "  epi_labels: " + yaml_list(s["epi_labels"].split(",")), ""]

    L += ["display:", '  %s: "%s"' % (combo, DISPLAY[combo]), "",
          "## output path ##",
          'outroot: "%s"' % s["outroot"], "",
          "## ref path ##"]

    if "infercnv" in tools:
        L += ["infercnv:",
              '  gene_pos_file: "%s"' % REF["gene_pos"],
              "  HMM: %s" % ("true" if hmm else "false")]
    if "conicsmat" in tools:
        L += ["conicsmat:", '  chrom_arm_file: "%s"' % REF["chrom_arm"]]
    if "copyvae" in tools:
        L += ['copyvae_gene_map_csv: "%s"' % REF["gene_map"],
              'copyvae_cell_cycle_file: "%s"' % REF["cell_cycle"],
              "copyvae:",
              '  cell_cycle_file: "%s"' % REF["cell_cycle"],
              '  gene_map_csv: "%s"' % REF["gene_map"]]
    if "starch" in tools:
        L += ["starch:", "  n_clusters: 3", "  output_name: STARCH",
              '  repo_dir: "%s"' % REF["starch_repo"]]
    if "numbat" in tools:
        L += ["numbat:",
              '  bam_path: "%s"' % s["bam_path"],
              '  gmap_file: "%s"' % REF["numbat_gmap"],
              '  eagle_bin: "%s"' % REF["numbat_eagle"],
              '  snp_vcf: "%s"' % REF["numbat_vcf"],
              '  panel_dir: "%s"' % REF["numbat_panel"],
              "  UMI_tag: UB", "  cell_tag: CB",
              "  max_entropy: 0.8", "  min_LLR: 5", "  plot: true"]
    if "xclone" in tools:
        L += ["xclone:",
              '  gene_region: "%s"' % REF["xclone_genes"],
              '  genome_mode: "hg38_genes"']
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean", action="store_true",
                    help="remove previously generated per-dataset folders first")
    args = ap.parse_args()

    samples = {r["sample"]: r for r in read_tsv(os.path.join(CONFIG_DIR, "samples.tsv"))}
    spec = read_tsv(os.path.join(CONFIG_DIR, "bestcase_spec.tsv"))

    datasets = sorted({r["dataset"] for r in spec})
    if args.clean:
        for d in datasets:
            shutil.rmtree(os.path.join(CONFIG_DIR, d), ignore_errors=True)

    missing_bam = []
    n = 0
    for row in spec:
        s = samples[row["sample"]]
        if "numbat" in row["tools"].split(",") and not s["bam_path"]:
            missing_bam.append("%s/%s %s" % (row["dataset"], row["sample"], row["run_combo"]))
        name = "%s__%s%s.yaml" % (row["sample"], row["run_combo"],
                                  "__infercnv_hmm" if row["infercnv_hmm"] == "true" else "")
        out_dir = os.path.join(CONFIG_DIR, row["dataset"])
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, name), "w") as fh:
            fh.write(render(s, row))
        n += 1

    print("wrote %d configs under %s" % (n, os.path.normpath(CONFIG_DIR)))
    for d in datasets:
        print("  %-28s %d" % (d, len(os.listdir(os.path.join(CONFIG_DIR, d)))))
    if missing_bam:
        print("\nWARNING: numbat requested but samples.tsv has no bam_path:")
        for m in missing_bam:
            print("  -", m)


if __name__ == "__main__":
    main()
