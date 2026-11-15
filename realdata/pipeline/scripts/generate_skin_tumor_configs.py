#!/usr/bin/env python3
"""
Generate per-sample Snakemake configs for the skin_tumor (cSCC) Visium samples
from Patient3 / Patient4 / Patient6 (Kim et al.).

For each sample we emit one YAML at:
    config/generated_skin_tumor/config_<sample>_run.yaml
plus a list file:
    config/skin_tumor_configs.list
that the runner script consumes.

Per-patient single-cell references (built earlier by
public_data/skin_tumor/scRNA_GEO/prepare_cell2location_ref.py) are picked up
from cell2location_ref/{P3,P4,P6}/.

Tools requested:
    infercnv, conicsmat, copykat, scevan, numbat, xclone

If the sample's spaceranger folder has no BAM, numbat and xclone are dropped
(both need possorted_genome_bam.bam). conicsmat / xclone / numbat only build
targets on combos with a reference — i.e. decon_stromaimm_ref — that gating
lives in the Snakefile, so we just leave them in the tools list.

Run modes (run_combos):
    allspots_no_ref, decon_no_ref, decon_stromaimm_ref
"""

import os
import sys

PIPE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA_ROOT = "<DATA_ROOT>/public_data/skin_tumor"
SC_REF_ROOT = f"{DATA_ROOT}/scRNA_GEO/cell2location_ref"
OUTROOT = f"{DATA_ROOT}/run_CNV"
CONFIG_DIR = f"{PIPE_DIR}/config/generated_skin_tumor"
LIST_FILE = f"{PIPE_DIR}/config/skin_tumor_configs.list"

# (sample_id, patient, visium input_path)
# Patient3 sections come straight from a folder (no outs/, no BAM).
# Patient4/6 sections use spaceranger_rerun/outs (with BAM).
SAMPLES = [
    ("skin_P3_B1", "P3",
     f"{DATA_ROOT}/visium_from_kim/Patient3/V10F24-014_B1"),
    ("skin_P3_C1", "P3",
     f"{DATA_ROOT}/visium_from_kim/Patient3/V10F24-014_C1"),
    ("skin_P4_V1", "P4",
     f"{DATA_ROOT}/visium_from_kim/Patient4/V19T26-101_V1/spaceranger_rerun/outs"),
    ("skin_P4_V2", "P4",
     f"{DATA_ROOT}/visium_from_kim/Patient4/V19T26-101_V2/spaceranger_rerun/outs"),
    ("skin_P6_V3", "P6",
     f"{DATA_ROOT}/visium_from_kim/Patient6/V19T26-101_V3/spaceranger_rerun/outs"),
    ("skin_P6_V4", "P6",
     f"{DATA_ROOT}/visium_from_kim/Patient6/V19T26-101_V4/spaceranger_rerun/outs"),
]

# Cell type buckets for the Decon_StromaImm-ref selector. These are
# values of `level2_celltype` in the GSE144236 metadata. The R selector does
# case-insensitive substring matching against Decon.csv columns, so partial
# names with spaces are fine.
EPI_CT_BY_PATIENT = {
    "P3": "Keratinocyte",
    "P4": "Keratinocyte",
    "P6": "Tumor_KC_Basal,Tumor_KC_Cyc,Tumor_KC_Diff,TSK,Normal_KC_Basal",
}
STROMA_CT  = "Endothelial Cell,Fibroblast"
IMM_CT     = "Tcell,B Cell,NK,Mac,LC,CD1C,CLEC9A,MDSC,ASDC,PDC"

# cell2location config — uses level2_celltype as label, tum_norm as batch.
C2L_CT_COL    = "level2_celltype"
C2L_BATCH_KEY = "tum_norm"
C2L_PYBIN     = "<CONDA_ROOT>/envs/cell2loc_env/bin/python"

# Numbat reference resources (mirrors existing prostate_ff config).
NUMBAT_GMAP   = "<TOOLS_ROOT>/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
NUMBAT_EAGLE  = "<TOOLS_ROOT>/Eagle_v2.4.1/eagle"
NUMBAT_SNPVCF = "<REF_ROOT>/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"
NUMBAT_PANEL  = "<REF_ROOT>/numbat/1000G_hg38"


def render_config(sample, patient, input_path, has_bam, bam_path):
    base_tools = ["infercnv", "conicsmat", "copykat", "scevan"]
    if has_bam:
        tools = base_tools + ["numbat", "xclone"]
    else:
        tools = base_tools
    tools_str = "[" + ", ".join(f'"{t}"' for t in tools) + "]"

    sc_ref_folder = f"{SC_REF_ROOT}/{patient}"

    numbat_block = ""
    if has_bam:
        numbat_block = f"""
numbat:
  bam_path: "{bam_path}"
  gmap_file: "{NUMBAT_GMAP}"
  eagle_bin: "{NUMBAT_EAGLE}"
  snp_vcf:   "{NUMBAT_SNPVCF}"
  panel_dir: "{NUMBAT_PANEL}"
  UMI_tag: UB
  cell_tag: CB
  max_entropy: 0.8
  min_LLR: 5
  plot: true
"""

    xclone_block = """
xclone:
  gene_region: "<REF_ROOT>/xclone/annotate_genes_hg38_update_20230126.txt"
  genome_mode: "hg38_genes"
"""
    if not has_bam:
        xclone_block = ""

    return f"""## sample input ##
sample_type: "freshfrozen"
sample: "{sample}"
input_path: "{input_path}"
threads: 8

decon:
  epi_celltypes: "{EPI_CT_BY_PATIENT[patient]}"
  stroma_celltypes: "{STROMA_CT}"
  imm_celltypes: "{IMM_CT}"

cell2location:
  run: True
  sc_ref_folder: "{sc_ref_folder}"
  cell_type_col: "{C2L_CT_COL}"
  batch_key: "{C2L_BATCH_KEY}"
  python_bin: "{C2L_PYBIN}"
  gpu: False
  N_cells_per_location: 10
  detection_alpha: 20

params:
  epi_threshold: 0.5
  ref_top_frac: 0.10

run_combos: ["allspots_no_ref", "decon_no_ref", "decon_stromaimm_ref"]
tools: {tools_str}

display:
  allspots_no_ref: "AllSpots_no-ref"
  decon_no_ref: "Decon_no-ref"
  decon_stromaimm_ref: "Decon_StromaImm-ref"

## output path ##
outroot: "{OUTROOT}"

## ref path ##
infercnv:
  gene_pos_file: "<REF_ROOT>/gencode_v38_gene_pos.txt"
conicsmat:
  chrom_arm_file: "<TOOLS_ROOT>/CONICSmat/chromosome_arm_positions_grch38.txt"
{xclone_block}
refs:
  genome_build: hg38
  eagle_bin: {NUMBAT_EAGLE}
  eagle_dir: <TOOLS_ROOT>/Eagle_v2.4.1
  gmap_file: {NUMBAT_GMAP}
  snp_vcf:   <REF_ROOT>/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
  panel_dir: <REF_ROOT>/numbat/1000G_hg38
{numbat_block}"""


def main():
    os.makedirs(CONFIG_DIR, exist_ok=True)
    os.makedirs(OUTROOT, exist_ok=True)

    list_lines = []
    preflight_failed = False

    for sample, patient, input_path in SAMPLES:
        if not os.path.isdir(input_path):
            print(f"[FATAL] missing input_path: {input_path}", file=sys.stderr)
            preflight_failed = True
            continue

        sc_ref = f"{SC_REF_ROOT}/{patient}"
        for need in ("counts.mtx", "metadata.csv", "gene_name.csv"):
            if not os.path.exists(os.path.join(sc_ref, need)):
                print(f"[FATAL] missing {need} in {sc_ref}", file=sys.stderr)
                preflight_failed = True

        bam_path = os.path.join(input_path, "possorted_genome_bam.bam")
        has_bam = os.path.exists(bam_path)

        cfg_path = os.path.join(CONFIG_DIR, f"config_{sample}_run.yaml")
        with open(cfg_path, "w") as f:
            f.write(render_config(sample, patient, input_path, has_bam, bam_path))
        list_lines.append(cfg_path)
        bam_msg = "BAM=yes (numbat+xclone enabled)" if has_bam else "BAM=no (numbat+xclone skipped)"
        print(f"  wrote {cfg_path}  [{bam_msg}]")

    with open(LIST_FILE, "w") as f:
        f.write("\n".join(list_lines) + "\n")
    print(f"\nWrote list file: {LIST_FILE} ({len(list_lines)} configs)")

    if preflight_failed:
        sys.exit("Preflight check failed; fix the issues above before running snakemake.")


if __name__ == "__main__":
    main()
