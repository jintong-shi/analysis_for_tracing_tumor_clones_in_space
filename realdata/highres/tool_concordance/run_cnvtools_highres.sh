#!/usr/bin/env bash
# run_cnvtools_highres.sh -- CONICSmat + CopyKAT on the high-res ST samples.
#
# Usage:
#   bash realdata/highres/tool_concordance/run_cnvtools_highres.sh <platform> [tool]
#     platform : atera | novaseq | visiumhd
#     tool     : both (default) | conicsmat | copykat
#   THREADS=8 overridable via the environment.
#
# Runs the shared pipeline scripts realdata/pipeline/R/{run_conicsmat.R,run_copykat.R}
# UNMODIFIED, with the same parameters the Snakefile rules use, so the high-res
# samples are directly comparable to the Visium benchmark:
#   threads 8            (Snakefile GLOBAL_THREADS, config "threads: 8")
#   exp_thresh 0.1       (Snakefile conicsmat.exp_thresh default)
#   windowsize 101       (Snakefile GLOBAL_WINDOW_SIZE default)
#   threshold 4          (Snakefile conicsmat.clip_thresh default)
#   gene_pos  gencode_v38_gene_pos.txt
#   chrom_arm chromosome_arm_positions_grch38.txt
# CopyKAT takes no extra parameters; its seed (1234) is set inside run_copykat.R.
#
# Prerequisite: the matching prep script has been run, so <platform>/input/
# holds matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz + keep.txt + ref.txt:
#   realdata/highres/atera/atera_20_cnvtools_prep.R
#   realdata/highres/novaseq/novaseq_20_g3_cnvtools_prep.R
#   realdata/highres/visiumHD/infercnv/visiumhd_20_cnvtools_prep.R
#
# Numbat is deliberately absent: the NovaSeq X download ships no BAM/CRAM, and
# numbat needs a pileup over aligned reads for allele counts.
#
# Both tools cluster the tumour cells with a dense cell x cell distance matrix.
# At atera's 84,247 and novaseq's 90,000 tumour cells that is ~57 / ~65 GB for
# the distance matrix alone, before CONICSmat's dense normMat() copy. Run
# visiumhd (2,158 tumour spots) first to confirm the wiring, and watch RAM on
# the other two -- check `free -g` before launching.
# ------------------------------------------------------------------
set -uo pipefail

PROJ="<DATA_ROOT>"
HIGHRES=$PROJ/high_res_ST
RSCRIPT="<CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript"
RDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../pipeline/R" && pwd)

GENE_POS="<REF_ROOT>/gencode_v38_gene_pos.txt"
CHROM_ARM="<TOOLS_ROOT>/CONICSmat/chromosome_arm_positions_grch38.txt"

THREADS=${THREADS:-8}
EXP_THRESH=0.1
WINDOWSIZE=101
THRESHOLD=4

platform=${1:-}
tool=${2:-both}

case "$platform" in
  atera)    base=$HIGHRES/atera/processed_data/cnvtools/atera
            logs=$HIGHRES/atera/logs ;;
  novaseq)  base=$HIGHRES/novaseq/processed_data/cnvtools/IDC_Grade3_ds90k
            logs=$HIGHRES/novaseq/logs ;;
  visiumhd) base=$HIGHRES/visiumHD/processed_data/cnvtools/pseudo_r40_Smax20
            logs=$HIGHRES/visiumHD/logs ;;
  # 60k subsamples: both tools call hclust(), which caps out at 65,536 cells.
  atera_ds60k)   base=$HIGHRES/atera/processed_data/cnvtools/atera_ds60k
                 logs=$HIGHRES/atera/logs ;;
  novaseq_ds60k) base=$HIGHRES/novaseq/processed_data/cnvtools/IDC_Grade3_ds90k_ds60k
                 logs=$HIGHRES/novaseq/logs ;;
  atera_ds20k)   base=$HIGHRES/atera/processed_data/cnvtools/atera_ds20k
                 logs=$HIGHRES/atera/logs ;;
  novaseq_ds20k) base=$HIGHRES/novaseq/processed_data/cnvtools/IDC_Grade3_ds90k_ds20k
                 logs=$HIGHRES/novaseq/logs ;;
  *) echo "usage: $0 <atera|novaseq|visiumhd|{atera,novaseq}_ds{20,60}k> [both|conicsmat|copykat]" >&2; exit 2 ;;
esac

IN=$base/input
for f in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz keep.txt ref.txt; do
  [ -s "$IN/$f" ] || { echo "[$platform] missing $IN/$f -- run the prep script first" >&2; exit 1; }
done
[ -s "$CHROM_ARM" ] || { echo "missing chrom arm file: $CHROM_ARM" >&2; exit 1; }
[ -s "$GENE_POS" ]  || { echo "missing gene position file: $GENE_POS" >&2; exit 1; }

mkdir -p "$logs"
echo "[$platform] input=$IN  keep=$(wc -l < "$IN/keep.txt")  ref=$(wc -l < "$IN/ref.txt")  threads=$THREADS"

rc_copykat=skipped
rc_conics=skipped

if [ "$tool" = both ] || [ "$tool" = copykat ]; then
  out=$base/COPYKAT; mkdir -p "$out"
  echo "[$platform] CopyKAT -> $out"
  "$RSCRIPT" --vanilla "$RDIR/run_copykat.R" \
    --visium_folder      "$IN" \
    --output_dir         "$out" \
    --keep_barcodes_file "$IN/keep.txt" \
    --ref_barcodes_file  "$IN/ref.txt" \
    --threads            "$THREADS" \
    > "$logs/cnvtools_${platform}_copykat.log" 2>&1
  rc_copykat=$?
  echo "[$platform] CopyKAT exit=$rc_copykat (log: $logs/cnvtools_${platform}_copykat.log)"
fi

if [ "$tool" = both ] || [ "$tool" = conicsmat ]; then
  out=$base/CONICSmat; mkdir -p "$out"
  echo "[$platform] CONICSmat -> $out"
  "$RSCRIPT" --vanilla "$RDIR/run_conicsmat.R" \
    --visium_folder      "$IN" \
    --output_dir         "$out" \
    --keep_barcodes_file "$IN/keep.txt" \
    --ref_barcodes_file  "$IN/ref.txt" \
    --chrom_arm_file     "$CHROM_ARM" \
    --gene_pos_file      "$GENE_POS" \
    --exp_thresh         "$EXP_THRESH" \
    --windowsize         "$WINDOWSIZE" \
    --threshold          "$THRESHOLD" \
    --threads            "$THREADS" \
    > "$logs/cnvtools_${platform}_conicsmat.log" 2>&1
  rc_conics=$?
  echo "[$platform] CONICSmat exit=$rc_conics (log: $logs/cnvtools_${platform}_conicsmat.log)"
fi

echo "[$platform] done. copykat=$rc_copykat conicsmat=$rc_conics"
echo "[$platform] success markers:"
ls -1 "$base"/COPYKAT/COPYKAT.OK "$base"/CONICSmat/CONICSmat.OK 2>/dev/null || echo "  (none yet)"
