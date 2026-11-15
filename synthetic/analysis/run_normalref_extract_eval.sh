#!/usr/bin/env bash
# Extract + normalize (01) then evaluate (02) the normal-reference synthetic batch,
# producing benchmark_metrics_normalref.csv.
# Every Rep_* under $FLAT must contain synthetic/, Standardized_Inputs/, Tool_Outputs/
# and event_log.tsv (missing event_log.tsv leaves narrow/broad rows NA), and CopyVAE
# reps must already have copy.npy from run_copyvae.py.
# Re-runnable: 02 resumes from an existing --out_file.
# Usage:  nohup bash run_normalref_extract_eval.sh > normalref_extract_eval.out 2>&1 &
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"    # this analysis/ folder

# ---- inputs / outputs --------------------------------------------------------
SYN_RUNS="<DATA_ROOT>/synthetic/numbat_dirichlet_fix_test/runs"
FLAT="${FLAT:-$SYN_RUNS/complete_batch_normalref}"             # flat batch of Rep_* sample dirs
RESULTS="${RESULTS:-$FLAT/results}"
OUT02="${OUT02:-$RESULTS/benchmark_metrics_normalref.csv}"
EXTRACT_R="$SCRIPT_DIR/extract_and_normalize.R"
EVAL_R="$SCRIPT_DIR/evaluate_results.R"
COPYVAE_PY="$SCRIPT_DIR/postprocess_copyvae.py"

LOG_DIR="${LOG_DIR:-$RESULTS/extract_eval_logs}"; mkdir -p "$LOG_DIR" "$RESULTS"
TS="$(date +%Y%m%d_%H%M%S)"
L01A="$LOG_DIR/01a_extract_refbased_$TS.log"
L01B="$LOG_DIR/01b_extract_reffree_$TS.log"
L01C="$LOG_DIR/01c_copyvae_postprocess_$TS.log"
L02="$LOG_DIR/02_eval_$TS.log"

# ---- tool lists --------------------------------------------------------------
TOOLS01_REFBASED="InferCNV,CopyKAT,InferCNVPy,SCEVAN,CONICSmat,XClone,Numbat"
TOOLS01_REFFREE="STARCH,CalicoST"
# 02 takes the Tool_Outputs subdir names as produced above, plus CopyVAE.
TOOLS02="InferCNV_HMM17_TumorFreeNormal,CopyKAT_TumorFreeNormal,InferCNVPy_TumorFreeNormal,SCEVAN_TumorFreeNormal,CONICSmat_TumorFreeNormal,XClone_TumorFreeNormal,Numbat_TumorFreeNormal,STARCH,CalicoST,CopyVAE"

# ---- conda envs --------------------------------------------------------------
CONDA_SH="${CONDA_SH:-<CONDA_ROOT>/etc/profile.d/conda.sh}"
ENV_EXTRACT="${ENV_EXTRACT:-<CONDA_ROOT>/envs/numbat-env}"     # 01a/01b: R extract (EnsDb.Hsapiens.v86)
ENV_COPYVAE="${ENV_COPYVAE:-<CONDA_ROOT>/envs/mysynST_env}"    # 01c: numpy + pandas
ENV_EVAL="${ENV_EVAL:-<CONDA_ROOT>/envs/r4.4}"                 # 02: mclust/pROC/arrow/...
# ENV_* are env PREFIXES (paths), not names: each step calls $ENV/bin/<interpreter>
# explicitly rather than trusting PATH, which a profile pinning another conda can win.
# shellcheck disable=SC1090
source "$CONDA_SH" || { echo "FATAL: cannot source $CONDA_SH" >&2; exit 1; }
act(){   set +u; conda activate "$1" || { set -u; echo "FATAL: conda activate $1 failed" >&2; exit 1; }; set -u; }
deact(){ set +u; conda deactivate || true; set -u; }           # conda scripts are not nounset-safe
for e in "$ENV_EXTRACT" "$ENV_COPYVAE" "$ENV_EVAL"; do
  [ -d "$e" ] || { echo "FATAL: env prefix not a directory: $e" >&2; exit 1; }
done

# ---- check the batch ---------------------------------------------------------
n=0; missing_log=0
for d in "$FLAT"/Rep_*; do
  [ -d "$d" ] || continue
  n=$((n+1))
  [ -f "$d/event_log.tsv" ] || { echo "  WARNING: no event_log.tsv in $(basename "$d") — narrow/broad rows will be NA"; missing_log=$((missing_log+1)); }
done
echo "=== normal-ref extract + eval | host $(hostname) | $(date '+%F %T') ==="
echo "  flat batch : $FLAT  ($n reps, $missing_log without event_log.tsv)"
echo "  out02      : $OUT02"
[ "$n" -gt 0 ] || { echo "FATAL: no Rep_* under $FLAT" >&2; exit 1; }

# ---- 01a/01b. extract + normalize -------------------------------------------
act "$ENV_EXTRACT"
echo "[01a] extract ref-based (_TumorFreeNormal) -> $L01A"
"$ENV_EXTRACT/bin/Rscript" "$EXTRACT_R" \
  --root_dir "$FLAT" \
  --tool_out_name Tool_Outputs \
  --tool_out_suffix _TumorFreeNormal \
  --infercnv_dir_name InferCNV_HMM17 \
  --tools "$TOOLS01_REFBASED" \
  --skip_copyvae > "$L01A" 2>&1
RC01A=$?; echo "  rc_01a: $RC01A"

echo "[01b] extract ref-free (no suffix) -> $L01B"
"$ENV_EXTRACT/bin/Rscript" "$EXTRACT_R" \
  --root_dir "$FLAT" \
  --tool_out_name Tool_Outputs \
  --tool_out_suffix "" \
  --tools "$TOOLS01_REFFREE" \
  --skip_copyvae > "$L01B" 2>&1
RC01B=$?; echo "  rc_01b: $RC01B"
deact

# ---- 01c. CopyVAE post-process ----------------------------------------------
# postprocess_copyvae.py exits 0 (skips) for reps lacking copy.npy.
act "$ENV_COPYVAE"
echo "[01c] CopyVAE post-process -> $L01C"
: > "$L01C"
cv_done=0
for rep in "$FLAT"/Rep_*; do
  [ -d "$rep" ] || continue
  "$ENV_COPYVAE/bin/python" "$COPYVAE_PY" --batch_path "$rep" >> "$L01C" 2>&1
  [ -f "$rep/Tool_Outputs/CopyVAE/processed_log2_spots_x_genes.csv" ] && cv_done=$((cv_done+1))
done
echo "  01c: CopyVAE processed csv present in $cv_done/$n reps (see $L01C)"
deact

if [ "$RC01A" -ne 0 ] || [ "$RC01B" -ne 0 ]; then
  echo "[chain] 01 had a non-zero rc (a=$RC01A b=$RC01B); check logs before trusting 02." >&2
fi
[ "$cv_done" -gt 0 ] || echo "[chain] WARNING: 0 reps have CopyVAE processed csv; CopyVAE will be Skipped in 02." >&2

# ---- 02. evaluate ------------------------------------------------------------
act "$ENV_EVAL"
echo "[02] evaluate (incl CopyVAE) -> $L02"
"$ENV_EVAL/bin/Rscript" "$EVAL_R" \
  --root_dir "$FLAT" \
  --out_file "$OUT02" \
  --tool_out_name Tool_Outputs \
  --std_in_name Standardized_Inputs \
  --tools "$TOOLS02" \
  --metrics "correlation,auc_gain,auc_loss,max_f1" \
  --max_candidates 200 --event_cutoff_mb 3 \
  --workers "${WORKERS:-8}" > "$L02" 2>&1
RC02=$?; echo "  rc_02: $RC02"
deact

echo "=== DONE $(date '+%F %T') ; metrics -> $OUT02 ==="
[ -f "$OUT02" ] && {
  echo "  rows: $(($(wc -l < "$OUT02")-1))"
  echo "  CopyVAE rows: $(grep -c -i copyvae "$OUT02")"
}
