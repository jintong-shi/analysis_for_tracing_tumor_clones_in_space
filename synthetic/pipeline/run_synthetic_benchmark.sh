#!/usr/bin/env bash
# End-to-end driver for the synthetic benchmark: a batch of Rep_* replicates ->
# 10 CNV callers -> benchmark_metrics_normalref.csv.
# Every Rep_* under $BATCH must be a sample dir with synthetic/, baf/ and event_log.tsv.
# Not resumable, except infercnv (skips reps with infercnv.OK) and eval.
# Usage:  nohup bash run_synthetic_benchmark.sh > synthetic_benchmark.out 2>&1 &
#         STEPS=infercnv,eval bash run_synthetic_benchmark.sh   # subset of $ALL_STEPS
set -uo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"    # this pipeline/ folder
ANALYSIS_DIR="$(dirname "$PIPE_DIR")/analysis"
MASTER="$PIPE_DIR/master_pipeline.py"
PREP_TUMORFREE="$PIPE_DIR/prepare_tumorfree_ref_inputs.R"
INFERCNV_HMM17="$PIPE_DIR/run_infercnv_hmm17.R"
EVAL_SH="$ANALYSIS_DIR/run_normalref_extract_eval.sh"

# ---- inputs / outputs --------------------------------------------------------
SYN_RUNS="<DATA_ROOT>/synthetic/numbat_dirichlet_fix_test/runs"
BATCH="${BATCH:-$SYN_RUNS/complete_batch_normalref}"           # flat batch of Rep_* sample dirs
REF_DIR="${REF_DIR:-$SYN_RUNS/tumor_free_run/out/full}"        # the tumour-free reference slide
GENE_ORDER_FILE="${GENE_ORDER_FILE:-<REF_ROOT>/gencode_v38_gene_pos.txt}"
TUMOR_PROP_CUTOFF="${TUMOR_PROP_CUTOFF:-0.2}"                  # a ref spot must be below this
FORCE_PREP="${FORCE_PREP:-FALSE}"                              # TRUE re-writes existing std inputs

TOOL_OUT_NAME="Tool_Outputs"
STD_BASE="Standardized_Inputs"
STD_NORMAL="Standardized_Inputs_TumorFreeNormal"
NORMAL_SUFFIX="_TumorFreeNormal"
INFERCNV_OUT="InferCNV_HMM17${NORMAL_SUFFIX}"

RESULTS="${RESULTS:-$BATCH/results}"
LOG_DIR="${LOG_DIR:-$RESULTS/benchmark_logs}"; mkdir -p "$LOG_DIR" "$RESULTS"
TS="$(date +%Y%m%d_%H%M%S)"
STATUS="$LOG_DIR/run_status_$TS.tsv"

ALL_STEPS="prep,prep_tumorfree,numbat,copykat,scevan,conicsmat,infercnvpy,xclone,infercnv,starch,calicost,copyvae,eval"
STEPS="${STEPS:-$ALL_STEPS}"

want() { [[ ",$STEPS," == *",$1,"* ]]; }

# ---- conda envs --------------------------------------------------------------
# ENV_* are env PREFIXES (paths), not names: master_pipeline.py shells out to bare
# `Rscript` / `python3`, so each step prepends its own bin/ rather than relying on
# `conda activate`, which a profile pinning another conda's bin/ can override.
CONDA_SH="${CONDA_SH:-<CONDA_ROOT>/etc/profile.d/conda.sh}"
ENV_R="${ENV_R:-<CONDA_ROOT>/envs/stcnv-r-shared}"       # arrow, Seurat, infercnv, copykat, SCEVAN, CONICSmat
ENV_NUMBAT="${ENV_NUMBAT:-<CONDA_ROOT>/envs/numbat-env}"
ENV_INFERCNVPY="${ENV_INFERCNVPY:-<CONDA_ROOT>/envs/infercnvpy}"
ENV_XCLONE="${ENV_XCLONE:-<CONDA_ROOT>/envs/xclone}"
ENV_STARCH="${ENV_STARCH:-<CONDA_ROOT>/envs/starch}"
ENV_CALICOST="${ENV_CALICOST:-<CONDA_ROOT>/envs/calicost_env}"
ENV_COPYVAE="${ENV_COPYVAE:-<CONDA_ROOT>/envs/copyvae}"
# shellcheck disable=SC1090
source "$CONDA_SH" || { echo "FATAL: cannot source $CONDA_SH" >&2; exit 1; }

# ---- checks ------------------------------------------------------------------
echo "=== synthetic benchmark | host $(hostname) | $(date '+%F %T') ==="
for f in "$MASTER" "$PREP_TUMORFREE" "$INFERCNV_HMM17" "$EVAL_SH"; do
  [ -f "$f" ] || { echo "FATAL: missing script: $f" >&2; exit 1; }
done
[ -d "$BATCH" ] || { echo "FATAL: batch is not a directory: $BATCH" >&2; exit 1; }

# Check inputs up front, only for the steps actually being run.
want infercnv       && { [ -f "$GENE_ORDER_FILE" ] || { echo "FATAL: no gene order file: $GENE_ORDER_FILE" >&2; exit 1; }; }
want prep_tumorfree && { [ -d "$REF_DIR" ]         || { echo "FATAL: reference slide is not a directory: $REF_DIR" >&2; exit 1; }; }

REP_IDS=""; n=0
for d in "$BATCH"/Rep_*; do
  [ -d "$d" ] || continue
  n=$((n+1))
  id="$(basename "$d" | sed -E 's/^Rep_([0-9]+)_.*/\1/')"
  # A non-matching name (e.g. plain "Rep_3") leaves id as the whole basename.
  [[ "$id" =~ ^[0-9]+$ ]] || { echo "FATAL: cannot parse a rep id from $(basename "$d") — expected Rep_<id>_*" >&2; exit 1; }
  REP_IDS="${REP_IDS:+$REP_IDS,}$id"
  [ -f "$d/event_log.tsv" ] || echo "  WARNING: no event_log.tsv in $(basename "$d") — narrow/broad rows will be NA"
done
[ "$n" -gt 0 ] || { echo "FATAL: no Rep_* under $BATCH" >&2; exit 1; }

# Report every missing env up front.
missing_env=0
for e in "$ENV_R" "$ENV_NUMBAT" "$ENV_INFERCNVPY" "$ENV_XCLONE" "$ENV_STARCH" "$ENV_CALICOST" "$ENV_COPYVAE"; do
  [ -d "$e" ] || { echo "  MISSING ENV: $e" >&2; missing_env=$((missing_env+1)); }
done
[ "$missing_env" -eq 0 ] || { echo "FATAL: $missing_env env prefix(es) missing; build them from realdata/pipeline/envs/" >&2; exit 1; }

echo "  batch    : $BATCH  ($n reps: $REP_IDS)"
echo "  ref slide: $REF_DIR"
echo "  steps    : $STEPS"
echo "  status   : $STATUS"
printf 'timestamp\tstep\trc\n' > "$STATUS"

# ---- step machinery ----------------------------------------------------------
act(){   set +u; conda activate "$1" || { set -u; echo "FATAL: conda activate $1 failed" >&2; exit 1; }; set -u; }
deact(){ set +u; conda deactivate || true; set -u; }            # conda scripts are not nounset-safe

# Run one step under one env, tee'd to its own log, and record the rc in $STATUS.
# A failing tool never aborts the driver; its rc is left in $LAST_RC.
LAST_RC=0
step() {
  local name="$1" env_prefix="$2"; shift 2
  LAST_RC=0
  want "$name" || return 0
  local log="$LOG_DIR/${name}_$TS.log"
  local saved_path="$PATH"
  echo ""
  echo "[$(date '+%F %T')] $name (env=$(basename "$env_prefix")) -> $log"
  act "$env_prefix"
  export PATH="$env_prefix/bin:$PATH"
  "$@" > "$log" 2>&1
  local rc=$?
  LAST_RC=$rc
  deact
  export PATH="$saved_path"          # else each step's bin/ stacks up on the next
  printf '%s\t%s\t%s\n' "$(date '+%F %T')" "$name" "$rc" >> "$STATUS"
  echo "  rc: $rc"
}

# A master_pipeline.py tool step. $2/$3 pick the reference strategy.
master_step() {
  python3 "$MASTER" \
    --root_dir "$BATCH" \
    --step "$1" \
    --tool_out_name "$TOOL_OUT_NAME" \
    --std_in_name "$2" \
    --tool_out_suffix "$3"
}
ref_tool()  { master_step "$1" "$STD_NORMAL" "$NORMAL_SUFFIX"; }
free_tool() { master_step "$1" "$STD_BASE"   ""; }

tumorfree_prep() {
  Rscript "$PREP_TUMORFREE" \
    --rep_root "$BATCH" --batch_root "$BATCH" --ref_dir "$REF_DIR" \
    --reps "$REP_IDS" --tumor_prop_cutoff "$TUMOR_PROP_CUTOFF" \
    --normal_std_in "$STD_NORMAL" --force "$FORCE_PREP"
  local rc=$?
  # The R side exits with a SIGSEGV at teardown (arrow/libstdc++) often enough
  # that rc alone is not a verdict; the summary csv only exists if it finished.
  if [ "$rc" -ne 0 ] && [ -s "$BATCH/tumorfree_ref_input_summary.csv" ]; then
    echo "[prep_tumorfree] rc=$rc but tumorfree_ref_input_summary.csv is present; continuing."
    return 0
  fi
  return "$rc"
}

infercnv_hmm17() {
  local fail=0
  for rep in "$BATCH"/Rep_*; do
    [ -d "$rep" ] || continue
    local out_dir="$rep/$TOOL_OUT_NAME/$INFERCNV_OUT"
    if [ -f "$out_dir/infercnv.OK" ]; then
      echo "[infercnv] skip $(basename "$rep") (infercnv.OK present)"
      continue
    fi
    mkdir -p "$out_dir"
    echo "[infercnv] $(basename "$rep") -> $out_dir"
    Rscript "$INFERCNV_HMM17" "$rep/$STD_NORMAL/InferCNV" "$out_dir" "$GENE_ORDER_FILE" --HMM TRUE
    [ -f "$out_dir/infercnv.OK" ] || { echo "[infercnv] FAILED: $(basename "$rep")"; fail=1; }
  done
  return "$fail"
}

# ---- run ---------------------------------------------------------------------
# The prep steps are prerequisites for all ten tools, so these two do abort.
step prep           "$ENV_R"          master_step prep "$STD_BASE" ""
[ "$LAST_RC" -eq 0 ] || { echo "FATAL: prep failed (rc=$LAST_RC) — see $LOG_DIR/prep_$TS.log" >&2; exit 1; }
step prep_tumorfree "$ENV_R"          tumorfree_prep
[ "$LAST_RC" -eq 0 ] || { echo "FATAL: prep_tumorfree failed (rc=$LAST_RC) — see $LOG_DIR/prep_tumorfree_$TS.log" >&2; exit 1; }

step numbat         "$ENV_NUMBAT"     ref_tool numbat
step copykat        "$ENV_R"          ref_tool copykat
step scevan         "$ENV_R"          ref_tool scevan
step conicsmat      "$ENV_R"          ref_tool conicsmat
step infercnvpy     "$ENV_INFERCNVPY" ref_tool infercnvpy
step xclone         "$ENV_XCLONE"     ref_tool xclone
step infercnv       "$ENV_R"          infercnv_hmm17

step starch         "$ENV_STARCH"     free_tool starch
step calicost       "$ENV_CALICOST"   free_tool calicost
step copyvae        "$ENV_COPYVAE"    free_tool copyvae

# The eval driver activates its own three envs, so it runs outside step().
if want eval; then
  echo ""
  echo "[$(date '+%F %T')] eval -> $ANALYSIS_DIR/run_normalref_extract_eval.sh"
  FLAT="$BATCH" RESULTS="$RESULTS" bash "$EVAL_SH"
  printf '%s\t%s\t%s\n' "$(date '+%F %T')" "eval" "$?" >> "$STATUS"
fi

echo ""
echo "=== DONE $(date '+%F %T') ==="
column -t -s $'\t' "$STATUS"
awk -F'\t' 'NR>1 && $3!=0 {n++} END {if (n) printf "\n%d step(s) returned non-zero — see %s\n", n, "'"$LOG_DIR"'"}' "$STATUS"
