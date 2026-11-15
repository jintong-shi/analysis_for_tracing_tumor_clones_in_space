#!/usr/bin/env bash
# n_clones sensitivity sweep for STARCH + CalicoST at n_clones = 2,3,5,6.
# n_clones = 4 is not re-run; its results are reused from $BASELINE_ARI.
# Each K gets a shadow root under $SWEEP_DIR/roots/nclones_K/ with synthetic/ and
# Standardized_Inputs/ symlinked from $FLAT, so the production batch is untouched.
# The STARCH/CalicoST invocations are mirrored from master_pipeline.py — keep in sync.
# Idempotent: finished (rep, K, tool) work is skipped on re-run.
# Usage:  nohup bash run_nclones_sensitivity.sh > sweep.out 2>&1 &
#         SMOKE=1 KVALS=2 bash run_nclones_sensitivity.sh      # one rep, one K
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"    # this analysis/ folder

# ---- inputs / substrate ------------------------------------------------------
SYN_RUNS="<DATA_ROOT>/synthetic/numbat_dirichlet_fix_test/runs"
FLAT="${FLAT:-$SYN_RUNS/complete_batch_normalref}"            # production flat batch (n=4 baseline + inputs)
SWEEP_DIR="${SWEEP_DIR:-$SYN_RUNS/nclones_sensitivity}"       # all sweep output lands here
BASELINE_ARI="${BASELINE_ARI:-$FLAT/results/subclone_ari_summary.csv}"  # existing n_clones=4 results (reused)
SWEEP_ROOTS="$SWEEP_DIR/roots"                                # per-K shadow roots
RESULTS="$SWEEP_DIR/results"; mkdir -p "$RESULTS"
STATUS="$SWEEP_DIR/status.tsv"

# ---- tool executables / resources (mirrored from master_pipeline.py) ---------
STARCH_DIR="${STARCH_DIR:-<TOOLS_ROOT>/STARCH}"
STARCH_EXEC="$STARCH_DIR/run_STARCH.py"
STARCH_REF="$STARCH_DIR/hgTables_hg19.txt"
CALICOST_PY="${CALICOST_PY:-<CONDA_ROOT>/envs/calicost_env/bin/python}"
CALICOST_SRC="${CALICOST_SRC:-<TOOLS_ROOT>/CalicoST/src}"
CALICOST_PARSE="$CALICOST_SRC/calicost/parse_input.py"
CALICOST_MAIN="$CALICOST_SRC/calicost/calicost_main.py"
EVAL_R="$SCRIPT_DIR/extract_subclone_and_ari.R"
SUMMARIZE_R="$SCRIPT_DIR/summarize_nclones_sensitivity.R"

# ---- conda envs --------------------------------------------------------------
# ENV_STARCH may live under a separate mamba prefix — override if it is elsewhere.
CONDA_SH="${CONDA_SH:-<CONDA_ROOT>/etc/profile.d/conda.sh}"
ENV_STARCH="${ENV_STARCH:-<CONDA_ROOT>/envs/starch}"
ENV_CALICOST="${ENV_CALICOST:-<CONDA_ROOT>/envs/calicost_env}"
ENV_RSHARED="${ENV_RSHARED:-<CONDA_ROOT>/envs/stcnv-r-shared}"
# shellcheck disable=SC1090
source "$CONDA_SH" || { echo "FATAL: cannot source $CONDA_SH" >&2; exit 1; }
act(){ set +u; conda activate "$1" || { set -u; echo "FATAL: conda activate $1 failed" >&2; exit 1; }; set -u; }  # conda scripts are not nounset-safe

# ---- knobs -------------------------------------------------------------------
KVALS="${KVALS:-2 3 5 6}"                                     # 4 reused from baseline, not swept
TOOLS="${TOOLS:-STARCH CalicoST}"
THREADS="${THREADS:-1}"
POOL="${POOL:-6}"                                             # max concurrent reps per (tool,K)

# ---- discover reps -----------------------------------------------------------
if [ -n "${REPS:-}" ]; then read -r -a REP_LIST <<< "$REPS"
else
  REP_LIST=()
  for d in "$FLAT"/Rep_*; do [ -d "$d" ] && REP_LIST+=("$(basename "$d")"); done
fi
[ "${SMOKE:-0}" = "1" ] && REP_LIST=("${REP_LIST[0]}")
[ "${#REP_LIST[@]}" -gt 0 ] || { echo "FATAL: no reps under $FLAT" >&2; exit 1; }

echo "=== n_clones sensitivity | host $(hostname) | $(date '+%F %T') ==="
echo "  flat batch : $FLAT"
echo "  sweep dir  : $SWEEP_DIR"
echo "  KVALS      : $KVALS   (n_clones=4 reused from $BASELINE_ARI)"
echo "  tools      : $TOOLS"
echo "  reps (${#REP_LIST[@]}) : ${REP_LIST[*]}"
[ -f "$BASELINE_ARI" ] || echo "  WARNING: baseline n=4 ARI missing ($BASELINE_ARI) — aggregation will omit n=4."
[ -f "$STATUS" ] || printf 'date\tn_clones\ttool\trep\texit\twall\tpeak_RAM_GB\tlog\n' > "$STATUS"

record(){ # $1=n $2=tool $3=rep $4=exit $5=log
  local log="$5" wall ram="NA" kb
  wall=$(grep -aoE "Elapsed \(wall clock\) time.*" "$log" 2>/dev/null | sed 's/.*): //' | tail -1)
  kb=$(grep -aoE "Maximum resident set size \(kbytes\): [0-9]+" "$log" 2>/dev/null | grep -oE "[0-9]+$" | tail -1)
  [ -n "${kb:-}" ] && ram=$(awk -v k="$kb" 'BEGIN{printf "%.2f", k/1024/1024}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date '+%F %T')" "$1" "$2" "$3" "$4" "${wall:-NA}" "$ram" "$log" >> "$STATUS"
}

# shadow root for a (K, rep): symlink synthetic + Standardized_Inputs, real Tool_Outputs
make_shadow(){ # $1=K $2=rep -> echoes rep root path
  local K="$1" rep="$2" sr="$SWEEP_ROOTS/nclones_$1/$2"
  mkdir -p "$sr/Tool_Outputs"
  for sub in synthetic Standardized_Inputs; do
    [ -e "$sr/$sub" ] || ln -s "$FLAT/$rep/$sub" "$sr/$sub"
  done
  echo "$sr"
}

# ---- per-tool runners --------------------------------------------------------
run_starch(){ # $1=K $2=rep
  local K="$1" rep="$2" sr; sr="$(make_shadow "$K" "$rep")"
  local csv="$FLAT/$rep/Standardized_Inputs/STARCH/gene_expression_matrix.csv"
  local out="$sr/Tool_Outputs/STARCH" log="$sr/tool_starch.log"
  [ -f "$out/labels_STITCH_output.csv" ] && { echo "[STARCH K$K] $rep done"; return 0; }
  [ -f "$csv" ] || { echo "[STARCH K$K] $rep MISSING csv"; record "$K" STARCH "$rep" NA "$log"; return 1; }
  mkdir -p "$out"
  /usr/bin/time -v python3 "$STARCH_EXEC" -c "$K" -t "$THREADS" -i "$csv" \
      --outdir "$out" -m "$STARCH_REF" > "$log" 2>&1
  record "$K" STARCH "$rep" "$?" "$log"
}

run_calicost(){ # $1=K $2=rep
  local K="$1" rep="$2" sr; sr="$(make_shadow "$K" "$rep")"
  local snp="$FLAT/$rep/Standardized_Inputs/CalicoST/snpinfo"
  local out="$sr/Tool_Outputs/CalicoST" log="$sr/tool_calicost.log"
  local base_cfg="$FLAT/$rep/Tool_Outputs/CalicoST/configfile0"
  ls "$out"/clone*_rectangle0_w1.0/clone_labels.tsv >/dev/null 2>&1 && { echo "[CalicoST K$K] $rep done"; return 0; }
  [ -f "$base_cfg" ] || { echo "[CalicoST K$K] $rep MISSING base config"; record "$K" CalicoST "$rep" NA "$log"; return 1; }
  [ -d "$snp" ]      || { echo "[CalicoST K$K] $rep MISSING snpinfo"; record "$K" CalicoST "$rep" NA "$log"; return 1; }
  mkdir -p "$out"
  # clone the rep's proven n=4 config; change ONLY n_clones, output_dir, snp paths
  local cfg="$out/configfile0"
  sed -e "s#^n_clones :.*#n_clones : $K#" \
      -e "s#^spaceranger_dir :.*#spaceranger_dir : $snp#" \
      -e "s#^snp_dir :.*#snp_dir : $snp#" \
      -e "s#^output_dir :.*#output_dir : $out#" \
      "$base_cfg" > "$cfg"
  ( export PYTHONPATH="$CALICOST_SRC${PYTHONPATH:+:$PYTHONPATH}"
    /usr/bin/time -v bash -c "\"$CALICOST_PY\" \"$CALICOST_PARSE\" -c \"$cfg\" && \
                              \"$CALICOST_PY\" \"$CALICOST_MAIN\" -c \"$cfg\"" > "$log" 2>&1 )
  record "$K" CalicoST "$rep" "$?" "$log"
}

pool(){ # $1=fn $2=K ; runs fn K rep for all reps, <=POOL at a time
  local fn="$1" K="$2"
  for rep in "${REP_LIST[@]}"; do
    while [ "$(jobs -rp | wc -l)" -ge "$POOL" ]; do sleep 10; done
    "$fn" "$K" "$rep" &
  done
  wait
}

want(){ case " $TOOLS " in *" $1 "*) return 0;; *) return 1;; esac; }

# ---- execute -----------------------------------------------------------------
for K in $KVALS; do
  echo "### n_clones = $K  ($(date '+%F %T'))"
  if want STARCH;   then echo "  -- STARCH --";   act "$ENV_STARCH";   pool run_starch   "$K"; fi
  if want CalicoST; then echo "  -- CalicoST --"; act "$ENV_CALICOST"; pool run_calicost "$K"; fi

  echo "  -- eval (04) n_clones=$K --"
  act "$ENV_RSHARED"
  Rscript "$EVAL_R" --root_dir "$SWEEP_ROOTS/nclones_$K" \
      --tool_out_name Tool_Outputs --std_in_name Standardized_Inputs \
      > "$SWEEP_ROOTS/nclones_$K/eval_04.log" 2>&1 \
    && echo "    -> $SWEEP_ROOTS/nclones_$K/results/subclone_ari_summary.csv" \
    || echo "    [WARN] 04 eval failed for K=$K (see eval_04.log)"
done

# ---- aggregate across n_clones (+ reuse n=4 baseline) ------------------------
echo "### aggregate  ($(date '+%F %T'))"
act "$ENV_RSHARED"
Rscript "$SUMMARIZE_R" \
    --sweep_dir "$SWEEP_DIR" --baseline_ari "$BASELINE_ARI" > "$SWEEP_DIR/summarize.log" 2>&1 \
  && { echo "=== DONE ; summaries in $RESULTS ==="; ls -1 "$RESULTS"; } \
  || echo "[WARN] aggregation failed (see summarize.log)"