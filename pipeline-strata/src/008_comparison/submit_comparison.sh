#!/bin/bash
# src/008_comparison/submit_comparison.sh
#
# Run this ON THE LOGIN NODE (bash, not sbatch). It submits two jobs:
#
#   1. a MAP array  -- one task per shard, each comparing a stride of the
#      proteins and writing <out>/shards/part_NNNNN.pkl
#   2. a REDUCE job -- --dependency=afterok on the array, which combines the
#      shards, runs the global tests, and writes the final tables
#
# The pool inside run_comparison.sh is capped at one node's cores. This is not:
# proteins are independent, so with one task per protein the wall time is the
# SLOWEST PROTEIN plus queue, not the sum over 261 of them.
#
# The reduce step cannot be sharded -- the Bonferroni denominator and the
# var(z)/F tests need every protein in one frame -- but it is seconds of work.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="${GWAS_PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"

_GWAS_ENV_LIB="${PROJ_ROOT}/src/lib/env.sh"
[[ -f "$_GWAS_ENV_LIB" ]] || _GWAS_ENV_LIB="${SCRIPT_DIR}/../lib/env.sh"
[[ -f "$_GWAS_ENV_LIB" ]] || {
  echo "ERROR: src/lib/env.sh not found. The tree is incompletely synced;" >&2
  echo "       rsync pipeline-strata/ again (see TRANSFER.md)." >&2
  exit 1
}
# shellcheck disable=SC1090
source "$_GWAS_ENV_LIB"
load_env_file "${PROJ_ROOT}/.env" || true
mkdir -p "${PROJ_ROOT}/src/logs"

STAGED="${GWAS_STAGING_DIR:-}"
OUT=""
ANNOT=""
PER_TASK=1                       # proteins per array task
THROTTLE="${GWAS_ARRAY_THROTTLE:-}"   # empty = no %limit, use the whole cluster
PARTITION="${GWAS_PARTITION:-high_p}"
MAP_TIME="2:00:00"
MAP_MEM="8G"
DRY=0
PASS=()

usage() {
  sed -n '2,20p' "${BASH_SOURCE[0]:-$0}"
  cat <<'EOF'

Options:
  --staged DIR        staged tree (default: $GWAS_STAGING_DIR)
  --out DIR           output directory (default: <project>/analysis/comparison)
  --annot FILE        <panel>.annot.tsv (required for the random-proteome null)
  --per-task N        proteins per array task (default 1 = maximum parallelism)
  --throttle N        cap concurrent tasks as --array=...%N (default: uncapped)
  --time HH:MM:SS     walltime per map task (default 2:00:00)
  --mem SIZE          memory per map task (default 8G)
  --dry-run           print the sbatch commands without submitting
  anything else is forwarded to compare_strata.py
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --staged)   STAGED="$2"; shift 2 ;;
    --out)      OUT="$2"; shift 2 ;;
    --annot)    ANNOT="$2"; shift 2 ;;
    --per-task) PER_TASK="$2"; shift 2 ;;
    --throttle) THROTTLE="$2"; shift 2 ;;
    --time)     MAP_TIME="$2"; shift 2 ;;
    --mem)      MAP_MEM="$2"; shift 2 ;;
    --dry-run)  DRY=1; shift ;;
    -h|--help)  usage; exit 0 ;;
    *)          PASS+=("$1"); shift ;;
  esac
done

[[ -n "$STAGED" ]] || { echo "ERROR: --staged required (or GWAS_STAGING_DIR)" >&2; exit 1; }
[[ -d "$STAGED" ]] || { echo "ERROR: not a directory: $STAGED" >&2; exit 1; }
[[ -n "$OUT" ]] || OUT="${PROJ_ROOT}/analysis/comparison"

# Count proteins exactly as compare_strata.py does, so the shard arithmetic
# here and the striding there agree.
N_PROT=$(find "$STAGED" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
(( N_PROT > 0 )) || { echo "ERROR: no protein directories under $STAGED" >&2; exit 1; }
(( PER_TASK >= 1 )) || { echo "ERROR: --per-task must be >= 1" >&2; exit 1; }

N_SHARDS=$(( (N_PROT + PER_TASK - 1) / PER_TASK ))
ARRAY_SPEC="0-$(( N_SHARDS - 1 ))"
[[ -n "$THROTTLE" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${THROTTLE}"

echo "=========================================="
echo " 008_comparison  (map/reduce)"
echo "=========================================="
echo "Staged     : $STAGED"
echo "Out        : $OUT"
echo "Proteins   : $N_PROT"
echo "Map tasks  : $N_SHARDS  (--array=${ARRAY_SPEC}, ${PER_TASK} protein(s) each)"
echo "Per task   : 1 cpu, ${MAP_MEM}, ${MAP_TIME}, partition ${PARTITION}"
echo "Annotation : ${ANNOT:-<none -- no random-proteome null>}"
echo "=========================================="

if [[ -z "$ANNOT" ]]; then
  echo "WARNING: without --annot the random proteome is not labelled, so" >&2
  echo "         questions 1 and 3 cannot be answered." >&2
fi

mkdir -p "$OUT/shards"

# A stale shard from an earlier run would be silently folded into the gather.
if compgen -G "$OUT/shards/part_*.pkl" >/dev/null; then
  echo "NOTE: clearing $(ls -1 "$OUT"/shards/part_*.pkl | wc -l | tr -d ' ') shard(s) from a previous run"
  (( DRY )) || rm -f "$OUT"/shards/part_*.pkl
fi

# run_comparison.sh resolves the project from SLURM_SUBMIT_DIR in preference to
# GWAS_PROJECT_ROOT, and SLURM_SUBMIT_DIR is just the cwd at submission. Submit
# from the project root so the two always agree; otherwise the array tasks look
# for .env and compare_strata.py under whatever directory this was launched in.
cd "$PROJ_ROOT"

MAP_SCRIPT="${SCRIPT_DIR}/run_comparison.sh"
COMMON=( --staged "$STAGED" --out "$OUT" )
[[ -n "$ANNOT" ]] && COMMON+=( --annot "$ANNOT" )
COMMON+=( "${PASS[@]+"${PASS[@]}"}" )

submit() {
  # stdout is captured by $(...) for the job id, so the dry-run trace has to go
  # to stderr or it would be swallowed and never shown.
  if (( DRY )); then
    { printf 'DRY-RUN: '; printf '%q ' "$@"; printf '\n\n'; } >&2
    return 0
  fi
  "$@"
}

MAP_ID=$(submit sbatch --parsable \
  --job-name=008_map \
  --array="$ARRAY_SPEC" \
  --partition="$PARTITION" \
  --time="$MAP_TIME" \
  --mem="$MAP_MEM" \
  --cpus-per-task=1 \
  --output="${PROJ_ROOT}/src/logs/008_map_%A_%a.log" \
  --error="${PROJ_ROOT}/src/logs/008_map_%A_%a.err" \
  "$MAP_SCRIPT" "${COMMON[@]}" --shard-from-array --n-shards "$N_SHARDS")

if (( DRY )); then MAP_ID="<mapid>"; fi
echo "Map array  : $MAP_ID"

GATHER_ID=$(submit sbatch --parsable \
  --job-name=008_gather \
  --dependency="afterok:${MAP_ID}" \
  --partition="$PARTITION" \
  --time=1:00:00 \
  --mem=32G \
  --cpus-per-task=1 \
  --output="${PROJ_ROOT}/src/logs/008_gather.log" \
  --error="${PROJ_ROOT}/src/logs/008_gather.err" \
  "$MAP_SCRIPT" "${COMMON[@]}" --gather --n-shards "$N_SHARDS")

if (( DRY )); then GATHER_ID="<gatherid>"; fi
echo "Gather     : $GATHER_ID  (afterok:${MAP_ID})"
echo
echo "Watch:   squeue -j ${MAP_ID},${GATHER_ID}"
echo "Results: $OUT/pqtl_comparison.tsv  (after gather)"
echo "Report:  ${PROJ_ROOT}/src/logs/008_gather.log"
