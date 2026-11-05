#!/bin/bash
#SBATCH --job-name=009_perm
#SBATCH --output=src/logs/009_perm_%A_%a.log
#SBATCH --error=src/logs/009_perm_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=high_p
#SBATCH --open-mode=truncate
#
# One array task = one block of permutations. Submit with, e.g.
#   sbatch --array=0-19 src/009_permutation/run_permutation.sh --per-task 50
# for 1,000 permutations. METHODS.md 8.3.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="${SLURM_SUBMIT_DIR:-${GWAS_PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}}"
_LIB="${PROJ_ROOT}/src/lib/env.sh"
[[ -f "$_LIB" ]] || { echo "ERROR: ${_LIB} not found; sync the tree." >&2; exit 1; }
# shellcheck disable=SC1090
source "$_LIB"; load_env_file "${PROJ_ROOT}/.env" || true

PDIR="${PROJ_ROOT}/analysis/permutation"
ANNOT=""; PER=50; SEED=20260918
while [[ $# -gt 0 ]]; do
  case "$1" in
    --perm-dir) PDIR="$2"; shift 2 ;;
    --annot)    ANNOT="$2"; shift 2 ;;
    --per-task) PER="$2"; shift 2 ;;
    --seed)     SEED="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done
if [[ -z "$ANNOT" ]]; then
  mapfile -t C < <(ls -1 "${GWAS_PHENO_DIR:-/nonexistent}"/*.annot.tsv 2>/dev/null)
  [[ ${#C[@]} -eq 1 ]] && ANNOT="${C[0]}"
fi
[[ -n "$ANNOT" && -f "$ANNOT" ]] || { echo "ERROR: --annot required (panel annotation)" >&2; exit 1; }

K="${SLURM_ARRAY_TASK_ID:-0}"
START=$(( K * PER ))
echo "block $K: permutations ${START}..$(( START + PER - 1 ))   seed $SEED"
echo "perm-dir $PDIR"
echo "annot    $ANNOT"

if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  ENV_DIR="${GWAS_ENV_DIR:-${PROJ_ROOT}/../env}"
  [[ -d "${ENV_DIR}/herit" ]] && conda activate "${ENV_DIR}/herit"
fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

python3 "${PROJ_ROOT}/src/009_permutation/permute.py" \
  --perm-dir "$PDIR" --annot "$ANNOT" \
  --start "$START" --n "$PER" --seed "$SEED"
