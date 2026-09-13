#!/bin/bash
#SBATCH --job-name=008_comparison
#SBATCH --output=src/logs/008_comparison.log
#SBATCH --error=src/logs/008_comparison.err
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=high_p
#SBATCH --open-mode=truncate

# Compare pQTL effects between strata on the staged results.
#
#   sbatch src/008_comparison/run_comparison.sh --staged <dir> --out <dir>
#   bash   src/008_comparison/run_comparison.sh --staged <dir> --out <dir>
#
# Submit from the project root: the SBATCH --output paths above are relative to
# the submit directory.
#
# Reference data (LD panel, gene map, HLA window, cis window) is read straight
# out of pipeline/params.yaml so this cannot drift from what the pipeline used.
# Runs in the `herit` conda env, which already carries plink, scipy, pandas and
# numpy. Plots additionally need matplotlib; without it everything else still
# runs and the tables are written.

set -euo pipefail

# Under SLURM the job runs from the slurmd spool directory and this script is a
# COPY placed there, so BASH_SOURCE does not point into the project. Everything
# is therefore resolved from SLURM_SUBMIT_DIR when it is set.
SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="${SLURM_SUBMIT_DIR:-${GWAS_PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}}"
[[ -f "${PROJ_ROOT}/.env" ]] && { set -a; source "${PROJ_ROOT}/.env"; set +a; }
mkdir -p "${PROJ_ROOT}/src/logs"

STAGED="${GWAS_STAGING_DIR:-}"
OUT=""
STRATA="ever,never"
REFERENCE="meta"
PDISC="5e-8"
ANNOT=""
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --staged) STAGED="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --strata) STRATA="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --p-discovery) PDISC="$2"; shift 2 ;;
    --annot) ANNOT="$2"; shift 2 ;;
    -h|--help) sed -n '10,20p' "${BASH_SOURCE[0]:-$0}"; exit 0 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done

[[ -n "$STAGED" ]] || { echo "ERROR: --staged required (or GWAS_STAGING_DIR)" >&2; exit 1; }
[[ -d "$STAGED" ]] || { echo "ERROR: not a directory: $STAGED" >&2; exit 1; }
[[ -n "$OUT" ]] || OUT="${PROJ_ROOT}/analysis/comparison"

PARAMS="${GWAS_PARAMS_FILE:-${PROJ_ROOT}/pipeline/params.yaml}"
yml() { sed -n "s/^$1:[[:space:]]*\"\{0,1\}\([^\"#]*\)\"\{0,1\}.*/\1/p" "$PARAMS" | head -1 | xargs; }

LD_REF="$(yml ld_ref_pfile)"
GENE_MAP="$(yml protein_gene_map)"
CIS_KB="$(yml cis_window_kb)"
HLA="$(yml hla_region)"
CLUMP_R2="$(yml herit_clump_r2)"

# Find the panel annotation so random-control proteins are recognised as the
# empirical null. GWAS_INCLUDE_PROTEINS is only a default and is usually NOT the
# panel a given run used -- the 13 Sep run passed --proteins smoking-panel.txt
# while .env still pointed at test.txt, so the controls went unlabelled and the
# calibration check was silently lost. Fall back to globbing the phenotype
# directory, which is where find_proteins.R writes them.
if [[ -z "$ANNOT" && -n "${GWAS_INCLUDE_PROTEINS:-}" ]]; then
  CAND="${GWAS_INCLUDE_PROTEINS%.txt}.annot.tsv"
  [[ -f "$CAND" ]] && ANNOT="$CAND"
fi
if [[ -z "$ANNOT" && -n "${GWAS_PHENO_DIR:-}" ]]; then
  mapfile -t CANDS < <(ls -1 "${GWAS_PHENO_DIR}"/*.annot.tsv 2>/dev/null)
  if [[ ${#CANDS[@]} -eq 1 ]]; then
    ANNOT="${CANDS[0]}"
  elif [[ ${#CANDS[@]} -gt 1 ]]; then
    echo "NOTE: several .annot.tsv files under ${GWAS_PHENO_DIR}; pass --annot to choose:" >&2
    printf '  %s\n' "${CANDS[@]}" >&2
  fi
fi
if [[ -z "$ANNOT" ]]; then
  echo "WARNING: no panel annotation found, so no control set. The calibration" >&2
  echo "         check (var(z)~1 in controls) will be unavailable. Pass --annot." >&2
fi

echo "Staged:     $STAGED"
echo "Out:        $OUT"
echo "Strata:     $STRATA   (reference: $REFERENCE, P < $PDISC)"
echo "LD ref:     ${LD_REF:-<none, distance-based pruning>}"
echo "Gene map:   ${GENE_MAP:-<none, no cis/trans>}"
echo "Annotation: ${ANNOT:-<none, no control set>}"
echo

if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  ENV_DIR="${GWAS_ENV_DIR:-${PROJ_ROOT}/../env}"
  [[ -d "${ENV_DIR}/herit" ]] && conda activate "${ENV_DIR}/herit"
fi

# compare_strata.py lives beside this script in the project, NOT beside the
# spool copy SLURM runs.
CMP_PY="${PROJ_ROOT}/src/008_comparison/compare_strata.py"
[[ -f "$CMP_PY" ]] || CMP_PY="${SCRIPT_DIR}/compare_strata.py"
[[ -f "$CMP_PY" ]] || { echo "ERROR: cannot find compare_strata.py (looked in ${PROJ_ROOT}/src/008_comparison and ${SCRIPT_DIR})" >&2; exit 1; }

python3 "$CMP_PY" \
  --staged "$STAGED" \
  --out "$OUT" \
  --strata "$STRATA" \
  --reference "$REFERENCE" \
  --p-discovery "$PDISC" \
  ${LD_REF:+--ld-ref "$LD_REF"} \
  ${GENE_MAP:+--gene-map "$GENE_MAP"} \
  ${CIS_KB:+--cis-window-kb "$CIS_KB"} \
  ${HLA:+--hla "$HLA"} \
  ${CLUMP_R2:+--clump-r2 "$CLUMP_R2"} \
  ${ANNOT:+--annot "$ANNOT"} \
  ${EXTRA[@]+"${EXTRA[@]}"}
