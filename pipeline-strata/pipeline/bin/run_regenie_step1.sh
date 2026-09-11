#!/bin/bash
set -euo pipefail

step1_pfile=""
pheno_file=""
cov_file=""
master_file=""
phenotype_id=""
study=""
group=""
bsize=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --step1_pfile) step1_pfile="$2"; shift 2 ;;
    --pheno_file) pheno_file="$2"; shift 2 ;;
    --cov_file) cov_file="$2"; shift 2 ;;
    --master_file) master_file="$2"; shift 2 ;;
    --phenotype_id) phenotype_id="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --bsize) bsize="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$step1_pfile" || -z "$pheno_file" || -z "$cov_file" || -z "$master_file" || -z "$phenotype_id" || -z "$study" || -z "$group" || -z "$bsize" || -z "$outdir" ]]; then
  echo "Usage: $0 --step1_pfile <prefix> --pheno_file <file> --cov_file <file> --master_file <file> --phenotype_id <id> --study <id> --group <group> --bsize <val> --outdir <dir>" >&2
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${step1_pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${step1_pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"

log_file="${outdir}/step1.log"
threads="${NXF_TASK_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS="${threads}"
export OPENBLAS_NUM_THREADS="${threads}"
export MKL_NUM_THREADS="${threads}"
export NUMEXPR_MAX_THREADS="${threads}"

{
  echo "COVARIATE FILE PREVIEW (first 2 lines):"
  head -n 2 "${cov_file}" 2>/dev/null || echo "WARNING: failed to read covariate file '${cov_file}'"
  echo ""
} > "$log_file"

cmd=(
  regenie
  --step 1
  --pgen "${step1_pfile}"
  --phenoFile "${pheno_file}"
  --covarFile "${cov_file}"
  --bsize "${bsize}"
  --threads "${threads}"
  --gz
  --force-step1
  --run-l1 "${master_file}"
  --l1-phenoList "${phenotype_id}"
  --keep-l0
  --out "${outdir}/step1"
)

printf 'Running: %q ' "${cmd[@]}" >> "$log_file"
printf '\n' >> "$log_file"

if ! "${cmd[@]}" >> "$log_file" 2>&1; then
    echo "ERROR: REGENIE Step 1 failed. Log follows:"
    cat "$log_file"
    exit 1
fi

if [[ -f "${outdir}/step1_pred.list" ]]; then
  mv "${outdir}/step1_pred.list" "${outdir}/pred.list"
fi

for f in "${outdir}"/step1_*.loco.gz; do
  [[ -e "$f" ]] || continue
  base="${f##*/step1_}"
  mv "$f" "${outdir}/loco_${base}"
done

# Update pred.list to use the renamed LOCO files and keep paths relative so the
# file stays valid when Nextflow stages it into downstream task directories.
if [[ -f "${outdir}/pred.list" ]]; then
  python3 - "${outdir}/pred.list" <<'PY'
import os
import re
import sys

pred_path = sys.argv[1]
tmp_path = pred_path + ".tmp"

pattern = re.compile(r'\S+\.loco\.gz')

def rewrite_token(match: re.Match[str]) -> str:
    token = match.group(0)
    base = os.path.basename(token)
    if base.startswith("step1_") and base.endswith(".loco.gz"):
        base = "loco_" + base[len("step1_"):]
    if base.startswith("loco_") and base.endswith(".loco.gz"):
        return f"./{base}"
    return token

with open(pred_path, "r", encoding="utf-8") as src, open(tmp_path, "w", encoding="utf-8") as dst:
    for raw_line in src:
        dst.write(pattern.sub(rewrite_token, raw_line))

os.replace(tmp_path, pred_path)
PY
fi

echo "REGENIE Step 1 for $study ($group) phenotype $phenotype_id complete."
