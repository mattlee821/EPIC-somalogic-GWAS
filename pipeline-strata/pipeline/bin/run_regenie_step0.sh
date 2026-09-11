#!/bin/bash
set -euo pipefail

step1_pfile=""
pheno_file=""
cov_file=""
study=""
group=""
bsize=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --step1_pfile) step1_pfile="$2"; shift 2 ;;
    --pheno_file) pheno_file="$2"; shift 2 ;;
    --cov_file) cov_file="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --bsize) bsize="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$step1_pfile" || -z "$pheno_file" || -z "$cov_file" || -z "$study" || -z "$group" || -z "$bsize" || -z "$outdir" ]]; then
  echo "Usage: $0 --step1_pfile <prefix> --pheno_file <file> --cov_file <file> --study <id> --group <group> --bsize <val> --outdir <dir>" >&2
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${step1_pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${step1_pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"
outdir="$(cd "$outdir" && pwd -P)"

threads="${NXF_TASK_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS="${threads}"
export OPENBLAS_NUM_THREADS="${threads}"
export MKL_NUM_THREADS="${threads}"
export NUMEXPR_MAX_THREADS="${threads}"
l0_jobs="${threads}"
if (( l0_jobs < 2 )); then
  l0_jobs=2
fi
master_prefix="${outdir}/regenie_step0"
master_file="${master_prefix}.master"
split_log="${outdir}/step0_split.log"

split_cmd=(
  regenie
  --step 1
  --pgen "${step1_pfile}"
  --phenoFile "${pheno_file}"
  --covarFile "${cov_file}"
  --bsize "${bsize}"
  --threads "${threads}"
  --gz
  --force-step1
  --out "${outdir}/step0_split"
  --split-l0 "${master_prefix},${l0_jobs}"
)

printf 'Running: %q ' "${split_cmd[@]}" > "$split_log"
printf '\n' >> "$split_log"
if ! "${split_cmd[@]}" >> "$split_log" 2>&1; then
  echo "ERROR: REGENIE Step 0 split-l0 failed. Log follows:"
  cat "$split_log"
  exit 1
fi

for job_id in $(seq 1 "${l0_jobs}"); do
  l0_log="${outdir}/step0_l0_job${job_id}.log"
  run_l0_cmd=(
    regenie
    --step 1
    --pgen "${step1_pfile}"
    --phenoFile "${pheno_file}"
    --covarFile "${cov_file}"
    --bsize "${bsize}"
    --threads "${threads}"
    --gz
    --force-step1
    --out "${outdir}/step0_l0_job${job_id}"
    --run-l0 "${master_file},${job_id}"
  )

  printf 'Running: %q ' "${run_l0_cmd[@]}" > "$l0_log"
  printf '\n' >> "$l0_log"
  if ! "${run_l0_cmd[@]}" >> "$l0_log" 2>&1; then
    echo "ERROR: REGENIE Step 0 run-l0 job ${job_id} failed. Log follows:"
    cat "$l0_log"
    exit 1
  fi
done

echo "REGENIE Step 0 for $study ($group) complete using ${l0_jobs} level-0 shards."
