#!/bin/bash
set -euo pipefail

pfile=""
keep=""
extract=""
study=""
window=""
step=""
r2=""
maf=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) pfile="$2"; shift 2 ;;
    --keep) keep="$2"; shift 2 ;;
    --extract) extract="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --window) window="$2"; shift 2 ;;
    --step) step="$2"; shift 2 ;;
    --r2) r2="$2"; shift 2 ;;
    --maf) maf="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$pfile" || -z "$keep" || -z "$extract" || -z "$study" || -z "$window" || -z "$step" || -z "$r2" || -z "$maf" || -z "$outdir" ]]; then
  echo "Usage: $0 --pfile <prefix> --keep <file> --extract <file> --study <val> --window <kb> --step <val> --r2 <val> --maf <val> --outdir <dir>" >&2
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"
log_file="${outdir}/ld_pruning.log"
exec > >(tee -a "$log_file") 2>&1

# Constrain plink2 to the SLURM allocation: by default plink2 reserves ~half the
# node's RAM (seen reserving ~190 GB), which both wastes time and risks a cgroup
# OOM kill; cap it to the requested memory and use the allocated CPUs.
threads="${SLURM_CPUS_PER_TASK:-2}"
mem_args=()
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "${SLURM_MEM_PER_NODE}" -gt 2048 ] 2>/dev/null; then
  mem_args=(--memory $(( SLURM_MEM_PER_NODE - 1024 )))
fi

echo "LD Pruning for $study using PLINK2 prefix $pfile (threads=${threads})..."

# --maf is applied AFTER --keep, i.e. computed on this study's analysis subset
# (not the full psam). This drops variants that are monomorphic / ultra-rare in
# the deduped sample, which would otherwise survive the upstream (full-sample)
# variant QC and crash REGENIE step 1 with a "low variance (=0.000000)" error.
plink2 --pfile "$pfile" \
  --keep "$keep" \
  --extract "$extract" \
  --maf "$maf" \
  --threads "$threads" ${mem_args[@]+"${mem_args[@]}"} \
  --indep-pairwise "${window}kb" "$step" "$r2" \
  --out "${outdir}/ld_pruned"

plink2 --pfile "$pfile" \
  --keep "$keep" \
  --extract "${outdir}/ld_pruned.prune.in" \
  --maf "$maf" \
  --threads "$threads" ${mem_args[@]+"${mem_args[@]}"} \
  --make-pgen \
  --out "${outdir}/step1_input"

echo "Pruned step1 input files written to ${outdir}/step1_input.*"
