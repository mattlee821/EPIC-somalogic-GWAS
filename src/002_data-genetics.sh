#!/bin/bash
#SBATCH --job-name=002_data-genetics
#SBATCH --output=src/logs/002_data-genetics_%a.log
#SBATCH --error=src/logs/002_data-genetics_%a.err
#SBATCH --array=0-22
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
# Script: src/002_data-genetics.sh
# Purpose: Stage imputation-pipeline output into a single, analysis-ready PLINK2
#          dataset per study, with all upstream exclusions applied in one pass.
#
# Input layout (one sub-directory per study under GWAS_SOURCE_GENETICS_DIR):
#   <study>/
#     <study>.tar.gz            -> <study>/<study>_chr{1..22,X}.{pgen,pvar}
#                                  + a single shared <study>.psam
#     report-*.html             (provenance, not consumed here)
#     review/
#       ancestry.exclude        sample IDs (FID IID) - ancestry/PCA outliers   -> REMOVE
#       related.exclude         sample IDs (FID IID) - KING relateds           -> REMOVE
#       hwe.exclude             variant IDs          - HWE failures (controls) -> EXCLUDE
#       het.het                 per-sample F coefficient (raw)                 -> REPORT ONLY
#       king.kin0               KING pairs                                     (reference)
#       pca.eigenvec/.eigenval  10 PCs computed upstream                       -> STAGED for reuse
#
# Output: one directory per study under GWAS_GENETICS_DIR/<study>/ containing
#   <study>.pgen/.pvar/.psam    one merged all-chromosome PLINK2 dataset
#   <study>.eigenvec/.eigenval  upstream PCs, staged for use as covariates
#   <study>.het_flagged.txt     samples beyond mean +/- 3 SD of F (NOT removed)
# The pfile prefix for the pipeline samplesheet is therefore
#   GWAS_GENETICS_DIR/<study>/<study>
#
# Heterozygosity policy: the imputation pipeline deliberately does not remove
# het outliers (it only acts on ancestry + relatedness). We mirror that here:
# het outliers are reported for inspection only, never excluded.
#
# Array job: one task per study. Submit:
#   sbatch src/002_data-genetics.sh
# Single study locally (no SLURM):
#   bash src/002_data-genetics.sh Glbd_01
# All studies serially (no SLURM):
#   bash src/002_data-genetics.sh --all

set -euo pipefail

# --- environment -------------------------------------------------------------
# Under SLURM the job runs from the slurmd spool dir, so look in SLURM_SUBMIT_DIR
# (where sbatch was launched) before the PWD/script-relative locations.
for env_file in "${SLURM_SUBMIT_DIR:-/nonexistent}/.env" ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

SRC_BASE="${GWAS_SOURCE_GENETICS_DIR:?Error: GWAS_SOURCE_GENETICS_DIR not set in .env}"
DIR_OUT="${GWAS_GENETICS_DIR:?Error: GWAS_GENETICS_DIR not set in .env}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
HET_SD="${GWAS_HET_SD:-3}"   # SD multiplier for het flag-only report

# Tell plink2 how much memory it may use (leave ~2 GB headroom). SLURM exports
# SLURM_MEM_PER_NODE in MB when --mem is set; if absent, plink2 auto-detects.
PLINK_MEM_ARGS=()
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "${SLURM_MEM_PER_NODE}" -gt 2048 ] 2>/dev/null; then
  PLINK_MEM_ARGS=(--memory $(( SLURM_MEM_PER_NODE - 2048 )))
fi

SRC_BASE="${SRC_BASE%/}"
DIR_OUT="${DIR_OUT%/}"

# --- locate plink2 -----------------------------------------------------------
# plink2 is provided by the prebuilt 'plink2' conda env created by
# src/001_source.sh under <root>/env/plink2 (it is NOT in nf_master). Activate the
# prefix env first, falling back to a legacy named 'plink2' env if present.
if ! command -v plink2 >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    plink2_prefix="${GWAS_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/env/plink2"
    if [ -d "$plink2_prefix" ]; then
      conda activate "$plink2_prefix" 2>/dev/null || true
    fi
    command -v plink2 >/dev/null 2>&1 || conda activate plink2 2>/dev/null || true
  fi
fi
if ! command -v plink2 >/dev/null 2>&1; then
  echo "ERROR: plink2 not found on PATH. Run 'sbatch src/001_source.sh' to build" >&2
  echo "       the env/plink2 conda env, or 'conda activate env/plink2' before running." >&2
  exit 1
fi

# --- discover studies from the source directory ------------------------------
# A directory qualifies as a study if it contains <study>/<study>.tar.gz.
STUDIES=()
for d in "$SRC_BASE"/*/; do
  [ -d "$d" ] || continue
  s="$(basename "$d")"
  [ -f "${d}${s}.tar.gz" ] || continue
  STUDIES+=("$s")
done
IFS=$'\n' STUDIES=($(sort <<<"${STUDIES[*]:-}")); unset IFS
if [ "${#STUDIES[@]}" -eq 0 ]; then
  echo "ERROR: no study sub-directories (with <study>.tar.gz) found under ${SRC_BASE}" >&2
  exit 1
fi

# =============================================================================
# process_study <study_name>
# =============================================================================
process_study() {
  local study="$1"
  local src_study="${SRC_BASE}/${study}"
  local tarball="${src_study}/${study}.tar.gz"
  local review="${src_study}/review"
  local out_study="${DIR_OUT}/${study}"   # final per-study directory
  local geno="${out_study}/_extract"      # temporary extraction area (removed at the end)

  echo "START: study=${study}"

  [ -f "$tarball" ] || { echo "ERROR: tarball not found: ${tarball}" >&2; return 1; }
  [ -d "$review" ]  || { echo "ERROR: review dir not found: ${review}" >&2; return 1; }

  rm -rf "$out_study"
  mkdir -p "$geno"

  # 1. Extract the archive directly under the study directory. The tarball has a
  #    top-level <study>/ folder, so --strip-components=1 places the per-chromosome
  #    files at ${geno}/<study>_chr*.{pgen,pvar} and the shared ${geno}/<study>.psam.
  #    Use pigz to decompress when available (piped, so its flags are parsed by
  #    the shell rather than by tar's --use-compress-program); fall back to gzip.
  if command -v pigz >/dev/null 2>&1; then
    pigz -dc -p "${THREADS}" "$tarball" | tar -xf - -C "$geno" --strip-components=1
  else
    tar -xzf "$tarball" -C "$geno" --strip-components=1
  fi

  local shared_psam="${geno}/${study}.psam"
  [ -f "$shared_psam" ] || { echo "ERROR: shared psam missing after extraction: ${shared_psam}" >&2; return 1; }

  # 2. Build the merge list. Each per-chromosome fileset is given as three
  #    explicit tokens (pgen, pvar, shared psam) so we avoid the mode-keyword
  #    ambiguity of bare-prefix lists and need no per-chromosome psam symlinks.
  local merge_list="${geno}/merge_list.txt"
  : > "$merge_list"
  local pgen pfx n_chr=0
  for pgen in "${geno}/${study}"_chr*.pgen; do
    [ -e "$pgen" ] || continue
    pfx="${pgen%.pgen}"
    [ -f "${pfx}.pvar" ] || { echo "ERROR: missing ${pfx}.pvar" >&2; return 1; }
    printf '%s\t%s\t%s\n' "${pfx}.pgen" "${pfx}.pvar" "$shared_psam" >> "$merge_list"
    n_chr=$((n_chr + 1))
  done
  [ "$n_chr" -gt 0 ] || { echo "ERROR: no per-chromosome pgen files for ${study}" >&2; return 1; }
  sort -V -o "$merge_list" "$merge_list"
  echo "  chromosomes found: ${n_chr}"

  # 3. Build exclusion lists from review/
  #    Samples: ancestry + related (FID IID); Variants: hwe (variant IDs).
  local samples_remove="${geno}/samples_remove.txt"
  local variants_remove="${geno}/variants_remove.txt"
  cat "${review}/ancestry.exclude" "${review}/related.exclude" 2>/dev/null \
    | awk 'NF' | sort -u > "$samples_remove" || : > "$samples_remove"
  if [ -f "${review}/hwe.exclude" ]; then
    awk 'NF' "${review}/hwe.exclude" | sort -u > "$variants_remove"
  else
    : > "$variants_remove"
  fi
  echo "  samples to remove (ancestry+related): $(wc -l < "$samples_remove" | tr -d ' ')"
  echo "  variants to exclude (hwe):            $(wc -l < "$variants_remove" | tr -d ' ')"

  # 4. Heterozygosity: FLAG ONLY (mean +/- HET_SD * SD of F), never remove.
  #    het.het columns: FID IID O(HOM) E(HOM) N(NM) F  (F is the last column).
  local het_flagged="${out_study}/${study}.het_flagged.txt"
  if [ -f "${review}/het.het" ]; then
    awk -v sd="$HET_SD" '
      NR>1 { f=$NF; n++; sum+=f; sumsq+=f*f; fid[n]=$1; iid[n]=$2; val[n]=f }
      END {
        if (n < 2) { exit 0 }
        mean = sum/n
        var  = (sumsq - n*mean*mean)/(n-1)
        if (var < 0) var = 0
        s = sqrt(var)
        lo = mean - sd*s; hi = mean + sd*s
        for (i=1;i<=n;i++) if (val[i] < lo || val[i] > hi) print fid[i]"\t"iid[i]"\t"val[i]
      }' "${review}/het.het" > "$het_flagged" || : > "$het_flagged"
    echo "  het outliers flagged (NOT removed):   $(wc -l < "$het_flagged" | tr -d ' ')  [mean +/- ${HET_SD} SD of F]"
  else
    : > "$het_flagged"
    echo "  het.het not found; skipping het report"
  fi

  # 5. Merge per-chromosome filesets into one all-chromosome dataset.
  local merged="${geno}/${study}_merged"
  plink2 \
    --pmerge-list "$merge_list" \
    --make-pgen \
    --threads "$THREADS" \
    ${PLINK_MEM_ARGS[@]+"${PLINK_MEM_ARGS[@]}"} \
    --out "$merged"

  # 6. Apply sample + variant exclusions in a single pass -> final output.
  #    Only pass --remove/--exclude when the corresponding list is non-empty,
  #    so an empty review file (e.g. related.exclude) never trips plink2.
  local final="${out_study}/${study}"
  local -a filter_args=(--pfile "$merged" --make-pgen --threads "$THREADS" ${PLINK_MEM_ARGS[@]+"${PLINK_MEM_ARGS[@]}"} --out "$final")
  [ -s "$samples_remove" ]  && filter_args+=(--remove "$samples_remove")
  [ -s "$variants_remove" ] && filter_args+=(--exclude "$variants_remove")
  plink2 "${filter_args[@]}"

  # 7. Stage upstream PCs for reuse as covariates by the GWAS pipeline.
  cp -f "${review}/pca.eigenvec" "${final}.eigenvec"
  cp -f "${review}/pca.eigenval" "${final}.eigenval"

  # 8. Clean the per-study staging area.
  rm -rf "$geno"

  local n_samp n_var
  n_samp=$(($(wc -l < "${final}.psam" | tr -d ' ') - 1))
  n_var=$(wc -l < "${final}.pvar" | tr -d ' ')
  echo "  WROTE: ${final}.{pgen,pvar,psam}  (samples=${n_samp}, variant_lines=${n_var})"
  echo "END: study=${study}"
}

# =============================================================================
# Dispatch
# =============================================================================
if [ "${1:-}" = "--all" ]; then
  for s in "${STUDIES[@]}"; do process_study "$s"; done
elif [ -n "${1:-}" ]; then
  process_study "$1"
elif [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  idx="$SLURM_ARRAY_TASK_ID"
  if [ "$idx" -ge "${#STUDIES[@]}" ]; then
    echo "SKIP: array index ${idx} >= study count ${#STUDIES[@]}; nothing to do."
    exit 0
  fi
  process_study "${STUDIES[$idx]}"
else
  echo "Usage: bash src/002_data-genetics.sh <STUDY> | --all" >&2
  echo "       (or submit as a SLURM array: sbatch src/002_data-genetics.sh)" >&2
  echo "Discovered ${#STUDIES[@]} studies under ${SRC_BASE}:" >&2
  printf '  %s\n' "${STUDIES[@]}" >&2
  exit 1
fi
