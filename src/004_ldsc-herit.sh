#!/bin/bash
#SBATCH --job-name=004_ldsc-herit-setup
#SBATCH --output=src/logs/004_ldsc-herit.log
#SBATCH --error=src/logs/004_ldsc-herit.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=13
#SBATCH --partition=high_p

# src/004_ldsc-herit.sh -- ALL setup for the post-meta analysis modules
# (pipeline/modules/ldsc.nf and pipeline/modules/heritability.nf), including data
# downloading. Kept separate from 001_source.sh because two of its steps depend
# on outputs of earlier stages (so it cannot run at pure "source" time):
#   * the protein-gene map needs meta-feature.txt from 003_data-phenotype.R
#   * the LD reference needs the per-study genotypes from 002_data-genetics
#
# Data layout under ${GWAS_DATA_DIR:-/WHERE/epic-gwas/data}:
#   reference/   SHARED by both modules: LDSC source repo + EUR LD scores +
#                HapMap3 (no-MHC) SNP list.
#   herit/       HERITABILITY-only: protein_gene_map.tsv + combined LD reference.
#   (No ldsc/ dir: the LDSC module uses only the shared reference data.)
#
# The python/plink env (env/herit, env/ldsc) is built by src/001_source.sh.
#
# Usage:
#   bash   src/004_ldsc-herit.sh                 # everything (login node: needs internet)
#   bash   src/004_ldsc-herit.sh --reference     # just the shared downloads
#   bash   src/004_ldsc-herit.sh --protein-map   # just the SeqId->gene map (needs 003 + internet)
#   sbatch src/004_ldsc-herit.sh --ld-ref        # just the LD reference (needs 002 + 003)
#   bash   src/004_ldsc-herit.sh --force         # rebuild/redownload everything
#
# LD reference speed knobs:
#   GWAS_LDREF_JOBS=13 sbatch src/004_ldsc-herit.sh --ld-ref
#   The per-study subsets are cached under data/herit/ldref_parts/ so reruns
#   resume from completed studies instead of starting over.

set -euo pipefail

for env_file in "${SLURM_SUBMIT_DIR:-/nonexistent}/.env" ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env"; do
  if [ -f "$env_file" ]; then set -a; source "$env_file"; set +a; break; fi
done

# Under SLURM the script runs from the spool dir (a copy), so BASH_SOURCE does not
# point at the project. Resolve the root from GWAS_PROJECT_ROOT / SLURM_SUBMIT_DIR
# first (as 001_source.sh does), falling back to the script location.
if [[ -n "${GWAS_PROJECT_ROOT:-}" ]]; then
  PROJ_ROOT="${GWAS_PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  PROJ_ROOT="${SLURM_SUBMIT_DIR}"
else
  PROJ_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
fi
HERIT_ENV="${PROJ_ROOT}/env/herit"

DATA_DIR="${GWAS_DATA_DIR:-/WHERE/epic-gwas/data}"
REFERENCE_DIR="${DATA_DIR}/reference"   # shared (LDSC source + EUR refs + HM3)
HERIT_DIR="${DATA_DIR}/herit"           # heritability-only (gene map + LD reference)

LDSC_HOME="${REFERENCE_DIR}/ldsc"
LDSC_REF_DIR="${REFERENCE_DIR}/eur_w_ld_chr"
LDSC_HM3="${REFERENCE_DIR}/w_hm3.snplist"
LDSC_REPO="https://github.com/CBIIT/ldsc.git"
LDSC_BRANCH="ldsc39"
LDSC_ZENODO="https://zenodo.org/api/records/10515792/files"

FEATURE_FILE="${GWAS_FEATURE_FILE:-${DATA_DIR}/phenotype/meta-feature.txt}"
# Phenofile (003) defines the exact (IID, STUDY) pairs used in the GWAS -- the LD
# reference must use only these, not the full per-study .psam.
PHENO_FILE="${GWAS_PHENO_FILE:-${DATA_DIR}/phenotype/phenofile.txt}"
PROTEIN_MAP="${HERIT_DIR}/protein_gene_map.tsv"
LDREF="${HERIT_DIR}/ldref"
LDREF_MAF="${GWAS_LDREF_MAF:-0.01}"
LDREF_JOBS="${GWAS_LDREF_JOBS:-13}"
GENETICS_DIR="${GWAS_GENETICS_DIR:-}"

FORCE=0
DO_REF=0; DO_MAP=0; DO_LDREF=0; DO_ALL=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference)    DO_REF=1;   DO_ALL=0; shift ;;
    --protein-map)  DO_MAP=1;   DO_ALL=0; shift ;;
    --ld-ref)       DO_LDREF=1; DO_ALL=0; shift ;;
    --force)        FORCE=1; shift ;;
    -h|--help)      sed -n '9,33p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
if [[ "$DO_ALL" -eq 1 ]]; then DO_REF=1; DO_MAP=1; DO_LDREF=1; fi

command -v conda >/dev/null 2>&1 || { echo "ERROR: conda required" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---------------------------------------------------------------------------
# Shared reference data (reference/): LDSC source + EUR LD scores + HM3 list.
# ---------------------------------------------------------------------------
setup_reference() {
  echo ">>> [reference] LDSC source + EUR panel -> ${REFERENCE_DIR}"
  mkdir -p "$REFERENCE_DIR"

  # 1. LDSC source (CBIIT python-3 port, ldsc39).
  if [[ -d "${LDSC_HOME}/.git" && "$FORCE" -ne 1 ]]; then
    echo "    [skip] LDSC source present: ${LDSC_HOME} (use --force to update)"
  elif [[ -d "${LDSC_HOME}/.git" ]]; then
    git -C "$LDSC_HOME" fetch --depth 1 origin "$LDSC_BRANCH" \
      && git -C "$LDSC_HOME" checkout "$LDSC_BRANCH" \
      && echo "    [ok] LDSC source updated" || echo "    [warn] LDSC update failed" >&2
  else
    git clone --depth 1 --branch "$LDSC_BRANCH" "$LDSC_REPO" "$LDSC_HOME" \
      && echo "    [ok] LDSC source: ${LDSC_HOME}" \
      || echo "    [warn] git clone of LDSC failed (needs internet; run on a login node)." >&2
  fi

  # 2. EUR LD scores (Zenodo 10515792): LDscore.<chr>.l2.* -> strip prefix.
  if [[ -f "${LDSC_REF_DIR}/.download_complete" && "$FORCE" -ne 1 ]]; then
    echo "    [skip] EUR LD scores present: ${LDSC_REF_DIR}"
  else
    echo "    [..]   downloading EUR LD scores (Zenodo, ~13 MB)"
    mkdir -p "$LDSC_REF_DIR"
    local tmp_dir; tmp_dir="$(mktemp -d "${REFERENCE_DIR}/ld_XXXXXX")"
    if curl -fsSL "${LDSC_ZENODO}/1000G_Phase3_ldscores.tgz/content" \
         | tar -xzf - -C "$tmp_dir" --strip-components=1; then
      local f base
      for f in "$tmp_dir"/LDscore.*; do
        [[ -f "$f" ]] || continue
        base="$(basename "$f")"; mv -f "$f" "${LDSC_REF_DIR}/${base#LDscore.}"
      done
      touch "${LDSC_REF_DIR}/.download_complete"
      echo "    [ok] EUR LD scores: ${LDSC_REF_DIR}"
    else
      echo "    [warn] EUR LD score download failed (needs internet)." >&2
    fi
    rm -rf "$tmp_dir"
  fi

  # 3. HapMap3 SNP list -- the no-MHC variant (so LDSC excludes the MHC).
  if [[ -f "$LDSC_HM3" && "$FORCE" -ne 1 ]]; then
    echo "    [skip] HapMap3 (no-MHC) list present: ${LDSC_HM3}"
  else
    echo "    [..]   downloading HapMap3 (no-MHC) SNP list (Zenodo, ~12 MB)"
    if { printf 'SNP\n'; curl -fsSL "${LDSC_ZENODO}/hm3_no_MHC.list.txt/content"; } > "${LDSC_HM3}.tmp"; then
      mv -f "${LDSC_HM3}.tmp" "$LDSC_HM3"
      echo "    [ok] HapMap3 list: ${LDSC_HM3} ($(wc -l < "$LDSC_HM3" | tr -d ' ') lines)"
    else
      rm -f "${LDSC_HM3}.tmp"
      echo "    [warn] HapMap3 list download failed (needs internet)." >&2
    fi
  fi
}

# ---------------------------------------------------------------------------
# Heritability-only data (herit/): SeqId->gene hg38 map.
# ---------------------------------------------------------------------------
setup_protein_map() {
  echo ">>> [herit] SeqId -> gene hg38 map -> ${PROTEIN_MAP}"
  if [[ ! -d "$HERIT_ENV" ]]; then
    echo "    [skip] env/herit not built (run: bash src/001_source.sh)" >&2; return 0
  fi
  if [[ ! -f "$FEATURE_FILE" ]]; then
    echo "    [skip] feature file not found: ${FEATURE_FILE}"
    echo "           (written by 003_data-phenotype.R; run that first, then re-run --protein-map)"
    return 0
  fi
  if [[ -f "$PROTEIN_MAP" && "$FORCE" -ne 1 ]]; then
    echo "    [skip] protein map present: ${PROTEIN_MAP} (use --force)"; return 0
  fi
  mkdir -p "$HERIT_DIR"
  conda activate "$HERIT_ENV"
  if python "${PROJ_ROOT}/pipeline/bin/build_protein_gene_map.py" \
       --feature-file "$FEATURE_FILE" --out "$PROTEIN_MAP"; then
    echo "    [ok] ${PROTEIN_MAP}"
  else
    echo "    [warn] protein map build failed -- it queries MyGene.info, so run on a" >&2
    echo "           node with internet (login node)." >&2
  fi
  conda deactivate || true
}

# ---------------------------------------------------------------------------
# Heritability-only data (herit/): LD reference (bed/bim/fam, IDs = chr:pos) over
# exactly the GWAS samples, for PLINK --clump.
# ---------------------------------------------------------------------------
setup_ld_reference() {
  echo ">>> [herit] GWAS-sample LD reference -> ${LDREF}.{bed,bim,fam}"
  if [[ ! -d "$HERIT_ENV" ]]; then
    echo "    [skip] env/herit not built (run: bash src/001_source.sh)" >&2; return 0
  fi
  if [[ -z "$GENETICS_DIR" ]]; then
    echo "    [skip] GWAS_GENETICS_DIR not set; cannot build LD reference."; return 0
  fi
  if [[ ! -f "$PHENO_FILE" ]]; then
    echo "    [skip] phenofile not found: ${PHENO_FILE}"
    echo "           (written by 003_data-phenotype.R; it defines the GWAS study-sample pairs)"
    return 0
  fi
  if [[ -f "${LDREF}.bed" && "$FORCE" -ne 1 ]]; then
    echo "    [skip] LD reference present: ${LDREF}.bed (use --force)"; return 0
  fi

  mkdir -p "$HERIT_DIR"
  conda activate "$HERIT_ENV"
  if ! command -v plink2 >/dev/null 2>&1; then
    echo "    [warn] plink2 not found after activating ${HERIT_ENV}; LD reference not created." >&2
    conda deactivate || true
    return 0
  fi

  if ! [[ "$LDREF_JOBS" =~ ^[0-9]+$ ]] || [[ "$LDREF_JOBS" -lt 1 ]]; then
    echo "    [warn] invalid GWAS_LDREF_JOBS=${LDREF_JOBS}; using 1" >&2
    LDREF_JOBS=1
  fi

  local threads_total="${SLURM_CPUS_PER_TASK:-8}"
  if [[ "$LDREF_JOBS" -gt "$threads_total" ]]; then
    LDREF_JOBS="$threads_total"
  fi
  local threads_per_job=$(( threads_total / LDREF_JOBS ))
  [[ "$threads_per_job" -ge 1 ]] || threads_per_job=1

  local -a subset_mem_args=()
  local -a merge_mem_args=()
  if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" -gt 2048 ]] 2>/dev/null; then
    local subset_mem=$(( (SLURM_MEM_PER_NODE - 1024) / LDREF_JOBS ))
    [[ "$subset_mem" -ge 1024 ]] || subset_mem=1024
    subset_mem_args=(--memory "$subset_mem")
    merge_mem_args=(--memory $(( SLURM_MEM_PER_NODE - 1024 )))
  fi

  echo "    LD reference build jobs: ${LDREF_JOBS}; threads/job: ${threads_per_job}; merge threads: ${threads_total}"

  local keep_dir="${HERIT_DIR}/ldref_keep"
  local parts_dir="${HERIT_DIR}/ldref_parts"
  local logs_dir="${parts_dir}/logs"
  rm -rf "$keep_dir"
  mkdir -p "$keep_dir" "$parts_dir" "$logs_dir"

  # The GWAS analysis set is exactly the (IID, STUDY) pairs in the phenofile (003:
  # proteomics-matched, each person assigned to ONE retained study, under-powered
  # studies excluded) -- NOT the raw .psam, which holds many more samples. Write a
  # two-column plink2 keep file per study from it. The staged .psam files use
  # FID == IID, matching PREPARE_PHENOTYPES/LD_PRUNING keep files. Each IID is in
  # exactly one study, so there are no duplicate IIDs at merge time.
  if ! awk -F'\t' -v d="$keep_dir" '
      NR==1{ for(i=1;i<=NF;i++){ if($i=="IID") ic=i; if($i=="STUDY") sc=i }
             if(!ic||!sc) exit 3; next }
      $ic != "" && $sc != "" {
        key=$sc SUBSEP $ic
        if (!(key in seen)) {
          f=d"/"$sc".keep"
          print $ic "\t" $ic >> f
          seen[key]=1
        }
      }
    ' "$PHENO_FILE"; then
    echo "    [warn] phenofile has no IID/STUDY columns: ${PHENO_FILE}" >&2
    rm -rf "$keep_dir"; conda deactivate || true; return 0
  fi

  build_study_ref() {
    local keep="$1"
    local study pfile_prefix n_keep part_prefix plink_log ext
    study="$(basename "$keep" .keep)"
    pfile_prefix="${GENETICS_DIR}/${study}/${study}"
    part_prefix="${parts_dir}/${study}_ref"
    plink_log="${logs_dir}/${study}.plink2.log"
    n_keep=$(wc -l < "$keep" | tr -d ' ')

    for ext in pgen pvar psam; do
      if [[ ! -f "${pfile_prefix}.${ext}" ]]; then
        echo "    [warn] missing genotype file for study ${study}: ${pfile_prefix}.${ext}; skipping" >&2
        return 0
      fi
    done

    if [[ -f "${part_prefix}.done" && -f "${part_prefix}.pgen" && -f "${part_prefix}.pvar" && -f "${part_prefix}.psam" && "$FORCE" -ne 1 ]]; then
      echo "    [skip] ${study}: cached subset (${n_keep} GWAS samples)"
      return 0
    fi

    rm -f "${part_prefix}.pgen" "${part_prefix}.pvar" "${part_prefix}.psam" \
          "${part_prefix}.log" "${part_prefix}.done"
    echo "    ${study}: ${n_keep} GWAS samples"
    if plink2 --pfile "$pfile_prefix" \
         --keep "$keep" \
         --set-all-var-ids '@:#' --new-id-max-allele-len 200 missing \
         --rm-dup exclude-all \
         --maf "$LDREF_MAF" \
         --threads "$threads_per_job" ${subset_mem_args[@]+"${subset_mem_args[@]}"} \
         --make-pgen --out "$part_prefix" > "$plink_log" 2>&1; then
      touch "${part_prefix}.done"
    else
      echo "    [warn] plink2 subset failed for ${study}; skipping. Last log lines:" >&2
      tail -n 20 "$plink_log" | sed 's/^/      /' >&2
    fi
  }

  local keep active
  for keep in "${keep_dir}"/*.keep; do
    [[ -e "$keep" ]] || { echo "    [warn] no studies parsed from phenofile" >&2; break; }
    while true; do
      active=$(jobs -rp | wc -l | tr -d ' ')
      [[ "$active" -lt "$LDREF_JOBS" ]] && break
      sleep 5
    done
    build_study_ref "$keep" &
  done
  while [[ "$(jobs -rp | wc -l | tr -d ' ')" -gt 0 ]]; do
    sleep 5
  done
  wait || true

  local tmp; tmp="$(mktemp -d "${HERIT_DIR}/build_XXXXXX")"
  local refs_list="${tmp}/merge_list.txt"; : > "$refs_list"
  local study part_prefix
  for keep in "${keep_dir}"/*.keep; do
    [[ -e "$keep" ]] || continue
    study="$(basename "$keep" .keep)"
    part_prefix="${parts_dir}/${study}_ref"
    if [[ -f "${part_prefix}.done" && -f "${part_prefix}.pgen" && -f "${part_prefix}.pvar" && -f "${part_prefix}.psam" ]]; then
      echo "$part_prefix" >> "$refs_list"
    fi
  done

  local n_refs; n_refs=$(wc -l < "$refs_list" | tr -d ' ')
  if [[ "$n_refs" -lt 1 ]]; then
    echo "    [warn] no per-study reference filesets built; LD reference not created." >&2
    rm -rf "$tmp" "$keep_dir"; conda deactivate || true; return 0
  fi
  rm -f "${LDREF}.bed" "${LDREF}.bim" "${LDREF}.fam" "${LDREF}.log" "${LDREF}.nosex"

  # PLINK2 cannot currently do the sample-wise merge we need here when variants
  # overlap across filesets ("Non-concatenating --pmerge-list is under
  # development"). Convert each cached PGEN subset to BED, then use PLINK 1.9's
  # mature --merge-list implementation to combine samples across studies.
  if ! command -v plink >/dev/null 2>&1; then
    echo "    [warn] plink 1.9 not found after activating ${HERIT_ENV}; LD reference not merged." >&2
    rm -rf "$tmp" "$keep_dir"; conda deactivate || true; return 0
  fi

  local bed_prefixes="${tmp}/bed_prefixes.txt"; : > "$bed_prefixes"
  local part_prefix part_name bed_log
  while IFS= read -r part_prefix; do
    [[ -n "$part_prefix" ]] || continue
    part_name="$(basename "$part_prefix")"
    bed_log="${logs_dir}/${part_name}.makebed.log"
    if [[ -f "${part_prefix}.bed.done" && -f "${part_prefix}.bed" && -f "${part_prefix}.bim" && -f "${part_prefix}.fam" && "$FORCE" -ne 1 ]]; then
      echo "    [skip] ${part_name}: cached BED conversion"
    else
      rm -f "${part_prefix}.bed" "${part_prefix}.bim" "${part_prefix}.fam" "${part_prefix}.bed.done"
      echo "    converting ${part_name} to BED"
      if plink2 --pfile "$part_prefix" \
          --threads "$threads_per_job" ${subset_mem_args[@]+"${subset_mem_args[@]}"} \
          --make-bed --out "$part_prefix" > "$bed_log" 2>&1; then
        touch "${part_prefix}.bed.done"
      else
        echo "    [warn] BED conversion failed for ${part_name}; Last log lines:" >&2
        tail -n 20 "$bed_log" | sed 's/^/      /' >&2
        continue
      fi
    fi
    echo "$part_prefix" >> "$bed_prefixes"
  done < "$refs_list"

  n_refs=$(wc -l < "$bed_prefixes" | tr -d ' ')
  if [[ "$n_refs" -lt 1 ]]; then
    echo "    [warn] no BED filesets available after conversion; LD reference not created." >&2
    rm -rf "$tmp" "$keep_dir"; conda deactivate || true; return 0
  fi

  if [[ "$n_refs" -eq 1 ]]; then
    plink --bfile "$(cat "$bed_prefixes")" --maf "$LDREF_MAF" \
      --threads "$threads_total" ${merge_mem_args[@]+"${merge_mem_args[@]}"} \
      --make-bed --out "$LDREF"
  else
    local base_prefix merge_bed_list missnp
    base_prefix="$(head -n 1 "$bed_prefixes")"
    merge_bed_list="${tmp}/merge_bed_list.txt"
    tail -n +2 "$bed_prefixes" | awk '{ print $0 ".bed\t" $0 ".bim\t" $0 ".fam" }' > "$merge_bed_list"

    echo "    merging ${n_refs} cached per-study BED references"
    if ! plink --bfile "$base_prefix" \
        --merge-list "$merge_bed_list" \
        --maf "$LDREF_MAF" \
        --threads "$threads_total" ${merge_mem_args[@]+"${merge_mem_args[@]}"} \
        --make-bed --out "$LDREF"; then
      local missnp="${LDREF}-merge.missnp"
      if [[ ! -s "$missnp" ]]; then
        echo "    [warn] PLINK 1.9 merge failed and no ${missnp} file was produced." >&2
        rm -rf "$tmp" "$keep_dir"; conda deactivate || true; return 0
      fi
      local n_excl; n_excl=$(wc -l < "$missnp" | tr -d ' ')
      echo "    [warn] multiallelic variants across studies: ${n_excl} — excluding from all studies and retrying" >&2
      rm -f "${LDREF}.bed" "${LDREF}.bim" "${LDREF}.fam" "${LDREF}.log" "${LDREF}.nosex"

      # --exclude on the primary file only (plink1.9 default) is not enough: the
      # multiallelic variants are spread across ALL study files. Filter every BED.
      local filtered_list="${tmp}/filtered_prefixes.txt"; : > "$filtered_list"
      while IFS= read -r pfx; do
        [[ -n "$pfx" ]] || continue
        local pfx_name pfx_filt filt_log
        pfx_name="$(basename "$pfx")"
        pfx_filt="${tmp}/${pfx_name}_filt"
        filt_log="${logs_dir}/${pfx_name}.filt.log"
        if plink --bfile "$pfx" --exclude "$missnp" \
            --threads "$threads_per_job" ${subset_mem_args[@]+"${subset_mem_args[@]}"} \
            --make-bed --out "$pfx_filt" > "$filt_log" 2>&1; then
          echo "$pfx_filt" >> "$filtered_list"
        else
          echo "    [warn] missnp filter failed for ${pfx_name}; see ${filt_log}" >&2
        fi
      done < "$bed_prefixes"

      local n_filt; n_filt=$(wc -l < "$filtered_list" | tr -d ' ')
      if [[ "$n_filt" -lt 1 ]]; then
        echo "    [warn] no filtered BED sets available; LD reference not created." >&2
        rm -rf "$tmp" "$keep_dir"; conda deactivate || true; return 0
      fi

      local filt_base filt_merge_list
      filt_base="$(head -n 1 "$filtered_list")"
      filt_merge_list="${tmp}/merge_filt_list.txt"
      tail -n +2 "$filtered_list" | awk '{ print $0 ".bed\t" $0 ".bim\t" $0 ".fam" }' > "$filt_merge_list"

      if [[ "$n_filt" -eq 1 ]]; then
        plink --bfile "$filt_base" --maf "$LDREF_MAF" \
          --threads "$threads_total" ${merge_mem_args[@]+"${merge_mem_args[@]}"} \
          --make-bed --out "$LDREF"
      else
        plink --bfile "$filt_base" \
          --merge-list "$filt_merge_list" \
          --maf "$LDREF_MAF" \
          --threads "$threads_total" ${merge_mem_args[@]+"${merge_mem_args[@]}"} \
          --make-bed --out "$LDREF"
      fi
    fi
  fi
  echo "    [ok] ${LDREF}.{bed,bim,fam}  variants: $(wc -l < "${LDREF}.bim" | tr -d ' ')  samples: $(wc -l < "${LDREF}.fam" | tr -d ' ')"
  rm -rf "$tmp" "$keep_dir"
  conda deactivate || true
}

[[ "$DO_REF"   -eq 1 ]] && setup_reference
[[ "$DO_MAP"   -eq 1 ]] && setup_protein_map
[[ "$DO_LDREF" -eq 1 ]] && setup_ld_reference

echo ">>> 004 LDSC/heritability setup complete."
