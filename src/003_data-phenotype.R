# Script: src/003_data-phenotype.R
# Purpose: Prepare phenotype inputs for the GWAS pipeline.

# environment ====
load_dot_env <- function(env_file = ".env") {
  if (file.exists(env_file)) {
    lines <- readLines(env_file, warn = FALSE)
    lines <- lines[!grepl("^\\s*(#|$)", lines)]
    for (line in lines) {
      parts <- strsplit(line, "=")[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = "="))
        value <- gsub("^['\"]|['\"]$", "", value)
        env_list <- list(value)
        names(env_list) <- key
        do.call(Sys.setenv, env_list)
      }
    }
  }
}

project_root <- Sys.getenv("GWAS_PROJECT_ROOT")
search_paths <- c(
  if (project_root != "") project_root else NULL,
  ".",
  "..",
  "../.."
)

for (p in search_paths) {
  env_p <- file.path(p, ".env")
  if (file.exists(env_p)) {
    load_dot_env(env_p)
    break
  }
}

get_env <- function(key, default = NULL) {
  val <- Sys.getenv(key)
  if (val == "") return(default)
  return(val)
}


# data somalogic ====
data_soma_path <- get_env("GWAS_EXTERNAL_SOMA_RDS")
if (is.null(data_soma_path) || !file.exists(data_soma_path)) {
  stop("GWAS_EXTERNAL_SOMA_RDS not found. Check your .env file.")
}
data_soma <- readRDS(data_soma_path)
data_features <- data_soma$data
meta_sample <- data_soma$others
meta_feature <- data_soma$aux

# format
# Idepic       : soma/phenotype identifier (single-underscore form), key for proteins
# Idepic_Bio   : genotype identifier (matches the PLINK2 .psam IID), key for genetics
meta_sample <- meta_sample |>
  dplyr::mutate(IID = Idepic,
                IID_clean = gsub("_+", "_", Idepic),
                bio = Idepic_Bio)
data_features <- data_features |>
  dplyr::left_join(meta_sample |>
                     dplyr::select(SubjectID, IID, IID_clean),
                   by = "SubjectID") |>
  dplyr::select(SubjectID, IID, IID_clean, dplyr::everything()) |>
  dplyr::rename_with(
    ~ gsub("\\.", "-", gsub("^seq\\.", "", .x)),
    dplyr::starts_with("seq."))

protein_cols <- unique(meta_feature$SeqId)

# data genetic ====
# Read the analysis-ready, merged per-study PLINK2 sample files produced by
# 002_data-genetics.sh (one <study>.psam per study, exclusions already applied).
genetics_path <- get_env("GWAS_GENETICS_DIR")
if (is.null(genetics_path) || !dir.exists(genetics_path)) {
  stop("GWAS_GENETICS_DIR not found. Run 002_data-genetics.sh first and check .env.")
}

# Genotype IIDs are the Idepic_Bio form. Some are compound, e.g.
# "00____0000000A_00____00000001", where the soma-matching piece is the pure-digit
# "<digits>____<digits>" component. extract_bio() returns the Idepic_Bio that a
# genotype IID maps to (the IID itself if already a bio ID, else the matching
# component), or NA when the genotyped sample has no proteomics.
bio_set <- unique(stats::na.omit(meta_sample$bio))
bio_lookup <- meta_sample |>
  dplyr::distinct(bio, IID_clean)

extract_bio <- function(iid, bio_set) {
  hit <- ifelse(iid %in% bio_set, iid, NA_character_)
  need <- is.na(hit)
  if (any(need)) {
    cand <- regmatches(iid[need], gregexpr("[0-9]+_{4,}[0-9]+", iid[need]))
    hit[need] <- vapply(cand, function(ms) {
      m <- ms[ms %in% bio_set]
      if (length(m)) m[[1]] else NA_character_
    }, character(1))
  }
  hit
}

data_genetic <- Sys.glob(file.path(genetics_path, "*", "*.psam")) |>
  purrr::map_dfr(function(f) {
    STUDY_name <- basename(f) |> tools::file_path_sans_ext()
    dt <- data.table::fread(f)
    data.table::setnames(dt, sub("^#", "", names(dt)))
    # IID_geno is the genotype key that REGENIE must match. Also capture the
    # case/control disease status (PHENO / PHENO1): it is used as a covariate in
    # the combined GWAS, so samples missing it are dropped by REGENIE (handled by
    # the missing-PHENO exclusion below). has_pheno_col records whether this
    # study's .psam carries a PHENO column at all, so studies without one are
    # never treated as "missing".
    pheno_col <- intersect(c("PHENO", "PHENO1"), names(dt))
    data.frame(
      IID_geno      = as.character(dt$IID),
      STUDY         = STUDY_name,
      PHENO_geno    = if (length(pheno_col) > 0) as.character(dt[[pheno_col[1]]]) else NA_character_,
      has_pheno_col = length(pheno_col) > 0,
      stringsAsFactors = FALSE
    )
  }) |>
  dplyr::mutate(bio = extract_bio(IID_geno, bio_set)) |>
  dplyr::filter(!is.na(bio)) |>            # keep genotyped samples with proteomics
  dplyr::left_join(bio_lookup, by = "bio") |>   # adds IID_clean (Idepic) for joins
  dplyr::distinct(IID_geno, STUDY, .keep_all = TRUE)

# comprehensive sample completeness assessment (pre-dedup) ====
# Every sample entering the GWAS must have complete data for ALL covariates the
# association model conditions on. REGENIE removes any sample with a missing
# covariate value from every test, silently lowering the effective N (observed:
# 9,203 -> 9,174 from missing PHENO alone). We therefore assess all GWAS
# covariates up front -- before any dedup or study-retention counting -- and
# remove incomplete samples so they never influence study assignment or sizes,
# and every downstream N is consistent and explicit (keep_samples.txt, PCA,
# covariates, REGENIE).
#
# GWAS covariates assessed:
#   - PHENO    : case/control disease status (.psam; combined-analysis
#                adjustment), only where the study .psam carries a PHENO column.
#   - Sex, Age_Recr : from the SomaLogic sample metadata.
# Proteins are OUTCOMES (per-trait missingness is handled by REGENIE, so a sample
# missing one protein still contributes to the others) and PlateId is excluded
# from the covariate set, so neither triggers sample removal. PCs are computed
# downstream on the retained set and are always complete.

# Join the metadata covariates needed for the check (reused for the phenofile).
data_genetic <- data_genetic |>
  dplyr::left_join(
    meta_sample |>
      dplyr::distinct(IID_clean, .keep_all = TRUE) |>
      dplyr::select(IID_clean, Sex, Age_Recr, PlateId),
    by = "IID_clean"
  )

# Robust to haven_labelled (SPSS) columns from the soma RDS, whose as.character
# method errors: strip the class to the underlying atomic vector first.
is_missing <- function(x) {
  x <- unclass(x)
  attributes(x) <- NULL
  is.na(x) | trimws(as.character(x)) %in% c("", "NA")
}

miss_flags <- list(
  Sex      = is_missing(data_genetic$Sex),
  Age_Recr = is_missing(data_genetic$Age_Recr),
  # PHENO is only a covariate (and only "missing") where the study uses it.
  PHENO    = data_genetic$has_pheno_col & is_missing(data_genetic$PHENO_geno)
)

incomplete <- Reduce(`|`, miss_flags)

message("Sample completeness over GWAS covariates (Sex, Age_Recr, PHENO):")
for (v in names(miss_flags)) {
  message("  missing ", v, ": ", sum(miss_flags[[v]]))
}

if (any(incomplete)) {
  excl_tbl <- data_genetic[incomplete, , drop = FALSE] |>
    dplyr::count(STUDY, name = "n_excluded")
  message("Excluding ", sum(incomplete),
          " samples with any missing GWAS covariate, before dedup:")
  for (i in seq_len(nrow(excl_tbl))) {
    message("  ", excl_tbl$STUDY[i], ": ", excl_tbl$n_excluded[i])
  }
  data_genetic <- data_genetic[!incomplete, , drop = FALSE]
} else {
  message("No samples with missing GWAS covariates to exclude.")
}

# Deduplicate across studies and exclude under-powered studies.
#   1. A provisional "keep in smallest" assignment gives each study its maximum
#      reachable matched size; studies still below min_study_n are excluded.
#   2. Each person is then assigned to the SMALLEST RETAINED study they appear in
#      (duplicates reassigned to smaller studies, but never into a study that would
#      be excluded). People present only in excluded studies are dropped.
min_study_n <- as.integer(get_env("GWAS_MIN_STUDY_N", "100"))

prov_size <- data_genetic |> dplyr::count(STUDY, name = "study_n")

retained_studies <- data_genetic |>
  dplyr::left_join(prov_size, by = "STUDY") |>
  dplyr::group_by(bio) |>
  dplyr::arrange(study_n, STUDY, .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::count(STUDY, name = "n") |>
  dplyr::filter(n >= min_study_n) |>
  dplyr::pull(STUDY)

message("Retained ", length(retained_studies), " studies with >= ", min_study_n,
        " matched samples; excluded ",
        dplyr::n_distinct(data_genetic$STUDY) - length(retained_studies), ".")

retained_rank <- prov_size |>
  dplyr::filter(STUDY %in% retained_studies) |>
  dplyr::arrange(study_n, STUDY) |>
  dplyr::mutate(rk = dplyr::row_number())

data_genetic <- data_genetic |>
  dplyr::filter(STUDY %in% retained_studies) |>
  dplyr::left_join(dplyr::select(retained_rank, STUDY, rk), by = "STUDY") |>
  dplyr::group_by(bio) |>
  dplyr::arrange(rk, .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(-rk)

# make phenofile ====
# The phenofile IID must be the genotype IID (IID_geno) so REGENIE can match it to
# the PLINK2 .psam. The Sex / Age_Recr / PlateId covariates were already joined
# above for the completeness check; here we only add the protein traits (by
# IID_clean / Idepic). PHENO_geno / has_pheno_col are diagnostic only and are
# dropped before writing (PHENO must not leak into the phenofile, where numeric
# columns become traits).
data_pheno <- data_genetic |>
  dplyr::left_join(data_features |>
                     dplyr::select(IID_clean, dplyr::all_of(protein_cols)),
                   by = "IID_clean") |>
  dplyr::mutate(IID = IID_geno) |>
  dplyr::select(IID, STUDY, dplyr::everything(),
                -IID_geno, -IID_clean, -bio, -PHENO_geno, -has_pheno_col)

# write ====
# Phenotype outputs live alongside the genetics data, under <data>/phenotype/,
# where <data> is the parent of GWAS_GENETICS_DIR (genetics_path).
pheno_dir <- file.path(dirname(genetics_path), "phenotype")
if (!dir.exists(pheno_dir)) {
  dir.create(pheno_dir, recursive = TRUE)
}

phenofile <- file.path(pheno_dir, "phenofile.txt")
data.table::fwrite(
  x = data_pheno,
  file = phenofile,
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
message("Wrote ", phenofile, " with ", length(protein_cols), " proteins.")

path <- file.path(pheno_dir, "meta-feature.txt")
data.table::fwrite(
  x = meta_feature,
  file = path,
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

id <- c(
  "15602-43", # IL6R
  "4139-71", # IL6R
  "8092-29", # IL6R
  "5231-79", # PCSK9
  "8246-9", # PCSK9
  "10620-21", # MSMB
  "3580-25", # SERPINA1
  "4337-49", # CRP
  "4721-54", # TFF3
  "8323-163", # TFF3
  "18878-15" # GREM1
)

writeLines(id, con = file.path(pheno_dir, "test.txt"))
message("Wrote test.txt with 10 traits.")

writeLines(protein_cols[1:100], con = file.path(pheno_dir, "test-100.txt"))
message("Wrote test-100.txt with 100 traits.")

n_batches <- 6
batch_idx <- cut(seq_along(protein_cols), n_batches, labels = FALSE)
batches <- split(protein_cols, batch_idx)
for (i in seq_len(n_batches)) {
  filename <- file.path(pheno_dir, paste0("batch", i, ".txt"))
  writeLines(batches[[i]], con = filename)
  n <- length(batches[[i]])
  message("Wrote batch", i, ".txt with ", n, " traits.")
}

# samplesheet ====
# (merged from the former src/004_samplesheet.R) Build the Nextflow samplesheet
# directly from the retained studies and their merged PLINK2 filesets. SomaLogic
# proteins are quantitative, so analysis is normally 'combined' only; a
# cases/controls split is added only if a study's .psam carries a PHENO column
# with both arms large enough (>= GWAS_MIN_ARM_N), counted over the deduped
# samples kept for this study (data_pheno, in memory -- no phenofile re-read).
samplesheet_out <- get_env("GWAS_SAMPLESHEET", "pipeline/samplesheet.csv")
group_col       <- "PHENO"
cases_val       <- "2"
min_arm         <- as.integer(get_env("GWAS_MIN_ARM_N", "100"))
pheno_family    <- c("PHENO", "PHENO1", "PHENOTYPE", "PHENOTYPE1")

read_sample_header <- function(sample_file) {
  hdr_lines <- readLines(sample_file, warn = FALSE)
  hdr_lines <- hdr_lines[!grepl("^\\s*$", hdr_lines)]
  hdr_lines <- hdr_lines[!grepl("^##", hdr_lines)]
  if (length(hdr_lines) == 0) return(character())
  strsplit(trimws(gsub("^#", "", hdr_lines[1])), "\\s+")[[1]]
}

get_abs_path <- function(p) {
  if (grepl("^/", p)) return(p)
  normalizePath(file.path(getwd(), p), mustWork = FALSE)
}

samplesheet_rows <- list()
for (study in sort(unique(data_pheno$STUDY))) {
  pfile_prefix <- file.path(genetics_path, study, study)
  for (ext in c(".pgen", ".pvar", ".psam")) {
    if (!file.exists(paste0(pfile_prefix, ext))) {
      stop("Missing PLINK2 file ", paste0(study, ext), " under ", genetics_path)
    }
  }
  sample_file <- paste0(pfile_prefix, ".psam")
  
  group_col_use  <- group_col
  pheno_col_psam <- NA_character_
  hdr_cols <- read_sample_header(sample_file)
  if (group_col_use %in% hdr_cols) {
    pheno_col_psam <- group_col_use
  } else if (identical(group_col_use, "PHENO") && "PHENO1" %in% hdr_cols) {
    pheno_col_psam <- "PHENO1"
  } else {
    group_col_use <- ""
  }
  
  groups <- "combined"
  if (group_col_use != "" && !is.na(pheno_col_psam)) {
    psam_dt <- data.table::fread(sample_file)
    data.table::setnames(psam_dt, sub("^#", "", names(psam_dt)))
    kept_iids <- data_pheno$IID[data_pheno$STUDY == study]
    ph <- as.character(psam_dt[IID %in% kept_iids][[pheno_col_psam]])
    missing_codes <- c(cases_val, "0", "-9", "NA", "")
    n_cases    <- sum(ph == cases_val, na.rm = TRUE)
    n_controls <- sum(!ph %in% missing_codes & !is.na(ph))
    if (n_cases    >= min_arm && n_controls >= 1) groups <- c(groups, "cases")
    if (n_controls >= min_arm && n_cases    >= 1) groups <- c(groups, "controls")
    cat(sprintf("  %s: cases=%d controls=%d -> groups=%s\n",
                study, n_cases, n_controls, paste(groups, collapse = ",")))
  }
  
  group_column_out <- if (is.na(pheno_col_psam)) {
    ""
  } else if (toupper(pheno_col_psam) %in% pheno_family) {
    "PHENO"
  } else {
    pheno_col_psam
  }
  
  samplesheet_rows[[study]] <- data.frame(
    study_id = study,
    pfile = get_abs_path(pfile_prefix),
    sample_file = get_abs_path(sample_file),
    group_column = group_column_out,
    cases_value = cases_val,
    analysis_groups = paste(groups, collapse = ","),
    stringsAsFactors = FALSE
  )
}

final_df <- do.call(rbind, samplesheet_rows)
rownames(final_df) <- NULL
ss_dir <- dirname(samplesheet_out)
if (!dir.exists(ss_dir)) dir.create(ss_dir, recursive = TRUE)
data.table::fwrite(final_df, file = samplesheet_out, sep = ",", quote = FALSE, row.names = FALSE)
message("Wrote samplesheet: ", samplesheet_out, " (", nrow(final_df), " studies).")
