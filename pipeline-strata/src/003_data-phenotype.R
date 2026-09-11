# Script: src/003_data-phenotype.R  (pipeline-strata)
# Purpose: Prepare phenotype inputs for the STRATIFIED extended analysis.
#
# Differs from the frozen publication tree in four ways:
#   1. Derives the analysis strata columns (SMOKE_CAT, BMI_CAT) from the raw
#      EPIC metadata. The recoding lives in ONE place -- strata_spec below.
#   2. Writes to GWAS_PHENO_DIR, an explicit path, never derived from the
#      shared genetics directory (which would overwrite the published
#      analysis's own phenofile).
#   3. Emits a samplesheet whose analysis_groups are strata names, sized
#      against the strata defined in pipeline/params.yaml.
#   4. Does NOT treat a missing strata value as a reason to drop a sample: an
#      unmapped smoker still contributes to the unstratified analysis.

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


# strata definition ====
# ---------------------------------------------------------------------------
# THIS IS THE STAGE AT WHICH ANALYSIS VARIABLES ARE DEFINED.
#
# Each entry derives ONE new column in the phenofile from the raw EPIC sample
# metadata (data_soma$others). pipeline/params.yaml then refers to that column
# by name and picks which of its levels form each stratum.
#
# To add a variable for a new stratification:
#   1. Add an entry here: the output column name, the raw source column, and a
#      recode() mapping raw values onto the labels you want.
#   2. Add matching strata to pipeline/params.yaml, e.g.
#        alcohol_hi: { column: "ALC_CAT", values: ["HIGH"] }
#   3. Re-run this script; it rewrites the samplesheet with the strata each
#      study is large enough to run.
#
# Rules that make this safe:
#   * recode() returning NA means "in no stratum". That is how exclusions are
#     expressed -- there is no catch-all group to fall into.
#   * optional = TRUE lets a variable be declared before its source column is
#     confirmed: if the column is absent the variable is skipped with a
#     warning instead of aborting the run.
#   * A raw-to-mapped crosstab is printed for every variable, and the script
#     stops if a required variable maps nothing.
# ---------------------------------------------------------------------------
# Flatten a metadata column to plain character, RESOLVING SPSS VALUE LABELS.
#
# The soma RDS comes from SPSS, so columns like smoke_statVV are haven_labelled:
# the storage is integer codes and the labels live in attr(x, "labels"). This
# matters because table(df$smoke_statVV) in an interactive session prints the
# LABELS ("Never", "Former", ...) while the underlying vector holds 1, 2, 3 --
# so a recode written against what you saw on screen silently matches nothing.
# Simply unclassing (as the completeness check does, to dodge haven's erroring
# as.character method) throws the labels away and leaves the bare codes.
#
# Here we translate codes back to their labels when a label dictionary exists,
# and fall back to the code as a string when it does not. Recodes in
# strata_spec can therefore be written against the human-readable levels.
as_plain <- function(x) {
  # FACTORS FIRST. as.character() on a factor yields the level LABELS;
  # unclass() yields the integer level codes. The soma metadata stores
  # smoke_statVV as a factor, so unclassing it turned Never/Former/Smoker into
  # 1/2/3 and every recode matched nothing. table() shows a factor's labels,
  # which is why the coding looked like strings in an interactive session.
  if (is.factor(x)) return(trimws(as.character(x)))

  # haven_labelled (SPSS): storage is integer codes, labels in attr "labels".
  labels <- attr(x, "labels", exact = TRUE)

  raw <- unclass(x)
  attributes(raw) <- NULL
  out <- trimws(as.character(raw))

  if (!is.null(labels) && length(labels) > 0 && !is.null(names(labels))) {
    lookup <- stats::setNames(trimws(names(labels)),
                              trimws(as.character(unname(labels))))
    hit <- unname(lookup[out])
    out <- ifelse(is.na(hit), out, hit)
  }
  out
}

# Human-readable description of a column's levels, for the run log. When there
# is no dictionary at all this prints the observed raw values, so a coding
# mismatch is diagnosable from this one line without cross-referencing the
# crosstab in the other log file.
describe_labels <- function(x) {
  if (is.factor(x)) {
    return(paste0("factor with levels: ", paste(levels(x), collapse = ", ")))
  }
  labels <- attr(x, "labels", exact = TRUE)
  if (!is.null(labels) && length(labels) > 0 && !is.null(names(labels))) {
    return(paste0("SPSS value labels: ",
                  paste(sprintf("%s=%s", unname(labels), names(labels)), collapse = ", ")))
  }
  plain <- unclass(x)
  attributes(plain) <- NULL
  observed <- sort(unique(trimws(as.character(plain))))
  paste0("no label dictionary; observed values: ",
         paste(utils::head(observed, 12), collapse = ", "),
         if (length(observed) > 12) paste0(" ... (", length(observed), " distinct)") else "")
}

strata_spec <- list(

  # Smoking. Verified against the data:
  #   table(data_soma$others$smoke_statVV)
  #     Never 8046 | Former 4911 | Smoker 4535 | Unknown 102
  # Ever = Former + Smoker. Unknown (102) maps to NA and is therefore in
  # NEITHER arm -- it is not "not-never". Expressing that was impossible in the
  # frozen pipeline, whose control arm was the complement of the case value.
  SMOKE_CAT = list(
    source   = c("smoke_statVV", "Smoke_Stat", "smoke_stat"),
    optional = FALSE,
    recode   = function(x) {
      v <- tolower(as_plain(x))
      dplyr::case_when(
        v == "never"                  ~ "NEVER",
        v %in% c("former", "smoker")  ~ "EVER",
        TRUE                          ~ NA_character_
      )
    }
  ),

  # BMI. Verified against the data: Bmi_C is CONTINUOUS despite the `_C`
  # suffix -- 1,846 distinct values spanning 15.63 to 67.86, and the threshold
  # below classified all 10,191 samples with none excluded. Kept optional so an
  # unexpected coding would skip it rather than block the smoking strata.
  # WHO cut-off: obese is BMI >= 30, so 30.0 itself belongs to HI.
  BMI_CAT = list(
    source   = c("bmi_cVV", "bmiVV", "Bmi_C", "BMI", "bmi"),
    optional = TRUE,
    recode   = function(x) {
      v <- suppressWarnings(as.numeric(as_plain(x)))
      dplyr::case_when(
        is.na(v)         ~ NA_character_,
        v <= 10 | v > 80 ~ NA_character_,   # implausible, treat as missing
        v <  30          ~ "LO",
        TRUE             ~ "HI"
      )
    }
  )
)

# Read the strata names, and the minimum stratum size, straight out of
# pipeline/params.yaml so the samplesheet can never disagree with what the
# pipeline will actually run. Uses the yaml package when it is installed and a
# small purpose-built parser otherwise, so this adds no hard dependency.
read_strata_params <- function(params_path) {
  if (!file.exists(params_path)) {
    stop("Cannot find params.yaml at ", params_path,
         " -- set GWAS_PARAMS_FILE in .env.")
  }

  if (requireNamespace("yaml", quietly = TRUE)) {
    y <- yaml::read_yaml(params_path)
    groups <- lapply(y$strata, function(sp) {
      list(column = if (is.null(sp$column)) "" else as.character(sp$column),
           values = if (is.null(sp$values)) character() else as.character(unlist(sp$values)))
    })
    return(list(groups = groups,
                min_n = as.integer(y$min_stratum_n %||% 100L),
                pca_group = as.character(y$pca_group %||% "all")))
  }

  # Fallback: the strata block has a fixed two-level shape that this script and
  # params.yaml both control, so a targeted parser is safe here.
  lines <- readLines(params_path, warn = FALSE)
  start <- grep("^strata:\\s*$", lines)
  if (length(start) == 0) stop("No 'strata:' block found in ", params_path)
  start <- start[1]
  rest <- lines[(start + 1):length(lines)]
  stop_at <- which(grepl("^[^[:space:]#]", rest))
  if (length(stop_at)) rest <- rest[seq_len(stop_at[1] - 1)]

  groups <- list(); current <- NULL
  for (ln in rest) {
    if (grepl("^  [A-Za-z0-9_]+:\\s*$", ln)) {
      current <- trimws(sub(":.*$", "", ln))
      groups[[current]] <- list(column = "", values = character())
    } else if (!is.null(current) && grepl("^    column:", ln)) {
      groups[[current]]$column <- gsub('^\\s*column:\\s*|"|\x27', "", ln)
    } else if (!is.null(current) && grepl("^    values:", ln)) {
      raw <- gsub('^\\s*values:\\s*|\\[|\\]|"|\x27', "", ln)
      vals <- trimws(strsplit(raw, ",")[[1]])
      groups[[current]]$values <- vals[nzchar(vals)]
    }
  }

  scalar <- function(key, default) {
    hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
    if (!length(hit)) return(default)
    trimws(gsub(paste0("^", key, ":\\s*|#.*$|\x22|\x27"), "", hit[1]))
  }

  list(groups = groups,
       min_n = as.integer(scalar("min_stratum_n", "100")),
       pca_group = scalar("pca_group", "all"))
}

`%||%` <- function(a, b) if (is.null(a)) b else a


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

# Per-study PLINK2 filesets. Two layouts exist and both are supported:
#   flat    <genetics>/<study>.psam          (current: data/genetics/)
#   nested  <genetics>/<study>/<study>.psam  (older /scratch staging from 002)
# The flat layout is probed first. Detecting this rather than assuming matters:
# the wrong glob silently finds zero files and every downstream count is 0.
psam_files    <- Sys.glob(file.path(genetics_path, "*.psam"))
genetics_flat <- length(psam_files) > 0
if (!genetics_flat) {
  psam_files <- Sys.glob(file.path(genetics_path, "*", "*.psam"))
}
if (length(psam_files) == 0) {
  stop("No .psam files found under ", genetics_path,
       "\nLooked for both <dir>/*.psam and <dir>/*/*.psam. ",
       "Check GWAS_GENETICS_DIR in .env.")
}
message("Genetics layout: ", if (genetics_flat) "flat" else "nested",
        " -- ", length(psam_files), " study fileset(s) under ", genetics_path)

# Resolve a study's PLINK2 prefix under whichever layout is in use.
study_pfile_prefix <- function(study) {
  if (genetics_flat) file.path(genetics_path, study)
  else file.path(genetics_path, study, study)
}

data_genetic <- psam_files |>
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
  dplyr::mutate(bio = extract_bio(IID_geno, bio_set))

# Per-study proteomics match report. A study whose genotype IIDs match no
# Idepic_Bio is dropped entirely by the filter below, which previously showed up
# only as a missing study much later. Report it here, with example IIDs, so the
# cause is visible: an empty .psam, an unexpected IID format, or a study that
# genuinely has no proteomics.
match_report <- data_genetic |>
  dplyr::group_by(STUDY) |>
  dplyr::summarise(n_read = dplyr::n(),
                   n_matched = sum(!is.na(bio)),
                   .groups = "drop") |>
  dplyr::mutate(pct = ifelse(n_read > 0, 100 * n_matched / n_read, NA_real_)) |>
  dplyr::arrange(pct)

message("\nProteomics match per study (genotype IID -> Idepic_Bio):")
for (i in seq_len(nrow(match_report))) {
  message(sprintf("  %-16s %6d read, %6d matched (%5.1f%%)",
                  match_report$STUDY[i], match_report$n_read[i],
                  match_report$n_matched[i], match_report$pct[i]))
}

unmatched_studies <- match_report$STUDY[match_report$n_matched == 0]
for (st in unmatched_studies) {
  examples <- utils::head(data_genetic$IID_geno[data_genetic$STUDY == st], 5)
  warning("Study '", st, "' matched NO proteomics samples and is excluded. ",
          "Read ", match_report$n_read[match_report$STUDY == st], " genotype IIDs; ",
          if (length(examples)) paste0("examples: ", paste(examples, collapse = ", "))
          else "the .psam is EMPTY",
          ". Compare these against Idepic_Bio in the soma metadata -- a format ",
          "mismatch here silently removes the whole study.",
          call. = FALSE, immediate. = TRUE)
}

data_genetic <- data_genetic |>
  dplyr::filter(!is.na(bio)) |>            # keep genotyped samples with proteomics
  dplyr::left_join(bio_lookup, by = "bio") |>   # adds IID_clean (Idepic) for joins
  dplyr::distinct(IID_geno, STUDY, .keep_all = TRUE)

# strata variables ====
# Join the metadata columns the GWAS needs -- the model covariates (Sex,
# Age_Recr, PlateId) and each strata variable's source column -- then derive the
# strata variables themselves. The completeness check that uses the covariates
# follows further down.

# Resolve each variable's source column against the SomaLogic sample metadata.
# A required variable that cannot be resolved stops the run and prints the
# available column names; an optional one is skipped with a warning.
strata_source <- vapply(names(strata_spec), function(nm) {
  candidates <- strata_spec[[nm]]$source
  hit <- candidates[candidates %in% names(meta_sample)]
  if (length(hit) > 0) return(hit[[1]])

  msg <- paste0("Variable '", nm, "': none of the candidate source columns (",
                paste(candidates, collapse = ", "), ") exist in the SomaLogic ",
                "sample metadata.")
  if (isTRUE(strata_spec[[nm]]$optional)) {
    warning(msg, " Declared optional, so it is SKIPPED. Any stratum in ",
            "params.yaml using it will be empty and will not run.",
            call. = FALSE, immediate. = TRUE)
    return(NA_character_)
  }
  stop(msg, "\nAvailable columns:\n  ",
       paste(sort(names(meta_sample)), collapse = "\n  "),
       "\nAdd the correct name to strata_spec$", nm, "$source near the top of this script.")
}, character(1))

strata_active <- names(strata_source)[!is.na(strata_source)]
message("\nGenotyped samples read from .psam files: ", nrow(data_genetic),
        " across ", dplyr::n_distinct(data_genetic$STUDY), " study file(s) under ", genetics_path)
# 23 filesets but 22 studies in the last run: one .psam contributed no rows.
# Name it rather than leaving the discrepancy to be spotted in the counts.
studies_expected <- sort(unique(tools::file_path_sans_ext(basename(psam_files))))
studies_read     <- sort(unique(data_genetic$STUDY))
studies_missing  <- setdiff(studies_expected, studies_read)
if (length(studies_missing)) {
  warning(length(studies_missing), " fileset(s) contributed no samples and are ",
          "absent from the analysis: ", paste(studies_missing, collapse = ", "),
          ". Check whether the .psam is empty or unreadable.",
          call. = FALSE, immediate. = TRUE)
}

if (nrow(data_genetic) < 500) {
  warning("Only ", nrow(data_genetic), " genotyped samples were found under ", genetics_path,
          ". The full EPIC analysis set is ~9,000. Check GWAS_GENETICS_DIR points at the ",
          "complete per-study PLINK2 filesets and that scratch has not been purged.",
          call. = FALSE, immediate. = TRUE)
}

message("Strata variables active: ",
        paste(sprintf("%s <- %s", strata_active, strata_source[strata_active]), collapse = ", "))
if (length(strata_active) < length(strata_source)) {
  message("Strata variables skipped: ",
          paste(setdiff(names(strata_source), strata_active), collapse = ", "))
}

data_genetic <- data_genetic |>
  dplyr::left_join(
    meta_sample |>
      dplyr::distinct(IID_clean, .keep_all = TRUE) |>
      dplyr::select(IID_clean, Sex, Age_Recr, PlateId,
                    dplyr::all_of(unname(strata_source[strata_active]))),
    by = "IID_clean"
  )

# Derive the analysis variables. These are deliberately NOT part of the
# completeness check below: a sample with unknown smoking status is excluded
# from the ever/never strata but still contributes to the unstratified
# analysis, so dropping it outright would shrink every stratum for no reason.
strata_derived <- character(0)
for (nm in strata_active) {
  raw    <- data_genetic[[strata_source[[nm]]]]
  mapped <- strata_spec[[nm]]$recode(raw)

  message("\n", nm, " (from ", strata_source[[nm]], ")")
  message("  value labels: ", describe_labels(raw))
  message("  raw -> mapped:")
  print(table(raw    = as_plain(raw),
              mapped = ifelse(is.na(mapped), "<excluded>", mapped),
              useNA  = "ifany"))

  n_unmapped <- sum(is.na(mapped))
  pct <- 100 * n_unmapped / length(mapped)
  message(sprintf("  %d/%d (%.1f%%) excluded from every %s stratum.",
                  n_unmapped, length(mapped), pct, nm))
  if (all(is.na(mapped))) {
    msg <- paste0("Variable '", nm, "' mapped NOTHING. The recode in strata_spec ",
                  "does not match how '", strata_source[[nm]], "' is coded -- see ",
                  "the value labels and crosstab above.")
    # A REQUIRED variable mapping nothing is fatal. An OPTIONAL one is skipped:
    # an unexpected coding in a secondary variable must not block the strata
    # that are working, so it is dropped and its strata simply never run.
    if (isTRUE(strata_spec[[nm]]$optional)) {
      warning(msg, " Declared optional, so it is SKIPPED.", call. = FALSE, immediate. = TRUE)
      next
    }
    stop(msg)
  }

  data_genetic[[nm]] <- mapped
  strata_derived <- c(strata_derived, nm)
  if (pct > 20) {
    warning(sprintf("%s: %.1f%% excluded. Check the recode against the data dictionary.",
                    nm, pct), call. = FALSE, immediate. = TRUE)
  }
}

if (length(strata_derived) == 0) {
  stop("No strata variables were derived. Nothing can be stratified.")
}
message("\nStrata variables written to the phenofile: ",
        paste(strata_derived, collapse = ", "))
skipped <- setdiff(names(strata_spec), strata_derived)
if (length(skipped)) {
  message("Skipped (their strata in params.yaml will count 0 and will not run): ",
          paste(skipped, collapse = ", "))
}

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
# Explicit, never derived from the genetics directory. The frozen tree derives
# it as dirname(GWAS_GENETICS_DIR)/phenotype; since both trees read the SAME
# genetics directory, deriving it here would write this analysis's phenofile
# directly over the published one.
pheno_dir <- get_env("GWAS_PHENO_DIR")
if (is.null(pheno_dir) || pheno_dir == "") {
  stop("GWAS_PHENO_DIR is not set. This tree must write its phenotype files to a ",
       "directory of its own -- see pipeline-strata/.env.example.")
}
frozen_pheno_dir <- file.path(dirname(genetics_path), "phenotype")
if (normalizePath(pheno_dir, mustWork = FALSE) ==
    normalizePath(frozen_pheno_dir, mustWork = FALSE)) {
  stop("GWAS_PHENO_DIR (", pheno_dir, ") is the frozen publication analysis's own ",
       "phenotype directory. Writing there would overwrite published inputs. ",
       "Point it somewhere else, e.g. ", frozen_pheno_dir, "-strata")
}
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
# Build the Nextflow samplesheet from the retained studies. Unlike the frozen
# tree, the stratifying column and its value coding are NOT written here: they
# live in params.yaml, and this script only decides WHICH strata each study is
# big enough to run. Reading params.yaml back means the samplesheet can never
# disagree with what the pipeline will actually do.
samplesheet_out <- get_env("GWAS_SAMPLESHEET", "pipeline/samplesheet.csv")
params_file     <- get_env("GWAS_PARAMS_FILE", "pipeline/params.yaml")

sp          <- read_strata_params(params_file)
strata_defs <- sp$groups
min_n       <- sp$min_n
pca_group   <- sp$pca_group

if (!pca_group %in% names(strata_defs)) {
  stop("pca_group '", pca_group, "' in ", params_file,
       " is not a defined stratum (", paste(names(strata_defs), collapse = ", "), ").")
}
if (nzchar(strata_defs[[pca_group]]$column)) {
  stop("pca_group '", pca_group, "' must be an unfiltered stratum (column: \"\") in ",
       params_file, ".")
}

message("\nStrata from ", params_file, ": ",
        paste(names(strata_defs), collapse = ", "),
        "  (min ", min_n, " per stratum; PCA reference: ", pca_group, ")")

# Count a stratum within one study's retained, deduped samples.
count_stratum <- function(study_rows, def) {
  if (!nzchar(def$column)) return(nrow(study_rows))
  if (!def$column %in% names(study_rows)) {
    # The variable was declared optional and skipped. Count 0 so the stratum
    # simply never runs, rather than aborting a run whose other strata are fine.
    return(0L)
  }
  sum(study_rows[[def$column]] %in% def$values, na.rm = TRUE)
}

get_abs_path <- function(p) {
  if (grepl("^/", p)) return(p)
  normalizePath(file.path(getwd(), p), mustWork = FALSE)
}

samplesheet_rows <- list()
for (study in sort(unique(data_pheno$STUDY))) {
  pfile_prefix <- study_pfile_prefix(study)
  for (ext in c(".pgen", ".pvar", ".psam")) {
    if (!file.exists(paste0(pfile_prefix, ext))) {
      stop("Missing PLINK2 file ", paste0(study, ext), " under ", genetics_path)
    }
  }
  sample_file <- paste0(pfile_prefix, ".psam")
  study_rows  <- data_pheno[data_pheno$STUDY == study, , drop = FALSE]

  counts <- vapply(names(strata_defs),
                   function(g) count_stratum(study_rows, strata_defs[[g]]),
                   integer(1))

  # The PCA stratum is mandatory: PCs are computed on it and shared with every
  # other stratum, so a study that omitted it would produce no PCs and silently
  # yield zero results.
  keep <- names(counts)[counts >= min_n]
  if (!pca_group %in% keep) {
    if (counts[[pca_group]] < min_n) {
      message(sprintf("  %s: SKIPPED -- PCA stratum '%s' has only %d samples (min %d).",
                      study, pca_group, counts[[pca_group]], min_n))
      next
    }
    keep <- c(pca_group, keep)
  }
  keep <- unique(c(pca_group, keep))

  cat(sprintf("  %s: %s -> groups=%s\n", study,
              paste(sprintf("%s=%d", names(counts), counts), collapse = " "),
              paste(keep, collapse = ",")))

  samplesheet_rows[[study]] <- data.frame(
    study_id = study,
    pfile = get_abs_path(pfile_prefix),
    sample_file = get_abs_path(sample_file),
    # Pipe-separated on purpose: Nextflow's splitCsv does not honour quoted
    # fields unless quote: is passed, so a comma here is a parsing hazard.
    analysis_groups = paste(keep, collapse = "|"),
    stringsAsFactors = FALSE
  )
}

if (length(samplesheet_rows) == 0) {
  stop("No study met the minimum stratum size (", min_n, "). Nothing to analyse.")
}

final_df <- do.call(rbind, samplesheet_rows)
rownames(final_df) <- NULL
ss_dir <- dirname(samplesheet_out)
if (!dir.exists(ss_dir)) dir.create(ss_dir, recursive = TRUE)
data.table::fwrite(final_df, file = samplesheet_out, sep = ",", quote = FALSE, row.names = FALSE)
message("Wrote samplesheet: ", samplesheet_out, " (", nrow(final_df), " studies).")
