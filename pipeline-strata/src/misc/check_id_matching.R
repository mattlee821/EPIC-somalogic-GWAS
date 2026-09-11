#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# check_id_matching.R  --  DIAGNOSTIC ONLY, writes nothing.
#
# Why: Neuro_01 read 4,415 genotype IIDs and matched 0 proteomics samples, so
# 003 dropped the whole study. Its IIDs look like "00_00000003" -- a SINGLE
# underscore -- while extract_bio() in 003 only recognises Idepic_Bio directly
# or a "<digits>____<digits>" component (FOUR or more underscores). That is the
# suspected cause, and this script tests it against the real data.
#
# It replays the ID matching of 003 for every study in one or more genetics
# directories, and additionally tries every OTHER identifier space in the soma
# metadata, so each study's .psam can be attributed to the ID space it actually
# uses. Run it against the current data/genetics and, if you still have it, the
# older nested /scratch tree, to see what changed.
#
# Usage:
#   Rscript src/misc/check_id_matching.R
#   Rscript src/misc/check_id_matching.R --genetics /path/a,/path/b
#   Rscript src/misc/check_id_matching.R --study Neuro_01     # verbose on one
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr)
})

# --- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
focus_study <- arg_val("--study", NA_character_)

# --- .env (same loader as 003) ---------------------------------------------
load_dot_env <- function(env_file) {
  if (!file.exists(env_file)) return(invisible(FALSE))
  lines <- readLines(env_file, warn = FALSE)
  lines <- lines[!grepl("^\\s*(#|$)", lines)]
  for (line in lines) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) >= 2) {
      key <- trimws(parts[1])
      value <- gsub("^['\"]|['\"]$", "", trimws(paste(parts[-1], collapse = "=")))
      do.call(Sys.setenv, setNames(list(value), key))
    }
  }
  invisible(TRUE)
}
for (p in c(Sys.getenv("GWAS_PROJECT_ROOT"), ".", "..", "../..")) {
  if (nzchar(p) && file.exists(file.path(p, ".env"))) { load_dot_env(file.path(p, ".env")); break }
}
get_env <- function(k, d = NULL) { v <- Sys.getenv(k); if (v == "") d else v }

genetics_dirs <- strsplit(arg_val("--genetics", get_env("GWAS_GENETICS_DIR")), ",")[[1]]
genetics_dirs <- trimws(genetics_dirs[nzchar(trimws(genetics_dirs))])
if (!length(genetics_dirs)) stop("No genetics directory. Pass --genetics or set GWAS_GENETICS_DIR.")

# --- soma metadata: build every candidate ID space --------------------------
rds <- get_env("GWAS_EXTERNAL_SOMA_RDS")
if (is.null(rds) || !file.exists(rds)) stop("GWAS_EXTERNAL_SOMA_RDS not found: ", rds)
cat("Reading soma metadata:", rds, "\n")
meta <- readRDS(rds)$others

collapse_us <- function(x) gsub("_+", "_", x)

id_spaces <- list(
  Idepic_Bio        = unique(na.omit(as.character(meta$Idepic_Bio))),
  Idepic            = unique(na.omit(as.character(meta$Idepic))),
  Idepic_collapsed  = unique(na.omit(collapse_us(as.character(meta$Idepic)))),
  Bio_collapsed     = unique(na.omit(collapse_us(as.character(meta$Idepic_Bio))))
)

cat("\nIdentifier spaces in the soma metadata (n unique, examples):\n")
for (nm in names(id_spaces)) {
  cat(sprintf("  %-18s %6d   e.g. %s\n", nm, length(id_spaces[[nm]]),
              paste(head(id_spaces[[nm]], 3), collapse = ", ")))
}

# --- extract_bio, copied verbatim from 003 ----------------------------------
bio_set <- id_spaces$Idepic_Bio
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

read_psam_iids <- function(f) {
  dt <- data.table::fread(f, showProgress = FALSE)
  data.table::setnames(dt, sub("^#", "", names(dt)))
  as.character(dt$IID)
}

# --- per-directory, per-study report ----------------------------------------
for (gdir in genetics_dirs) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("GENETICS DIR: ", gdir, "\n", sep = "")
  if (!dir.exists(gdir)) { cat("  ** does not exist, skipping **\n"); next }

  psams <- Sys.glob(file.path(gdir, "*.psam"))
  layout <- "flat"
  if (!length(psams)) { psams <- Sys.glob(file.path(gdir, "*", "*.psam")); layout <- "nested" }
  if (!length(psams)) { cat("  ** no .psam files found **\n"); next }
  cat("  layout: ", layout, "   filesets: ", length(psams), "\n\n", sep = "")

  rows <- list()
  for (f in psams) {
    study <- tools::file_path_sans_ext(basename(f))
    if (!is.na(focus_study) && study != focus_study) next
    iids <- read_psam_iids(f)
    iids_collapsed <- collapse_us(iids)

    n_pipeline <- sum(!is.na(extract_bio(iids, bio_set)))
    counts <- vapply(names(id_spaces), function(sp) {
      max(sum(iids %in% id_spaces[[sp]]), sum(iids_collapsed %in% id_spaces[[sp]]))
    }, integer(1))

    best <- names(counts)[which.max(counts)]
    rows[[study]] <- data.frame(
      study = study, n_read = length(iids),
      pipeline = n_pipeline,
      Idepic_Bio = counts[["Idepic_Bio"]], Idepic = counts[["Idepic"]],
      Idepic_collapsed = counts[["Idepic_collapsed"]], Bio_collapsed = counts[["Bio_collapsed"]],
      best_space = best, best_n = max(counts),
      underscores = paste(sort(unique(nchar(gsub("[^_]", "", head(iids, 200))))), collapse = "/"),
      example = paste(head(iids, 2), collapse = " "),
      stringsAsFactors = FALSE)
  }

  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$pipeline / pmax(tab$n_read, 1)), ]

  cat(sprintf("%-16s %7s %9s %11s %8s %10s %9s  %-17s %s\n",
      "STUDY", "READ", "PIPELINE", "Idepic_Bio", "Idepic", "Idepic_col", "Bio_col", "BEST SPACE", "UNDERSCORES"))
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    flag <- if (r$pipeline == 0 && r$best_n > 0) "  <== RECOVERABLE" else
            if (r$pipeline == 0) "  <== no match anywhere" else ""
    cat(sprintf("%-16s %7d %9d %11d %8d %10d %9d  %-17s %s%s\n",
        r$study, r$n_read, r$pipeline, r$Idepic_Bio, r$Idepic,
        r$Idepic_collapsed, r$Bio_collapsed,
        sprintf("%s(%d)", r$best_space, r$best_n), r$underscores, flag))
  }

  cat("\nTotals: read =", sum(tab$n_read),
      "| matched by pipeline =", sum(tab$pipeline),
      "| best-case across spaces =", sum(tab$best_n), "\n")

  broken <- tab[tab$pipeline == 0 & tab$best_n > 0, ]
  if (nrow(broken)) {
    cat("\n", strrep("-", 78), "\n", sep = "")
    cat("STUDIES THE PIPELINE DROPS BUT ANOTHER ID SPACE WOULD RECOVER:\n")
    for (i in seq_len(nrow(broken))) {
      r <- broken[i, ]
      cat(sprintf("\n  %s: %d read, 0 matched by extract_bio, but %d match '%s'\n",
                  r$study, r$n_read, r$best_n, r$best_space))
      cat("    .psam IIDs   : ", r$example, "\n", sep = "")
      cat("    that space   : ", paste(head(id_spaces[[r$best_space]], 2), collapse = " "), "\n", sep = "")
      cat("    underscore counts in IID: ", r$underscores, "\n", sep = "")
    }
    cat("\n  extract_bio() accepts an Idepic_Bio directly, or a '<digits>____<digits>'\n")
    cat("  component with FOUR OR MORE underscores. An IID with a single underscore\n")
    cat("  matches neither, so the study is silently dropped.\n")
  }
}
cat("\nDone. Nothing was written.\n")
