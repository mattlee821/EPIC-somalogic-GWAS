#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# find_proteins.R -- which SomaScan aptamers exist for a list of genes?
#
# Searches data$aux (the feature annotation) by gene symbol, then falls back to
# a substring search of Target for anything not found, and finally confirms each
# SeqId is actually a column in the phenofile -- which is what --proteins is
# matched against. An aptamer present in aux but absent from the phenofile
# cannot be analysed, so both checks matter.
#
# Usage:
#   Rscript src/misc/find_proteins.R
#   Rscript src/misc/find_proteins.R --genes "CHRNA5,BDNF,FTO"
#   Rscript src/misc/find_proteins.R --out data/phenotype-strata/smoking-panel.txt
# ---------------------------------------------------------------------------

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}

# Default panel: nicotinic receptor cluster, nicotine metabolism, neuro /
# addiction, and the canonical BMI loci.
default_genes <- c("CHRNA5", "CHRNA3", "CHRNB4", "CYP2A6", "EGLN2", "BDNF",
                   "DBH", "SLC6A4", "DRD2", "ANKK1", "FTO", "MC4R", "TMEM18")

genes <- arg_val("--genes", paste(default_genes, collapse = ","))
# "DRD2/ANKK1" style input is split into separate symbols.
genes <- unique(toupper(trimws(unlist(strsplit(genes, "[,/]")))))
genes <- genes[nzchar(genes)]

out_file <- arg_val("--out", NA_character_)

# --- .env ------------------------------------------------------------------
load_dot_env <- function(f) {
  if (!file.exists(f)) return(invisible(FALSE))
  for (line in grep("^\\s*(#|$)", readLines(f, warn = FALSE), invert = TRUE, value = TRUE)) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) >= 2) do.call(Sys.setenv, setNames(
      list(gsub("^['\"]|['\"]$", "", trimws(paste(parts[-1], collapse = "=")))), trimws(parts[1])))
  }
}
for (p in c(Sys.getenv("GWAS_PROJECT_ROOT"), ".", "..", "../..")) {
  if (nzchar(p) && file.exists(file.path(p, ".env"))) { load_dot_env(file.path(p, ".env")); break }
}
get_env <- function(k, d = NULL) { v <- Sys.getenv(k); if (v == "") d else v }

rds <- get_env("GWAS_EXTERNAL_SOMA_RDS")
if (is.null(rds) || !file.exists(rds)) stop("GWAS_EXTERNAL_SOMA_RDS not found: ", rds)
aux <- readRDS(rds)$aux
cat("Feature annotation:", nrow(aux), "aptamers\n")

# Gene symbol fields can hold several genes; split into tokens for exact matching.
tokens <- strsplit(toupper(as.character(aux$EntrezGeneSymbol)), "[ ,;|/]+")

hits <- lapply(genes, function(g) {
  exact <- which(vapply(tokens, function(tk) g %in% tk, logical(1)))
  if (length(exact)) return(list(idx = exact, how = "symbol"))
  # Fall back to the Target description: catches symbol aliases and renamings.
  loose <- grep(g, toupper(as.character(aux$Target)), fixed = TRUE)
  if (length(loose)) return(list(idx = loose, how = "Target substring"))
  list(idx = integer(0), how = NA_character_)
})
names(hits) <- genes

# --- phenofile columns -----------------------------------------------------
pheno <- get_env("GWAS_PHENO_FILE")
pheno_cols <- character(0)
if (!is.null(pheno) && file.exists(pheno)) {
  hdr <- readLines(pheno, n = 1L)
  pheno_cols <- strsplit(hdr, "\t")[[1]]
  cat("Phenofile:", pheno, "--", length(pheno_cols), "columns\n")
} else {
  cat("NOTE: GWAS_PHENO_FILE not found; skipping the phenofile check.\n")
}

cat("\n", strrep("=", 76), "\n", sep = "")
found <- character(0)
for (g in genes) {
  h <- hits[[g]]
  if (!length(h$idx)) { cat(sprintf("%-10s  NOT ON THE PANEL\n", g)); next }
  cat(sprintf("%-10s  %d aptamer(s)  [matched by %s]\n", g, length(h$idx), h$how))
  for (i in h$idx) {
    sid <- as.character(aux$SeqId[i])
    in_ph <- if (!length(pheno_cols)) "?" else if (sid %in% pheno_cols) "yes" else "NO"
    cat(sprintf("    %-12s %-46s symbol=%-12s in phenofile: %s\n",
                sid, substr(as.character(aux$Target[i]), 1, 46),
                as.character(aux$EntrezGeneSymbol[i]), in_ph))
    if (in_ph != "NO") found <- c(found, sid)
  }
}
found <- unique(found)

cat("\n", strrep("=", 76), "\n", sep = "")
missing <- genes[vapply(hits, function(h) length(h$idx) == 0, logical(1))]
cat(sprintf("Requested %d genes: %d on the panel, %d absent.\n",
            length(genes), length(genes) - length(missing), length(missing)))
if (length(missing)) cat("Absent: ", paste(missing, collapse = ", "), "\n", sep = "")
cat("Usable aptamers (present in the phenofile): ", length(found), "\n", sep = "")

if (!is.na(out_file) && length(found)) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(found, out_file)
  cat("\nWrote ", out_file, " with ", length(found), " SeqIds.\n", sep = "")
  cat("Run it with:\n  sbatch src/005_run.sh --out <dir> --proteins ", normalizePath(out_file, mustWork = FALSE), "\n", sep = "")
} else if (is.na(out_file)) {
  cat("\nPass --out <file> to write these SeqIds as a --proteins list.\n")
}
