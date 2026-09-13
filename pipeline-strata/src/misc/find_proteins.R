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
# Random aptamers appended as null controls: they carry no hypothesis, so they
# show how the ever/never comparison behaves when nothing is expected. Seeded so
# the same panel is reproducible across runs.
n_random <- as.integer(arg_val("--random", "0"))
seed     <- as.integer(arg_val("--seed", "42"))

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

# --- random null controls ---------------------------------------------------
random_ids <- character(0)
if (n_random > 0) {
  if (!length(pheno_cols)) {
    stop("--random needs the phenofile to draw from; GWAS_PHENO_FILE was not found.")
  }
  # Draw only from real aptamers that are present in the phenofile, so the
  # controls are genuine measured proteins rather than metadata columns.
  pool <- setdiff(intersect(as.character(aux$SeqId), pheno_cols), found)
  if (length(pool) < n_random) {
    stop("Only ", length(pool), " aptamers available to sample ", n_random, " from.")
  }
  set.seed(seed)
  random_ids <- sort(sample(pool, n_random))

  cat("\n", strrep("=", 76), "\n", sep = "")
  cat("RANDOM NULL CONTROLS (seed ", seed, "):\n", sep = "")
  for (sid in random_ids) {
    i <- match(sid, as.character(aux$SeqId))
    cat(sprintf("    %-12s %-46s symbol=%s\n", sid,
                substr(as.character(aux$Target[i]), 1, 46),
                as.character(aux$EntrezGeneSymbol[i])))
  }
}

panel <- c(found, random_ids)

if (!is.na(out_file) && length(panel)) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  # The pipeline reads this as a bare SeqId list -- every non-empty line becomes
  # a protein name, so it can carry no comments. The provenance of each entry
  # goes in a companion file the pipeline never reads.
  writeLines(panel, out_file)

  annot <- sub("\\.txt$", "", out_file)
  annot <- paste0(annot, ".annot.tsv")
  idx <- match(panel, as.character(aux$SeqId))
  write.table(
    data.frame(
      seqid  = panel,
      source = c(rep("targeted", length(found)), rep("random_control", length(random_ids))),
      gene   = as.character(aux$EntrezGeneSymbol)[idx],
      target = as.character(aux$Target)[idx],
      stringsAsFactors = FALSE
    ),
    annot, sep = "\t", quote = FALSE, row.names = FALSE
  )

  cat("\nWrote ", out_file, "\n", sep = "")
  cat("  ", length(found), " targeted + ", length(random_ids),
      " random = ", length(panel), " SeqIds\n", sep = "")
  cat("Wrote ", annot, " (provenance; not read by the pipeline)\n", sep = "")
  cat("\nRun it with:\n  sbatch src/005_run.sh --out <dir> --proteins ",
      normalizePath(out_file, mustWork = FALSE), "\n", sep = "")
} else if (is.na(out_file)) {
  cat("\nPass --out <file> to write these SeqIds as a --proteins list.\n")
}
