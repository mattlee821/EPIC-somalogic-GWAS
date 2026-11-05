#!/usr/bin/env Rscript
# Purpose: Generate the sample de-duplication comparison figure.
#
# Usage:
#   Rscript src/misc/plot_dedup.R
#   Rscript src/misc/plot_dedup.R \
#     --out analysis/images/dedup_fig.png --width 14 --height 5
#
# Requires: data.table, ggplot2 (both in the r_gwas conda env).
# Run with:
#   conda activate /WHERE/Somalogic/work/GWAS/env/r_gwas
#
# Reads GWAS_EXTERNAL_SOMA_RDS, GWAS_GENETICS_DIR, GWAS_MIN_STUDY_N from .env.
# Reports final N per strategy to stdout.

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# -- args ---------------------------------------------------------------------
args          <- commandArgs(trailingOnly = TRUE)
out_path      <- "analysis/images/dedup_fig.png"
fig_w         <- 14
fig_h         <- 5
smiss_thresh  <- 0.05   # exclude samples with >5% missing dosage
plink2_bin    <- Sys.which("plink2")
i <- 1
while (i <= length(args)) {
  if (args[i] == "--out") {
    out_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--width") {
    fig_w <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--height") {
    fig_h <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--smiss-threshold") {
    smiss_thresh <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--plink2") {
    plink2_bin <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

# -- .env ---------------------------------------------------------------------
load_dot_env <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*(#|$)", lines)]
  for (line in lines) {
    kv <- strsplit(line, "=")[[1]]
    if (length(kv) >= 2) {
      key <- trimws(kv[1])
      val <- trimws(paste(kv[-1], collapse = "="))
      val <- gsub("^['\"]|['\"]$", "", val)
      do.call(Sys.setenv, setNames(list(val), key))
    }
  }
}

proj_root <- Sys.getenv("GWAS_PROJECT_ROOT")
for (p in c(if (proj_root != "") proj_root, ".", "..", "../..")) {
  f <- file.path(p, ".env")
  if (file.exists(f)) {
    load_dot_env(f)
    if (file.exists(file.path(p, "server.env"))) load_dot_env(file.path(p, "server.env"))
    break
  }
}

get_env <- function(k, d = NULL) {
  v <- Sys.getenv(k)
  if (v == "") d else v
}

# -- SomaLogic bio -> IID mapping ---------------------------------------------
rds_path <- get_env("GWAS_EXTERNAL_SOMA_RDS")
if (is.null(rds_path) || !file.exists(rds_path))
  stop("GWAS_EXTERNAL_SOMA_RDS not found. Check .env.")

message("Loading SomaLogic RDS: ", rds_path)
soma <- readRDS(rds_path)

others_dt <- as.data.table(soma$others)
others_dt[, IID_clean := gsub("_+", "_", Idepic)]
others_dt[, bio       := Idepic_Bio]

bio_set    <- unique(na.omit(others_dt$bio))
bio_lookup <- unique(others_dt[!is.na(bio), .(bio, IID_clean)])

# Compound IID_geno values (e.g. "00____00000002") embed the bio ID;
# extract_bio recovers it via regex.
extract_bio <- function(iid) {
  hit  <- ifelse(iid %in% bio_set, iid, NA_character_)
  need <- is.na(hit)
  if (any(need)) {
    cand      <- regmatches(
      iid[need],
      gregexpr("[0-9]+_{4,}[0-9]+", iid[need])
    )
    hit[need] <- vapply(cand, function(ms) {
      m <- ms[ms %in% bio_set]
      if (length(m)) m[[1]] else NA_character_
    }, character(1))
  }
  hit
}

# -- load PSAMs ---------------------------------------------------------------
genetics_path <- get_env("GWAS_GENETICS_DIR")
if (is.null(genetics_path) || !dir.exists(genetics_path))
  stop("GWAS_GENETICS_DIR not found. Run 002_data-genetics.sh first.")

psam_files <- Sys.glob(file.path(genetics_path, "*", "*.psam"))
if (length(psam_files) == 0)
  stop("No .psam files found under: ", genetics_path)

message("Reading ", length(psam_files), " PSAM files...")

read_psam <- function(f) {
  study <- tools::file_path_sans_ext(basename(f))
  dt    <- fread(f)
  setnames(dt, sub("^#", "", names(dt)))
  pheno_col <- intersect(c("PHENO", "PHENO1"), names(dt))
  pheno_col <- if (length(pheno_col)) pheno_col[1] else NA_character_
  data.table(
    IID_geno  = as.character(dt$IID),
    STUDY     = study,
    PHENO     = if (!is.na(pheno_col)) {
      as.character(dt[[pheno_col]])
    } else {
      NA_character_
    },
    has_pheno = !is.na(pheno_col)
  )
}

data_raw <- rbindlist(lapply(psam_files, read_psam))
data_raw[, bio := extract_bio(IID_geno)]
data_raw <- data_raw[!is.na(bio)]
data_raw <- merge(data_raw, bio_lookup, by = "bio", all.x = TRUE)
data_raw <- unique(data_raw, by = c("IID_geno", "STUDY"))

message(
  "Matched genotype-proteomics records: ", nrow(data_raw),
  " across ", uniqueN(data_raw$STUDY), " studies"
)

# -- comprehensive sample completeness assessment (pre-dedup) -----------------
# Mirror 003_data-phenotype.R exactly: remove any sample missing data for ANY
# GWAS covariate (Sex, Age_Recr, PHENO) before dedup, so the figure's matched
# set -- and therefore its N -- matches the GWAS. REGENIE drops samples with any
# missing covariate from every test, so an incomplete sample is never analysed.
# Proteins are outcomes (per-trait missingness handled by REGENIE) and PlateId
# is not a covariate, so neither triggers removal. PHENO (disease status, from
# the .psam) only counts where the study .psam carries a PHENO column.
cov_lookup <- unique(others_dt, by = "IID_clean")[, .(IID_clean, Sex, Age_Recr)]
data_raw <- merge(data_raw, cov_lookup, by = "IID_clean", all.x = TRUE)

# Robust to haven_labelled (SPSS) columns from the soma RDS, whose as.character
# method errors: strip the class to the underlying atomic vector first.
is_missing <- function(x) {
  x <- unclass(x)
  attributes(x) <- NULL
  is.na(x) | trimws(as.character(x)) %in% c("", "NA")
}

incomplete <- is_missing(data_raw$Sex) |
  is_missing(data_raw$Age_Recr) |
  (data_raw$has_pheno & is_missing(data_raw$PHENO))

message("Sample completeness over GWAS covariates (Sex, Age_Recr, PHENO):")
message("  missing Sex: ",      sum(is_missing(data_raw$Sex)))
message("  missing Age_Recr: ", sum(is_missing(data_raw$Age_Recr)))
message("  missing PHENO: ",    sum(data_raw$has_pheno & is_missing(data_raw$PHENO)))

n_excl <- sum(incomplete)
if (n_excl > 0) {
  message("Excluding ", n_excl,
          " samples with any missing GWAS covariate, before dedup:")
  excl_tab <- data_raw[incomplete, .N, by = STUDY][order(-N)]
  for (i in seq_len(nrow(excl_tab))) {
    message("  ", excl_tab$STUDY[i], ": ", excl_tab$N[i])
  }
  excl_iids <- data_raw$IID_geno[incomplete]
  excl_out  <- sub("\\.png$", "_covariate_excluded.txt", out_path)
  writeLines(excl_iids, excl_out)
  message("  IIDs written to: ", excl_out)
  data_raw <- data_raw[!incomplete]
  message("  Records after exclusion: ", nrow(data_raw))
} else {
  message("No samples with missing GWAS covariates to exclude.")
}

has_pheno <- any(!is.na(data_raw$PHENO))
if (has_pheno) {
  message("PHENO column detected -- figure will show case/control split")
} else {
  message("No PHENO column -- figure will show total counts only")
}

# -- dedup logic --------------------------------------------------------------
min_n <- as.integer(get_env("GWAS_MIN_STUDY_N", "100"))

# rank_desc = FALSE : each person assigned to their smallest retained study
# rank_desc = TRUE  : each person assigned to their largest retained study
run_dedup <- function(dt, rank_desc) {
  prov_size <- dt[, .(study_n = .N), by = STUDY]
  
  # Provisional assignment to determine which studies survive the min_n cut
  tmp <- merge(dt, prov_size, by = "STUDY")
  ord <- if (rank_desc) c(-1L, 1L) else c(1L, 1L)
  setorderv(tmp, c("study_n", "STUDY"), order = ord)
  prov_assigned  <- tmp[, .SD[1], by = bio]
  prov_retained  <- prov_assigned[, .N, by = STUDY][N >= min_n, STUDY]
  
  # Re-rank retained studies for final assignment
  retained_size <- prov_size[STUDY %in% prov_retained]
  setorderv(retained_size, c("study_n", "STUDY"), order = ord)
  retained_size[, rk := .I]
  
  final <- merge(
    dt[STUDY %in% prov_retained],
    retained_size[, .(STUDY, rk)],
    by = "STUDY"
  )
  setorder(final, rk)
  final <- final[, .SD[1], by = bio]
  final[, rk := NULL]
  final
}

# run_dedup_all: same ranking but no min_n filter -- used for the figure only
# so that all studies (including those below the retention threshold) are shown.
run_dedup_all <- function(dt, rank_desc) {
  prov_size <- dt[, .(study_n = .N), by = STUDY]
  ord <- if (rank_desc) c(-1L, 1L) else c(1L, 1L)
  tmp <- merge(dt, prov_size, by = "STUDY")
  setorderv(tmp, c("study_n", "STUDY"), order = ord)
  assigned <- tmp[, .SD[1], by = bio]
  assigned[, study_n := NULL]
  assigned
}

d_small <- run_dedup(data_raw, rank_desc = FALSE)
d_small[, strategy := "keep_in_smallest"]

d_large <- run_dedup(data_raw, rank_desc = TRUE)
d_large[, strategy := "keep_in_largest"]

d_small_all <- run_dedup_all(data_raw, rank_desc = FALSE)
d_large_all <- run_dedup_all(data_raw, rank_desc = TRUE)

# Figure datasets: retained studies use actual dedup N (matches GWAS exactly);
# non-retained studies show samples not already claimed by any retained study.
d_fig_small <- rbindlist(list(d_small, d_small_all[!bio %in% d_small$bio]), fill = TRUE)
d_fig_small[, strategy := "keep_in_smallest"]

d_fig_large <- rbindlist(list(d_large, d_large_all[!bio %in% d_large$bio]), fill = TRUE)
d_fig_large[, strategy := "keep_in_largest"]

message(
  "\n-- Strategy: reassign duplicates to smaller studies (pipeline) ----------"
)
message("   Studies retained : ", uniqueN(d_small$STUDY))
message("   Participants (N) : ", nrow(d_small))
message("   Per study:")
per_study <- d_small[, .N, by = STUDY][order(-N)]
invisible(apply(per_study, 1, function(r) {
  message("     ", r["STUDY"], ": ", r["N"])
}))

message("\n-- Strategy: remove duplicates from smaller studies ---------------")
message("   Studies retained : ", uniqueN(d_large$STUDY))
message("   Participants (N) : ", nrow(d_large))

# -- count by study and group -------------------------------------------------
# Samples classified as Case (PHENO == "2") or Control (all others, including
# PHENO = NA / "0" / "-9"). No samples are dropped — figure N must equal the
# Participants (N) reported above.
CASES_VAL <- "2"

label_group <- function(dt) {
  if (has_pheno) {
    dt[, group := ifelse(
      !is.na(PHENO) & PHENO == CASES_VAL, "Case", "Control"
    )]
  } else {
    dt[, group := "Total"]
  }
  dt
}

tab <- rbindlist(list(
  label_group(copy(d_fig_small)),
  label_group(copy(d_fig_large))
))
tab <- tab[, .(n = .N), by = .(strategy, STUDY, group)]
tab[, group := factor(group, levels = c("Case", "Control", "Total"))]

# -- figure -------------------------------------------------------------------
fill_cols <- if (has_pheno) {
  c("Case" = "#2c5f8a", "Control" = "#9dc3e6")
} else {
  c("Total" = "#2c5f8a")
}

strat_labels <- c(
  keep_in_smallest = "Duplicates reassigned to smaller studies",
  keep_in_largest  = "Duplicates removed from smaller studies"
)

build_panel <- function(strat_code) {
  tot <- tab[strategy == strat_code, .(N = sum(n)), by = STUDY]
  setorder(tot, -N)
  lvls      <- tot$STUDY
  n_large   <- sum(tot$N >= min_n)
  sum_large <- sum(tot$N[tot$N >= min_n])
  sum_all   <- sum(tot$N)
  ymax      <- max(tot$N)
  
  df <- copy(tab[strategy == strat_code])
  df[,  STUDY := factor(STUDY, levels = lvls)]
  tot[, STUDY := factor(STUDY, levels = lvls)]
  
  ggplot(df, aes(STUDY, n, fill = group)) +
    geom_col(width = 0.8) +
    geom_text(
      data     = df[n >= 15],
      aes(label = n),
      position = position_stack(vjust = 0.5),
      size = 2.4, colour = "white", fontface = "bold"
    ) +
    geom_text(
      data        = tot,
      aes(STUDY, N, label = N),
      inherit.aes = FALSE,
      vjust = -0.4, size = 2.8, fontface = "bold"
    ) +
    geom_vline(
      xintercept = n_large + 0.5,
      linetype = "dashed", colour = "firebrick"
    ) +
    annotate(
      "text",
      x = n_large + 0.7, y = ymax * 0.94,
      label = paste0("< ", min_n), hjust = 0,
      colour = "firebrick", size = 3.2
    ) +
    scale_fill_manual(values = fill_cols, name = NULL) +
    scale_y_continuous(
      limits = c(0, ymax * 1.12),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title    = strat_labels[[strat_code]],
      subtitle = paste0(
        "Total: N = ", format(sum_all, big.mark = ","),
        " | Retained (>= ", min_n, "): N = ",
        format(sum_large, big.mark = ",")
      ),
      x = NULL, y = "Matched samples"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position    = "bottom"
    )
}

p1 <- build_panel("keep_in_smallest")
p2 <- build_panel("keep_in_largest")

# -- save (base grid combines the two panels) ---------------------------------
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

png(out_path, width = fig_w, height = fig_h, units = "in", res = 150)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

message("\nFigure saved: ", out_path)
