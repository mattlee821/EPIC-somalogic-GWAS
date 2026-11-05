# Quantify the transform-artefact risk per protein from protein_distribution.tsv.
# The mechanism needs BOTH a stratum mean shift and raw skewness, so the risk
# statistic is their product. See METHODS.md 9.1.
f <- commandArgs(TRUE)[1]
v <- if (length(commandArgs(TRUE)) > 1) commandArgs(TRUE)[2] else "SMOKE_CAT"
d <- read.delim(f, stringsAsFactors = FALSE)
d <- d[d$stratum_var == v, ]
lv <- unique(d$stratum)
stopifnot(length(lv) == 2)
a <- d[d$stratum == lv[1], ]; b <- d[d$stratum == lv[2], ]
b <- b[match(a$protein, b$protein), ]

shift <- abs((a$mean - b$mean) / ((a$sd + b$sd) / 2))
skew  <- pmax(abs(a$skewness), abs(b$skewness))
risk  <- shift * skew
ok    <- is.finite(risk)

q <- function(x, l) {
  p <- quantile(x, c(.5, .9, .95, .99, 1), na.rm = TRUE)
  cat(sprintf("  %-26s median %.3f  p90 %.3f  p95 %.3f  p99 %.3f  max %.3f\n",
              l, p[1], p[2], p[3], p[4], p[5]))
}
cat(sprintf("%s: %s vs %s, %d proteins\n\n", v, lv[1], lv[2], sum(ok)))
q(shift[ok], "|mean shift| (SD)");  q(skew[ok], "|skew| (max over arms)")
q(risk[ok],  "risk = shift x skew")

# MEAN SHIFT is the clean observable and the necessary condition: a protein with
# no stratum mean shift is immune to the transform artefact whatever its shape.
# The skewness column is measured on the ALREADY-TRANSFORMED phenotype and is
# contaminated by the residualisation and by study subsetting -- smoking and BMI
# return near-identical skewness (0.548 vs 0.561) despite very different shifts --
# so the product is indicative, not calibrated. Filter on shift.
cat("\nProteins excluded by |mean shift| (the recommended filter):\n")
for (thr in c(0.05, 0.10, 0.15, 0.20, 0.30)) {
  n <- sum(shift > thr, na.rm = TRUE)
  cat(sprintf("  |shift| > %-5.2f : %5d excluded (%5.1f%%), %5d retained\n",
              thr, n, 100*n/sum(ok), sum(ok) - n))
}

cat("\nProteins excluded by the risk product (indicative only):\n")
for (thr in c(0.05, 0.1, 0.15, 0.2, 0.3)) {
  n <- sum(risk > thr, na.rm = TRUE)
  cat(sprintf("  risk > %-5.2f : %5d excluded (%5.1f%%), %5d retained\n",
              thr, n, 100*n/sum(ok), sum(ok) - n))
}

cat("\nWorst-case var(z) inflation at each shift quantile, assuming the protein\n")
cat("is maximally skewed on the raw scale (lognormal sigma=1), at a z_main of 40:\n")
xs <- c(0.027, 0.10, 0.30, 0.63); ys <- c(-2.67, -14.97, -31.70, -39.31)
for (lbl in c(0.5, 0.9, 0.95, 0.99, 1)) {
  s <- unname(quantile(shift, lbl, na.rm = TRUE))
  r <- approx(xs, ys, xout = pmin(pmax(s, min(xs)), max(xs)))$y
  zd <- 40 * abs(r) / 100 / sqrt(2)
  cat(sprintf("  q%-5.2f shift %.3f -> slope diff %6.1f%%  z_diff %5.2f  var(z) +%6.1f\n",
              lbl, s, r, zd, zd^2))
}
cat("  (a Gaussian protein is immune at ANY shift; this is the ceiling)\n")
cat("\nWorst 15 proteins by risk:\n")
o <- order(-risk)[1:15]
print(data.frame(protein = a$protein[o], shift = round(shift[o], 3),
                 skew = round(skew[o], 3), risk = round(risk[o], 3)),
      row.names = FALSE)
