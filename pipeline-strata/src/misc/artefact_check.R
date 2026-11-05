#!/usr/bin/env Rscript
# METHODS.md 9.1, gating items 1-2: is the oriented mean(z) explained by the
# phenotype transform rather than by biology?
#
# The artefact needs BOTH a between-stratum mean shift AND raw skewness, and it
# produces a DIRECTIONAL shift in z once oriented by sign(beta_meta). So the
# test is whether per-protein oriented mean(z) tracks shift x skew.
#
# Usage: Rscript artefact_check.R <pqtl_comparison.tsv> <protein_distribution.tsv> <annot.tsv>
a <- commandArgs(TRUE)
cmp <- read.delim(a[1], stringsAsFactors = FALSE)
dst <- read.delim(a[2], stringsAsFactors = FALSE)
ann <- read.delim(a[3], stringsAsFactors = FALSE)

# ---- per-protein oriented mean(z), primary QC-passing loci only ------------
bref <- grep("^BETA_meta$", names(cmp), value = TRUE)
if (!length(bref)) stop("no BETA_meta column in ", a[1])
d <- cmp[cmp$QC_PASS %in% c(TRUE, "True", "TRUE") & cmp$TIER == "primary", ]
d$zo <- d$Z_DIFF * sign(d[[bref]])
per <- aggregate(zo ~ PROTEIN, data = d, FUN = function(v) c(m = mean(v), n = length(v)))
per <- data.frame(protein = per$PROTEIN, zo = per$zo[, "m"], n = per$zo[, "n"])

# ---- per-protein mean shift and skewness ----------------------------------
s <- dst[dst$stratum_var == "SMOKE_CAT", ]
E <- s[s$stratum == "EVER", ]; N <- s[s$stratum == "NEVER", ]
stopifnot(nrow(E) > 0, nrow(N) > 0)
N <- N[match(E$protein, N$protein), ]

# SIGNS MATTER. The artefact's direction is set by which arm sits in the
# compressed region of the inverse-rank scale. For a right-skewed protein the
# upper tail is compressed, so if EVER has the higher mean, EVER's slope is the
# smaller one and the oriented z goes NEGATIVE. Hence:
#
#     predicted oriented z  proportional to  -sign(mean_ever - mean_never) * skew
#
# A first version of this script used abs(shift) and max(|skew|), which destroys
# exactly the information the test depends on.
pd <- data.frame(protein = E$protein,
                 shift_s = (E$mean - N$mean) / ((E$sd + N$sd) / 2),
                 skew_s  = (E$skewness + N$skewness) / 2)
pd$shift <- abs(pd$shift_s)
pd$skew  <- abs(pd$skew_s)
pd$pred  <- -pd$shift_s * pd$skew_s          # the artefact's signed prediction
pd$risk  <- pd$shift * pd$skew

m <- merge(per, pd, by = "protein")
m <- merge(m, ann[, c("seqid", "source")], by.x = "protein", by.y = "seqid", all.x = TRUE)
m$targeted <- m$source == "targeted"

cat(sprintf("%d proteins with primary loci (%d targeted)\n\n",
            nrow(m), sum(m$targeted, na.rm = TRUE)))

# ---- the 11 targeted proteins, named --------------------------------------
tg <- m[which(m$targeted), c("protein","n","zo","shift","skew","risk","shift_s","skew_s","pred")]
tg <- tg[order(-abs(tg$zo)), ]
cat("TARGETED proteins, ordered by |oriented mean z|:\n")
print(data.frame(protein = tg$protein, loci = tg$n,
                 oriented_z = round(tg$zo, 3),
                 shift_signed = round(tg$shift_s, 3),
                 skew_signed = round(tg$skew_s, 3),
                 artefact_pred = round(tg$pred, 3)),
      row.names = FALSE)
cat("\n  artefact_pred is -shift_signed x skew_signed: if the transform is\n")
cat("  driving things, oriented_z should follow its SIGN and magnitude.\n")

cat(sprintf("\n  targeted   : median shift %.3f, median risk %.3f\n",
            median(tg$shift, na.rm = TRUE), median(tg$risk, na.rm = TRUE)))
rn <- m[which(!m$targeted), ]
cat(sprintf("  random     : median shift %.3f, median risk %.3f\n",
            median(rn$shift, na.rm = TRUE), median(rn$risk, na.rm = TRUE)))
cat(sprintf("  Wilcoxon on shift, targeted vs random: P = %.3g\n",
            tryCatch(wilcox.test(tg$shift, rn$shift)$p.value, error = function(e) NA)))

# ---- does oriented mean(z) track the artefact predictor? -------------------
cat("\nDIAGNOSTIC REGRESSION  (random proteome only -- the targeted set is the\n")
cat("hypothesis, so it cannot also be the test)\n\n")
w <- rn[is.finite(rn$zo) & is.finite(rn$pred), ]
fit <- lm(zo ~ pred, data = w, weights = n)
print(summary(fit)$coefficients)
cat(sprintf("\n  R^2 = %.4f on %d proteins\n", summary(fit)$r.squared, nrow(w)))
cat(sprintf("  Spearman(oriented z, artefact_pred) = %.3f, P = %.3g\n",
            cor(w$zo, w$pred, method = "spearman"),
            suppressWarnings(cor.test(w$zo, w$pred, method = "spearman")$p.value)))

# Use the random-set slope to predict what the artefact alone would give the
# targeted proteins, and compare with what was observed.
b <- coef(fit)
tg$expected <- b[1] + b[2] * tg$pred
cat("\nTARGETED: observed oriented z vs what the random-set artefact model predicts\n")
print(data.frame(protein = tg$protein, observed = round(tg$zo, 3),
                 predicted = round(tg$expected, 3),
                 residual = round(tg$zo - tg$expected, 3)), row.names = FALSE)
cat("\nA residual near zero means the transform explains it. A residual that\n")
cat("stays large is the part the artefact cannot account for.\n")
