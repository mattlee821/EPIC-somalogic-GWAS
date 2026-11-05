#!/usr/bin/env Rscript
# List every index variant behind the TARGETED result, with enough context to
# judge METHODS.md 7.2: are these loci at the smoking-behaviour genes themselves
# (where SNP -> Smoking is established, so the collider path is open), or
# elsewhere in the genome?
#
# Usage: Rscript targeted_loci.R <pqtl_comparison.tsv> <annot.tsv>
a <- commandArgs(TRUE)
cmp <- read.delim(a[1], stringsAsFactors = FALSE)
ann <- read.delim(a[2], stringsAsFactors = FALSE)
tg  <- ann$seqid[ann$source == "targeted"]

cat(sprintf("Targeted aptamers on the panel: %d\n", length(tg)))
have <- intersect(tg, unique(cmp$PROTEIN))
cat(sprintf("  with any index variant     : %d\n", length(have)))

d <- cmp[cmp$PROTEIN %in% tg, ]
prim <- d[d$TIER == "primary" & d$QC_PASS %in% c(TRUE, "True", "TRUE"), ]
cat(sprintf("  with a PRIMARY QC-passing locus: %d\n",
            length(unique(prim$PROTEIN))))
cat(sprintf("  aptamers with NOTHING          : %s\n\n",
            paste(setdiff(tg, unique(prim$PROTEIN)), collapse = ", ")))

# gene each aptamer targets, from the panel annotation
g <- setNames(ann$gene, ann$seqid)
prim$TARGET_GENE <- g[prim$PROTEIN]
prim$ORIENTED_Z  <- prim$Z_DIFF * sign(prim$BETA_meta)

keep <- c("PROTEIN","TARGET_GENE","GENE","SNP","CHR","POS","CIS","EA","OA",
          "BETA_ever","SE_ever","BETA_never","SE_never","DELTA","Z_DIFF",
          "ORIENTED_Z","P_DIFF","P_DIFF_BONF","EAF_meta","MAX_HET_ISQ")
keep <- keep[keep %in% names(prim)]
out <- prim[order(prim$PROTEIN, -abs(prim$ORIENTED_Z)), keep]
num <- sapply(out, is.numeric)
out[num] <- lapply(out[num], function(v) signif(v, 3))

cat("PRIMARY index variants behind the targeted result:\n")
print(out, row.names = FALSE)

cat("\ncis / trans split: ")
print(table(prim$CIS))
cat("\nChromosomes carrying these variants: ",
    paste(sort(unique(prim$CHR)), collapse = ", "), "\n", sep = "")
cat("\nKnown smoking-behaviour regions to check against (GRCh38):\n")
cat("  CHRNA5-A3-B4  chr15 : 78.5-78.9 Mb\n")
cat("  CYP2A6        chr19 : 40.8-41.0 Mb\n")
cat("  EGLN2/CYP2A6  chr19 : 40.7-41.1 Mb\n")
cat("  DBH           chr9  : 133.6-133.7 Mb\n")
cat("  BDNF          chr11 : 27.6-27.7 Mb\n")
cat("\nAn index variant inside one of these is a variant where SNP -> Smoking is\n")
cat("established, so the collider path of METHODS.md 7.2 is open at that locus.\n")
