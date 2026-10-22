#!/usr/bin/env Rscript
# ==============================================================================
# 07b_continental_f3.R -- is the target admixed with a non-European continent?
#
# Admixture f3(T; E, X) for E = modern West Europeans and X = one continental
# proxy each (Yoruba = Africa, Han = East Asia, Karitiana = Americas,
# Papuan = Oceania). A NEGATIVE estimate with z < -3 proves that T is a mix of
# something E-like and something X-like. All positive means no detectable
# admixture with that continent, and the Europe-only left panel of steps 07/08
# is adequate. If any f3 is significantly negative, add that continent to the
# LEFT (sources) and keep deep outgroups on the right.
#
# Runs on the f2 cache step 07 built (it holds the EU moderns and the five
# right pops); nothing is recomputed.
#
# Usage (after step 07):
#   Rscript aadr_eu/07b_continental_f3.R <merged_prefix> <target>
#   e.g. Rscript aadr_eu/07b_continental_f3.R ./AM_FF/merged_pw_aadr_final_ready AM
# Output: <run_dir>/continental_f3_<target>.tsv
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 07b_continental_f3.R <merged_prefix> <target>")
outdir <- dirname(args[1]); target <- args[2]
f2_dir <- file.path(outdir, "merged_pw_aadr_ext7_f2cache")
if (!dir.exists(f2_dir)) stop("no step-07 f2 cache at ", f2_dir, " -- run 07_extended_models.R first")

EUROPEANS <- c("English", "French", "Spanish", "Italian_North", "Orcadian", "Icelandic")
CONTINENT <- c(Yoruba = "Africa", Han = "East Asia", Karitiana = "Americas", Papuan = "Oceania")

have <- dimnames(read_f2(f2_dir))[[1]]
eu <- intersect(EUROPEANS, have); x <- intersect(names(CONTINENT), have)
if (!target %in% have) stop(target, " is not in the f2 cache")
if (!length(eu) || !length(x)) stop("f2 cache lacks the European or continental populations")
missing <- setdiff(c(EUROPEANS, names(CONTINENT)), have)
if (length(missing)) cat("  not in the cache (skipped):", paste(missing, collapse = ", "), "\n")

f2b <- read_f2(f2_dir, pops = c(target, eu, x))
r <- qp3pop(f2b, target, eu, x)
r <- data.frame(E = r$pop2, X = r$pop3, continent = CONTINENT[r$pop3],
                f3 = signif(r$est, 3), se = signif(r$se, 3), z = round(r$z, 2))
r$verdict <- ifelse(r$z < -3, "ADMIXED", ifelse(r$f3 < 0, "negative, n.s.", "no admixture"))
r <- r[order(r$z), ]
cat(sprintf("\n  Admixture f3(%s; E, X), negative z < -3 = %s is a mix of E and X\n\n", target, target))
print(r, row.names = FALSE)
cat(sprintf("\n  most negative z: %.2f (%s, %s)\n", r$z[1], r$E[1], r$X[1]))
out <- file.path(outdir, paste0("continental_f3_", target, ".tsv"))
write.table(r, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved:", out, "\n")
