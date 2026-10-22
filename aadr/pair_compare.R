# ==============================================================================
# pair_compare.R -- direct comparison of two personal genomes on shared SNPs
#
# Adds a second target (a step-01 pw_plink fileset) to an existing set-3 or
# set-4 fileset of the first target, keeps SNPs genotyped in BOTH targets and in
# the AADR panel, drops groups covering < 60 % of those SNPs, builds one f2
# cache, then:
#   1. qpWave(left = T1, T2; right = deep outgroups [+ all other groups]):
#      p < 0.05 = the two people are distinguishable against that right set
#   2. f4(Mbuti, X; T1, T2) for every group X:
#      z <= -3: X closer to T1;  z >= 3: X closer to T2
#   3. distal WHG/Steppe/EEF proportions of both on the same SNPs
#
# Caveat: T1 WGS vs T2 array means different genotyping error profiles; using
# only SNPs called in both keeps ascertainment identical, but a platform effect
# can still add noise to f4 involving T1 vs T2.
#
# Usage:
#   Rscript aadr/pair_compare.R <set_fileset_prefix> <T1> <T2_plink_prefix> <T2> <aadr_prefix> <outdir>
# e.g.
#   Rscript aadr/pair_compare.R PW_WGS_HO/set3/set3_PW PW \
#     ./JW_MH/pw_plink JW $AADR_GENO/v66.p1_HO.aadr.patch.PUB \
#     PW_WGS_HO/pair_JW_set3
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Usage: Rscript pair_compare.R <set_prefix> <T1> <T2_plink> <T2> <aadr_prefix> <outdir>")
base <- args[1]; t1 <- args[2]; t2p <- args[3]; t2 <- args[4]; aadr <- args[5]; out <- args[6]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
run <- function(cmd, a) { r <- system2(cmd, a, stdout = FALSE, stderr = FALSE)
  if (r != 0L) stop(sprintf("%s failed: %s", cmd, paste(a, collapse = " "))) }
pre <- file.path(out, "pair")

RIGHT_DEEP <- c("Mbuti", "Russia_UstIshim_IUP", "Russia_Kostenki_UP", "Han", "Papuan", "Karitiana",
                "Israel_Natufian", "Iran_GanjDareh_N", "Georgia_KotiasKlde_Mesolithic", "Morocco_Iberomaurusian")
DISTAL <- c("Luxembourg_Loschbour_Mesolithic", "Russia_Samara_EBA_Yamnaya", "Turkey_N")

# ---- SNPs: genotyped in T2, in the panel, present in the base fileset --------
panel <- read.table(paste0(aadr, ".snp"), colClasses = c("character", rep("NULL", 5)))[[1]]
writeLines(t2, paste0(pre, "_t2.txt"))
run("plink2", c("--bfile", t2p, "--geno", "0", "--write-snplist", "--out", paste0(pre, "_t2")))
t2snps <- intersect(readLines(paste0(pre, "_t2.snplist")), panel)
writeLines(t2snps, paste0(pre, "_keep.snplist"))
cat(sprintf("  %s genotyped panel SNPs: %d\n", t2, length(t2snps)))
run("plink2", c("--bfile", base, "--extract", paste0(pre, "_keep.snplist"), "--make-bed", "--out", paste0(pre, "_a")))
run("plink2", c("--bfile", t2p, "--extract", paste0(pre, "_keep.snplist"), "--make-bed", "--out", paste0(pre, "_b")))
r <- system2("p-link", c("--bfile", paste0(pre, "_a"), "--bmerge", paste0(pre, "_b", c(".bed", ".bim", ".fam")),
                         "--make-bed", "--out", pre), stdout = FALSE, stderr = FALSE)
miss <- paste0(pre, "-merge.missnp")
if (r != 0L && file.exists(miss)) {
  cat(sprintf("  %d allele conflicts excluded\n", length(readLines(miss))))
  run("plink2", c("--bfile", paste0(pre, "_a"), "--exclude", miss, "--make-bed", "--out", paste0(pre, "_a2")))
  run("plink2", c("--bfile", paste0(pre, "_b"), "--exclude", miss, "--make-bed", "--out", paste0(pre, "_b2")))
  run("p-link", c("--bfile", paste0(pre, "_a2"), "--bmerge", paste0(pre, "_b2", c(".bed", ".bim", ".fam")),
                  "--make-bed", "--out", pre))
} else if (r != 0L) stop("merge failed, see ", pre, ".log")

# ---- coverage per label, drop sparse ------------------------------------------
fam <- read.table(paste0(pre, ".fam"), stringsAsFactors = FALSE)
labels <- unique(fam[[1]])
nsnp <- length(readLines(paste0(pre, ".bim")))
cov <- sapply(labels, function(g) { writeLines(g, paste0(pre, "_k.txt"))
  run("plink2", c("--bfile", pre, "--keep-fam", paste0(pre, "_k.txt"), "--geno", "0.999999",
                  "--write-snplist", "--out", paste0(pre, "_c")))
  length(readLines(paste0(pre, "_c.snplist"))) / nsnp })
sparse <- setdiff(labels[cov < 0.6], c(t1, t2, RIGHT_DEEP, DISTAL))
use <- setdiff(labels, sparse)
cat(sprintf("  shared SNPs: %d; groups used %d; dropped as sparse: %s\n", nsnp, length(use),
            if (length(sparse)) paste(sparse, collapse = ", ") else "none"))
extract_f2(pre, outdir = paste0(pre, "_f2"), pops = use, maxmiss = 0, maxmem = 1000,
           overwrite = TRUE, format = "plink", verbose = FALSE)
f2b <- read_f2(paste0(pre, "_f2"), pops = use)
others <- setdiff(use, c(t1, t2))

cat("\n================================================================================\n")
cat(sprintf("  %s vs %s\n", t1, t2))
cat("================================================================================\n")
# ---- 1. qpWave ------------------------------------------------------------------
qw <- function(right) { w <- suppressMessages(qpwave(f2b, left = c(t1, t2), right = right, verbose = FALSE))
  w$rankdrop$p[w$rankdrop$f4rank == 0] }
p_deep <- qw(RIGHT_DEEP); p_all <- qw(others)
cat(sprintf("\n1. qpWave rank 0 (are %s and %s a clade?)\n", t1, t2))
cat(sprintf("   right = 10 deep outgroups:        p = %.3g\n", p_deep))
cat(sprintf("   right = all %d other groups:      p = %.3g\n", length(others), p_all))

# ---- 2. f4 scan ------------------------------------------------------------------
x <- setdiff(others, "Mbuti")
f <- f4(f2b, "Mbuti", x, t1, t2)
f <- f[order(f$z), c("pop2", "est", "se", "z")]; names(f) <- c("X", "f4", "se", "z")
f$closer_to <- ifelse(f$z <= -2, t1, ifelse(f$z >= 2, t2, ""))
cat(sprintf("\n2. f4(Mbuti, X; %s, %s): negative = X closer to %s, positive = closer to %s\n", t1, t2, t1, t2))
print(transform(f, f4 = signif(f4, 3), se = signif(se, 3), z = round(z, 2)), row.names = FALSE)
cat(sprintf("   |z| >= 3: %d of %d;  |z| >= 2: %d\n", sum(abs(f$z) >= 3), nrow(f), sum(abs(f$z) >= 2)))
write.table(f, file.path(out, "pair_f4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 3. distal ---------------------------------------------------------------------
cat("\n3. distal WHG + Steppe + EEF on the shared SNPs (right = deep outgroups)\n")
d <- lapply(c(t1, t2), function(t) { m <- suppressMessages(qpadm(f2b, left = DISTAL, right = RIGHT_DEEP, target = t, verbose = FALSE))
  w <- m$weights[match(DISTAL, m$weights$left), ]
  data.frame(target = t, p = m$rankdrop$p[m$rankdrop$f4rank == 2], WHG = w$weight[1], WHG_se = w$se[1],
             Steppe = w$weight[2], Steppe_se = w$se[2], EEF = w$weight[3], EEF_se = w$se[3]) })
d <- do.call(rbind, d); print(transform(d, p = signif(p, 3)), row.names = FALSE, digits = 3)
write.table(d, file.path(out, "pair_distal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s/pair_*.tsv\n", out))
