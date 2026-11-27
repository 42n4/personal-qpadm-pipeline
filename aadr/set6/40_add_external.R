# ==============================================================================
# 40_add_external.R -- merge an EXTERNAL PLINK package into a set-N fileset
#
# Written for the Poseidon package 2025_Gretzinger_SlavicPeriod (Nature 2025,
# doi:10.1038/s41586-025-09437-6), which publishes genotypes rather than reads:
# PLINK bed/bim/fam, snpSet 1240K, SNP ids identical to the AADR 1240K panel
# (and HO is a strict subset of it), so no liftover or id translation is needed.
# Nothing here is specific to that package beyond the defaults.
#
# What it does
#   1. keeps the requested groups and drops individuals below MIN_SNPS
#      (aDNA packages carry near-empty individuals that would collapse a
#      maxmiss=0 cache)
#   2. merges into <run_dir>/<FROM_SET>/set2_<target>_all with p-link, excluding
#      allele conflicts from BOTH sides on the retry (pitfall 9)
#   3. re-applies the "target-genotyped AND in the panel" SNP filter (pitfall 18)
#   4. writes <run_dir>/set6/set2_<target>_all so the set-5 test script can run
#      on it unchanged (SET=set6)
#
# The external .fam must already carry the group in the FID column, which is the
# Poseidon convention. Individual-level relabelling is NOT done here.
#
# Usage:
#   Rscript aadr/set6/40_add_external.R <run_dir> <aadr_prefix> <target> <external_prefix>
# Env: FROM_SET (default set5), TO_SET (default set6), MIN_SNPS (default 50000),
#      GROUPS (comma-separated; default = the Gretzinger list below)
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: Rscript 40_add_external.R <run_dir> <aadr_prefix> <target> <external_prefix>")
run_dir <- args[1]; aadr_prefix <- args[2]; target <- args[3]; ext <- args[4]

from_set <- Sys.getenv("FROM_SET", "set5")
to_set   <- Sys.getenv("TO_SET", "set6")
min_snps <- as.integer(Sys.getenv("MIN_SNPS", "50000"))
groups   <- strsplit(Sys.getenv("GROUPS",
  "Poland_MP,Poland_SP,NWUkraine_SP,EGermany_MP,EGermany_SP,NWBalkan_SP,NWBalkan_MP,Latvia_MA"), ",")[[1]]

run <- function(cmd, a) { r <- system2(cmd, a, stdout = FALSE, stderr = FALSE)
  if (r != 0L) stop(sprintf("%s failed: %s", cmd, paste(a, collapse = " "))) }
read_fam <- function(p) read.table(paste0(p, ".fam"), stringsAsFactors = FALSE,
                                   col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"))

base <- file.path(run_dir, from_set, paste0("set2_", target, "_all"))
if (!file.exists(paste0(base, ".bed"))) stop("no ", from_set, " fileset at ", base)
out <- file.path(run_dir, to_set); dir.create(out, showWarnings = FALSE, recursive = TRUE)
pre <- file.path(out, paste0("set2_", target))
cat(sprintf("  base: %s\n  external: %s\n  groups: %s\n", base, ext, paste(groups, collapse = ", ")))

# ---- 1. select groups, drop low-coverage individuals -------------------------
efam <- read_fam(ext)
have <- intersect(groups, unique(efam$FID))
if (length(setdiff(groups, have)))
  cat(sprintf("  WARNING: not in the external .fam, skipped: %s\n",
              paste(setdiff(groups, have), collapse = ", ")))
if (!length(have)) stop("none of the requested groups is in ", ext)

miss_stem <- paste0(pre, "_extmiss")
run("plink2", c("--bfile", ext, "--missing", "sample-only", "--out", miss_stem))
m <- read.table(paste0(miss_stem, ".smiss"), header = TRUE, comment.char = "", sep = "\t")
names(m)[1] <- "FID"
m$covered <- m$OBS_CT - m$MISSING_CT
keep <- m[m$FID %in% have & m$covered >= min_snps, c("FID", "IID")]
dropped <- sum(m$FID %in% have) - nrow(keep)
cat(sprintf("  external individuals: %d in the requested groups, %d kept (>= %d SNPs), %d dropped\n",
            sum(m$FID %in% have), nrow(keep), min_snps, dropped))
per <- aggregate(covered ~ FID, data = m[m$FID %in% have, ], FUN = function(x) c(n = length(x), mean = mean(x)))
print(data.frame(group = per$FID, n = per$covered[, "n"], mean_snps = round(per$covered[, "mean"])), row.names = FALSE)
keepf <- paste0(pre, "_extkeep.txt")
write.table(keep, keepf, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
esub <- paste0(pre, "_ext")
run("plink2", c("--bfile", ext, "--keep", keepf, "--chr", "1-22", "--make-bed", "--out", esub))

# ---- 2. merge (pitfall 7: p-link; pitfall 9: exclude conflicts both sides) ----
merged_all <- paste0(pre, "_all")
ret <- system2("p-link", c("--bfile", base, "--bmerge", paste0(esub, c(".bed", ".bim", ".fam")),
                           "--make-bed", "--out", merged_all), stdout = FALSE, stderr = FALSE)
missnp <- paste0(merged_all, "-merge.missnp")
if (ret != 0L && file.exists(missnp)) {
  cat(sprintf("  %d allele conflicts -- excluding from both sides and retrying\n",
              length(readLines(missnp))))
  run("plink2", c("--bfile", base, "--exclude", missnp, "--make-bed", "--out", paste0(merged_all, "_a")))
  run("plink2", c("--bfile", esub, "--exclude", missnp, "--make-bed", "--out", paste0(merged_all, "_b")))
  run("p-link", c("--bfile", paste0(merged_all, "_a"), "--bmerge",
                  paste0(merged_all, "_b", c(".bed", ".bim", ".fam")), "--make-bed", "--out", merged_all))
} else if (ret != 0L) stop("p-link --bmerge failed; see ", merged_all, ".log")

# ---- 3. target-genotyped AND in the panel (pitfall 18) -----------------------
tf <- paste0(pre, "_target.txt"); writeLines(target, tf)
run("plink2", c("--bfile", merged_all, "--keep-fam", tf, "--geno", "0",
                "--write-snplist", "--out", paste0(pre, "_tgt")))
panel_ids <- read.table(paste0(aadr_prefix, ".snp"), colClasses = c("character", rep("NULL", 5)))[[1]]
keep_ids <- intersect(readLines(paste0(pre, "_tgt.snplist")), panel_ids)
writeLines(keep_ids, paste0(pre, "_tgt.snplist"))
final <- paste0(pre, "_all")
run("plink2", c("--bfile", merged_all, "--extract", paste0(pre, "_tgt.snplist"),
                "--chr", "1-22", "--make-bed", "--out", paste0(pre, "_tmp")))
for (e in c(".bed", ".bim", ".fam")) file.rename(paste0(pre, "_tmp", e), paste0(final, e))

fam <- read_fam(final)
cat(sprintf("\n  %s: %d SNPs, %d individuals, %d groups\n", final,
            length(readLines(paste0(final, ".bim"))), nrow(fam), length(unique(fam$FID))))
print(sort(table(fam$FID[fam$FID %in% have]), decreasing = TRUE))
cat(sprintf("\n  Next: SET=%s Rscript aadr/set5/30_trzciniec_tests.R %s %s\n", to_set, run_dir, target))
