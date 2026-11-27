# ==============================================================================
# set3/20_build_f2.R -- fileset + f2 cache with medieval Polish REGIONS
#
# Reads the individuals listed by set3/regions.R from AADR, labels each with its
# region (FID = Reg_*), adds the reference periods and deep outgroups, merges
# with the step-03 fileset, and computes one f2 cache for script 21.
#
# Usage:  Rscript aadr/set3/20_build_f2.R <merged_prefix> <aadr_prefix> [target]
# Env:    MAXMISS (default 0), SPARSE_FRAC (0.6): labels covering less than this
#         fraction of the target's SNPs are dropped from the joint cache (strict);
#         ALLSNPS=1 keeps every label and uses pairwise overlap (maxmiss=1).
#         Results of the two modes go to set3/ and set3_allsnps/ respectively.
# Output: <run_dir>/set3/set3_<target>.{bed,bim,fam}, f2cache_<target>/,
#         set3_<target>_pops.txt, set3_region_members.tsv
#
# Same conversion and merge safeguards as set2/10_build_f2.R (vectorised PED,
# "." alleles filled from the AADR pair, p-link merge with conflict retry,
# target-genotyped SNPs only, autosomes, maxmem=1000).
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "regions.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 20_build_f2.R <merged_prefix> <aadr_prefix> [target]")
merged_prefix <- args[1]; aadr_prefix <- args[2]
target  <- if (length(args) >= 3) args[3] else "PW"
maxmiss <- as.numeric(Sys.getenv("MAXMISS", "0"))
d <- set3_dirs(merged_prefix, target)

run <- function(cmd, a) {
  ret <- system2(cmd, a, stdout = FALSE, stderr = FALSE)
  if (ret != 0L) stop(sprintf("%s failed (exit %d): %s", cmd, ret, paste(a, collapse = " ")))
}
read_fam <- function(p) read.table(paste0(p, ".fam"), stringsAsFactors = FALSE,
                                   col.names = c("FID","IID","PID","MID","SEX","PHEN"))

cat("================================================================================\n")
cat("  SET 3 / 20: medieval Polish regions -- fileset + f2 cache\n")
cat("================================================================================\n\n")

mem <- region_members(aadr_prefix)
write.table(mem, file.path(d$out, "set3_region_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Region sizes:\n"); print(table(mem$region)); cat("\n")

ind <- read.table(paste0(aadr_prefix, ".ind"), stringsAsFactors = FALSE)
refs <- c(RIGHT3, DISTAL3, REF_PERIODS)
miss <- setdiff(refs, ind[[3]])
if (any(miss %in% c(RIGHT3, DISTAL3))) stop("required populations missing: ", paste(intersect(miss, c(RIGHT3, DISTAL3)), collapse = ", "))
if (length(miss) > 0) cat(sprintf("  WARNING: reference periods not in panel, skipped: %s\n", paste(miss, collapse = ", ")))
refs <- setdiff(refs, miss)

# label lookup: region for regional individuals, AADR population for the rest
ref_ids <- ind[[1]][ind[[3]] %in% refs]
per <- period_members(aadr_prefix)
cat("  Period splits:\n"); print(table(per$label)); cat("\n")
label <- c(setNames(mem$region, mem$IID), setNames(per$label, per$IID),
           setNames(ind[[3]][ind[[3]] %in% refs], ref_ids))
all_labels <- unique(c(target, names(REGIONS)[names(REGIONS) %in% mem$region],
                       unique(per$label), refs))
writeLines(all_labels, paste0(d$prefix, "_pops.txt"))

if (file.exists(paste0(d$prefix, ".bed")) && all(all_labels %in% read_fam(d$prefix)$FID)) {
  cat("  fileset exists with all labels -- skipping build\n\n")
} else {
  # Keep only the target from the step-03 merge; every reference individual is
  # re-read from AADR below with its region/population label, so nothing is
  # duplicated (step 02 had already added Mbuti/Yoruba/Han etc.).
  keep_tgt <- paste0(d$prefix, "_keep_target.txt")
  writeLines(paste(target, target), keep_tgt)
  base <- paste0(d$prefix, "_target_only")
  run("plink2", c("--bfile", merged_prefix, "--keep", keep_tgt, "--make-bed", "--out", base))

  want <- names(label)
  gd <- read_packedancestrymap(pref = aadr_prefix, inds = want, verbose = FALSE)
  geno <- gd$geno; snp <- gd$snp; gi <- gd$ind
  cat(sprintf("  loaded %d SNPs x %d individuals\n", nrow(geno), ncol(geno)))
  stem <- paste0(d$prefix, "_extra")
  write.table(data.frame(snp$CHR, snp$SNP, snp$cm, snp$POS), paste0(stem, ".map"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  a1 <- snp$A1; a2 <- snp$A2
  con <- file(paste0(stem, ".ped"), "w")
  for (j in seq_len(ncol(geno))) {
    g  <- geno[, j]
    x1 <- ifelse(is.na(g), "0", ifelse(g >= 1L, a1, a2))
    x2 <- ifelse(is.na(g), "0", ifelse(g == 2L, a1, a2))
    writeLines(paste(c(label[[gi[[1]][j]]], gi[[1]][j], 0, 0, 0, -9,
                       as.vector(rbind(x1, x2))), collapse = "\t"), con)
  }
  close(con); rm(geno); invisible(gc())
  run("plink2", c("--pedmap", stem, "--make-bed", "--out", stem))

  bim <- read.table(paste0(stem, ".bim"), stringsAsFactors = FALSE, colClasses = "character")
  pair <- match(bim$V2, snp$SNP); p1 <- a1[pair]; p2 <- a2[pair]
  f5 <- bim$V5 %in% c(".", "0"); f6 <- bim$V6 %in% c(".", "0"); both <- f5 & f6
  bim$V6[both] <- p1[both]; bim$V5[both] <- p2[both]
  o5 <- f5 & !both; bim$V5[o5] <- ifelse(bim$V6[o5] == p1[o5], p2[o5], p1[o5])
  o6 <- f6 & !both; bim$V6[o6] <- ifelse(bim$V5[o6] == p1[o6], p2[o6], p1[o6])
  write.table(bim, paste0(stem, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  merged_all <- paste0(d$prefix, "_all")
  ret <- system2("p-link", c("--bfile", base, "--bmerge", paste0(stem, c(".bed", ".bim", ".fam")),
                             "--make-bed", "--out", merged_all), stdout = FALSE, stderr = FALSE)
  missnp <- paste0(merged_all, "-merge.missnp")
  if (ret != 0L && file.exists(missnp)) {
    cat(sprintf("  %d allele conflicts -- excluding and retrying\n", length(readLines(missnp))))
    run("plink2", c("--bfile", base, "--exclude", missnp, "--make-bed", "--out", paste0(merged_all, "_a")))
    run("plink2", c("--bfile", stem, "--exclude", missnp, "--make-bed", "--out", paste0(merged_all, "_b")))
    run("p-link", c("--bfile", paste0(merged_all, "_a"), "--bmerge",
                    paste0(merged_all, "_b", c(".bed", ".bim", ".fam")), "--make-bed", "--out", merged_all))
  } else if (ret != 0L) stop("p-link --bmerge failed; see ", merged_all, ".log")

  keepf <- paste0(d$prefix, "_target.txt"); writeLines(target, keepf)
  run("plink2", c("--bfile", merged_all, "--keep-fam", keepf, "--geno", "0",
                  "--write-snplist", "--out", paste0(d$prefix, "_tgt")))
  # Keep SNPs that are BOTH target-genotyped AND in the AADR panel. A consumer
  # array is genotyped at ~500K non-panel sites where every AADR group is missing;
  # keeping them left JW with 4,950 usable SNPs in set 2 (2026-10-22).
  panel_ids <- read.table(paste0(aadr_prefix, ".snp"), colClasses = c("character", rep("NULL", 5)))[[1]]
  tgt_ids <- readLines(paste0(d$prefix, "_tgt.snplist"))
  keep_ids <- intersect(tgt_ids, panel_ids)
  cat(sprintf("  target-genotyped SNPs: %d, of which in the panel: %d\n", length(tgt_ids), length(keep_ids)))
  writeLines(keep_ids, paste0(d$prefix, "_tgt.snplist"))
  run("plink2", c("--bfile", merged_all, "--extract", paste0(d$prefix, "_tgt.snplist"),
                  "--chr", "1-22", "--make-bed", "--out", d$prefix))
  cat(sprintf("  fileset: %s (%d SNPs, %d individuals)\n\n", d$prefix,
              length(readLines(paste0(d$prefix, ".bim"))), nrow(read_fam(d$prefix))))
}

fam <- read_fam(d$prefix)
absent <- setdiff(all_labels, fam$FID)
if (length(absent) > 0) stop("labels missing from fileset: ", paste(absent, collapse = ", "))

# ---- coverage per label -------------------------------------------------------
# Pseudo-haploid low-coverage groups cover different random SNP subsets; under
# maxmiss=0 their intersection collapses (first run: 1,304 SNPs on WGS x HO,
# driven by Per_Poland_BA_late_n3 at 52K covered SNPs).
nsnp <- length(readLines(paste0(d$prefix, ".bim")))
covf <- paste0(d$prefix, "_cov")
cov <- sapply(all_labels, function(g) {
  writeLines(g, paste0(covf, ".keep"))
  run("plink2", c("--bfile", d$prefix, "--keep-fam", paste0(covf, ".keep"),
                  "--geno", "0.999999", "--write-snplist", "--out", covf))
  length(readLines(paste0(covf, ".snplist")))
})
cov_df <- data.frame(label = all_labels, n = as.integer(table(fam$FID)[all_labels]),
                     covered = cov, frac = round(cov / nsnp, 3))
write.table(cov_df, paste0(d$prefix, "_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
allsnps <- Sys.getenv("ALLSNPS", "0") == "1"
sparse_frac <- as.numeric(Sys.getenv("SPARSE_FRAC", "0.6"))
protected <- c(target, RIGHT3, DISTAL3)
sparse <- setdiff(cov_df$label[cov_df$frac < sparse_frac], protected)
if (any(cov_df$frac[cov_df$label %in% protected] < sparse_frac))
  cat("  WARNING: a required outgroup/distal source is itself sparse -- kept anyway\n")
cat("  Coverage (fraction of the target's SNPs with >= 1 genotype):\n")
print(cov_df[order(cov_df$frac), c("label", "n", "frac")], row.names = FALSE)

if (allsnps) {
  use <- all_labels; mm <- 1
  cat("\n  ALLSNPS=1: every label kept; each f2 uses the SNPs that pair shares (maxmiss=1)\n")
} else {
  use <- setdiff(all_labels, sparse); mm <- maxmiss
  cat(sprintf("\n  strict mode: dropped as sparse (< %.0f%%): %s\n", 100 * sparse_frac,
              if (length(sparse)) paste(sparse, collapse = ", ") else "none"))
}
writeLines(use, paste0(d$prefix, "_pops.txt"))
writeLines(if (allsnps) character(0) else sparse, paste0(d$prefix, "_sparse.txt"))

cat("--- extract_f2 (one pass, block-split) ---\n")
extract_f2(d$prefix, outdir = d$f2, pops = use, maxmiss = mm,
           maxmem = 1000, overwrite = TRUE, format = "plink", verbose = TRUE)
cat("\n  Sample sizes:\n"); print(table(fam$FID)[use])
cat(sprintf("\n  DONE. f2 cache: %s\n", d$f2))
