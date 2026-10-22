# ==============================================================================
# set2/10_build_f2.R -- fileset + f2 cache for the set-2 comparisons
#
# Adds every population named in set2/pops.R to the step-03 merge
# (merged_pw_aadr_final_ready), restricts to SNPs the target is genotyped at,
# and computes ONE f2 cache that scripts 11-14 read.
#
# Usage:
#   Rscript aadr/set2/10_build_f2.R <merged_prefix> <aadr_prefix> [target]
#   e.g. Rscript aadr/set2/10_build_f2.R ./PW_WGS_HO/merged_pw_aadr_final_ready \
#          $AADR_GENO/v66.p1_HO.aadr.patch.PUB PW
#
# Environment:
#   MAXMISS  extract_f2 maxmiss (default 0 = SNP must be covered in EVERY
#            population). Ancient singletons (Malta, UstIshim) cut the SNP
#            count; the script prints it -- raise to e.g. 0.1 only knowingly.
#
# Output (in <run_dir>/set2/):
#   set2_<target>.{bed,bim,fam}  fileset, FID = population
#   f2cache_<target>/            f2 blocks for scripts 11-14
#   set2_<target>_pops.txt       the population list actually used
#
# Differences from steps 05-08 (deliberate):
#   - PED written per INDIVIDUAL (vectorised over SNPs), not per SNP: same
#     output, minutes instead of the per-SNP R loop.
#   - "." alleles in the extra fileset are filled from the AADR .snp pair before
#     merging, so p-link cannot mis-flag them as 3+ alleles (pitfall 9).
#   - extract_f2 always runs with maxmem=1000 (block-splitting, same result);
#     unsplit, step 07 reached 43.7 GB on WGS x HO.
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)

here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
# POPS_FILE lets another set (e.g. set4/pops.R) reuse these scripts with its own
# population lists and output directory; default = this set's pops.R.
source(Sys.getenv("POPS_FILE", file.path(here, "pops.R")))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 10_build_f2.R <merged_prefix> <aadr_prefix> [target]")
merged_prefix <- args[1]
aadr_prefix   <- args[2]
target        <- if (length(args) >= 3) args[3] else "PW"
maxmiss       <- as.numeric(Sys.getenv("MAXMISS", "0"))
d <- set2_dirs(merged_prefix, target)

run <- function(cmd, a) {
  ret <- system2(cmd, a, stdout = FALSE, stderr = FALSE)
  if (ret != 0L) stop(sprintf("%s failed (exit %d): %s", cmd, ret, paste(a, collapse = " ")))
}
read_fam <- function(p) read.table(paste0(p, ".fam"), stringsAsFactors = FALSE,
                                   col.names = c("FID","IID","PID","MID","SEX","PHEN"))

cat("================================================================================\n")
cat("  SET 2 / 10: build fileset + f2 cache\n")
cat("================================================================================\n")
cat(sprintf("  merged : %s\n  AADR   : %s\n  target : %s\n  maxmiss: %g\n\n",
            merged_prefix, aadr_prefix, target, maxmiss))

counts   <- panel_counts(aadr_prefix)
required <- c(RIGHT_DISTAL, WHG, EHG, EEF, STEPPE)
missing_req <- setdiff(required, names(counts))
if (length(missing_req) > 0)
  stop("Required reference populations not in this panel: ", paste(missing_req, collapse = ", "))
cands <- available(CANDIDATES, counts, "candidates")
mods  <- available(MODERNS, counts, "moderns")
all_pops <- unique(c(required, cands, mods))
writeLines(c(target, all_pops), paste0(d$prefix, "_pops.txt"))
cat(sprintf("  populations: %d reference, %d candidates, %d moderns\n\n",
            length(required), length(cands), length(mods)))

# ---- 1. fileset (idempotent) -------------------------------------------------
have_all <- file.exists(paste0(d$prefix, ".bed")) &&
            all(c(target, all_pops) %in% read_fam(d$prefix)$FID)
if (have_all) {
  cat("  fileset exists with all populations -- skipping build\n\n")
} else {
  existing <- unique(read_fam(merged_prefix)$FID)
  extra <- setdiff(all_pops, existing)
  cat(sprintf("  already in merge: %d, reading from AADR: %d\n", length(intersect(all_pops, existing)), length(extra)))
  # PITFALL (step 2): never packedancestrymap_to_plink() -- whole panel, OOM.
  gd <- read_packedancestrymap(pref = aadr_prefix, pops = extra, verbose = FALSE)
  geno <- gd$geno; snp <- gd$snp; ind <- gd$ind
  cat(sprintf("  loaded %d SNPs x %d individuals\n", nrow(geno), ncol(geno)))

  stem <- paste0(d$prefix, "_extra")
  write.table(data.frame(snp$CHR, snp$SNP, snp$cm, snp$POS), paste0(stem, ".map"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  a1 <- snp$A1; a2 <- snp$A2
  # geno counts copies of A1 (same convention as steps 02/07).
  con <- file(paste0(stem, ".ped"), "w")
  for (j in seq_len(ncol(geno))) {
    g  <- geno[, j]
    x1 <- ifelse(is.na(g), "0", ifelse(g >= 1L, a1, a2))
    x2 <- ifelse(is.na(g), "0", ifelse(g == 2L, a1, a2))
    writeLines(paste(c(ind[[1]][j], ind[[1]][j], 0, 0, 0, -9,
                       as.vector(rbind(x1, x2))), collapse = "\t"), con)
  }
  close(con)
  rm(geno); invisible(gc())
  run("plink2", c("--pedmap", stem, "--make-bed", "--out", stem))

  # Fill "." alleles from the panel pair (monomorphic / all-missing in this
  # subset); otherwise p-link may flag "." vs allele as 3+ alleles.
  bim <- read.table(paste0(stem, ".bim"), stringsAsFactors = FALSE, colClasses = "character")
  pair <- match(bim$V2, snp$SNP)
  if (anyNA(pair)) stop("extra .bim has IDs absent from the AADR .snp")
  p1 <- a1[pair]; p2 <- a2[pair]
  fill5 <- bim$V5 %in% c(".", "0"); fill6 <- bim$V6 %in% c(".", "0")
  both <- fill5 & fill6
  bim$V6[both] <- p1[both]; bim$V5[both] <- p2[both]
  one5 <- fill5 & !both
  bim$V5[one5] <- ifelse(bim$V6[one5] == p1[one5], p2[one5], p1[one5])
  one6 <- fill6 & !both
  bim$V6[one6] <- ifelse(bim$V5[one6] == p1[one6], p2[one6], p1[one6])
  write.table(bim, paste0(stem, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  filled alleles: %d sites with one, %d with both unobserved\n",
              sum(one5 | one6), sum(both)))

  # Merge with p-link 1.9 (plink2 cannot do non-concatenating merges, pitfall 7).
  merged_all <- paste0(d$prefix, "_all")
  ret <- system2("p-link", c("--bfile", merged_prefix, "--bmerge", paste0(stem, c(".bed", ".bim", ".fam")),
                             "--make-bed", "--out", merged_all), stdout = FALSE, stderr = FALSE)
  missnp <- paste0(merged_all, "-merge.missnp")
  if (ret != 0L && file.exists(missnp)) {
    cat(sprintf("  %d allele conflicts -- excluding from both sides and retrying\n",
                length(readLines(missnp))))
    run("plink2", c("--bfile", merged_prefix, "--exclude", missnp, "--make-bed", "--out", paste0(merged_all, "_a")))
    run("plink2", c("--bfile", stem, "--exclude", missnp, "--make-bed", "--out", paste0(merged_all, "_b")))
    run("p-link", c("--bfile", paste0(merged_all, "_a"), "--bmerge",
                    paste0(merged_all, "_b", c(".bed", ".bim", ".fam")), "--make-bed", "--out", merged_all))
  } else if (ret != 0L) stop("p-link --bmerge failed; see ", merged_all, ".log")

  # FID = population for the added individuals (pitfall 12).
  fam <- read_fam(merged_all)
  lut <- setNames(ind[[3]], ind[[1]])
  hit <- fam$IID %in% names(lut)
  fam$FID[hit] <- lut[fam$IID[hit]]
  write.table(fam, paste0(merged_all, ".fam"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Keep SNPs the target is genotyped at, autosomes only (pitfall 11).
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
absent <- setdiff(c(target, all_pops), fam$FID)
if (length(absent) > 0) stop("populations missing from fileset: ", paste(absent, collapse = ", "))
n_by <- table(fam$FID)

# ---- 2. f2 cache -------------------------------------------------------------
cat("--- extract_f2 (one pass, block-split) ---\n")
extract_f2(d$prefix, outdir = d$f2, pops = c(target, all_pops), maxmiss = maxmiss,
           maxmem = 1000, overwrite = TRUE, format = "plink", verbose = TRUE)

cat("\n  Sample sizes in the cache:\n")
print(n_by[c(target, all_pops)])
cat(sprintf("\n  DONE. f2 cache: %s\n  Next: scripts 11-14 with the same arguments.\n", d$f2))
