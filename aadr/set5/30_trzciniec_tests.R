# ==============================================================================
# 30_trzciniec_tests.R -- what must be added to Trzciniec to get the target?
#
# Runs on the set-5 fileset (set 4 + south-eastern candidates) built by
#   POPS_FILE=aadr/set5/pops.R Rscript aadr/set2/10_build_f2.R <merged_prefix> <aadr_prefix> <target>
#
# Relabels before any statistic:
#   - Trzciniec split by region: Gustorzyn (Kuyavia, n=6) vs the rest
#     (Lesser Poland / Swietokrzyskie / Lublin, n=55); "Poland_BA_Trzciniec"
#     stays the pooled label wherever the split is not the question
#   - the small forest-steppe Early Iron Age labels pooled into
#     Ukraine_EIA_ForestSteppe; Ukraine_Classical + the Shushaky Chernyakhov
#     individual pooled into Ukraine_Chernyakhiv_era (same site, 235-350 CE)
#
# Tests (numbering follows info/trzciniec-domieszki-pw-20261023-r2-pl.md, sec. 8):
#   A (8.2) f4(Mbuti, X; Trzciniec, T) for every candidate, computed from the
#           genotypes on the SNPs present in those four groups only (so the
#           low-coverage Vysotska and n=1 groups can be tested too)
#   B (8.6) Gustorzyn vs the rest of Trzciniec: f4(Mbuti, T; Trz_rest, Trz_Gustorzyn)
#           and f4(Mbuti, Trz_Gustorzyn; Trz_rest, T)
#   C (8.3) qpAdm on one joint f2 cache (maxmiss=0): every south-eastern group
#           alone, south-east + Trzciniec, south-east + every other candidate;
#           right = 10 deep outgroups + all unused candidates (rotating)
#   D (8.7) 3-source Trzciniec + {Central-European, south-east} + {Baltic};
#           reported only with all weights > 0 and every SE < 0.1
#
# Usage: Rscript aadr/set5/30_trzciniec_tests.R <run_dir> <target>
#   e.g. Rscript aadr/set5/30_trzciniec_tests.R ./PW_WGS_HO PW
# Output: <run_dir>/set5/30_{f4,gustorzyn,qpadm,qpadm3}_<target>.tsv
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 30_trzciniec_tests.R <run_dir> <target>")
run_dir <- args[1]; target <- args[2]
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "pops.R"))
# SET=set6 runs the same tests on the fileset that set6/40_add_external.R built
# (set 5 + the external Migration/Slavic-period groups).
set_name <- Sys.getenv("SET", "set5")
out <- file.path(run_dir, set_name)
base <- file.path(out, paste0("set2_", target, "_all"))
if (!file.exists(paste0(base, ".bed"))) stop("no fileset at ", base, " -- build ", set_name, " first")

GUSTORZYN <- c("poz687.SG", "poz690.SG", "poz693.SG", "poz694.SG", "poz700.SG", "poz702.SG")
POOL <- c("Ukraine_EIA1_ForestSteppe_Scythian" = "Ukraine_EIA_ForestSteppe",
          "Ukraine_EIA1_ForestSteppe_Scythian-1" = "Ukraine_EIA_ForestSteppe",
          "Ukraine_EIA2_ForestSteppe_Scythian" = "Ukraine_EIA_ForestSteppe",
          "Ukraine_EIA2_BuhDniester_ForestSteppe_Scythian" = "Ukraine_EIA_ForestSteppe",
          "Ukraine_Kyiv_EIA1_ForestSteppe_Scythian-1" = "Ukraine_EIA_ForestSteppe",
          "Ukraine_Classical" = "Ukraine_Chernyakhiv_era",
          "Ukraine_Shushaku_Classical_Chernyakhov" = "Ukraine_Chernyakhiv_era")
SOUTHEAST <- c("Ukraine_EIA_ForestSteppe", "Ukraine_MIA_Scythian", "Ukraine_Chernyakhiv_era",
               "Moldova_Scythian", "Ukraine_LIA_Scythian")
# set6 only: Gretzinger et al. 2025 Migration (MP) and Slavic (SP) period groups.
# Poland_MP is the Polish Iron-Age/Migration link the cremation rite had hidden;
# NWUkraine_SP is Volhynia, i.e. the south-eastern direction with real samples
# instead of Scythian proxies. Absent labels are dropped by present() below.
MIGRATION <- c("Poland_MP", "Poland_SP", "NWUkraine_SP",
               "EGermany_MP", "EGermany_SP", "NWBalkan_SP", "NWBalkan_MP", "Latvia_MA")
F4_ONLY <- c("Ukraine_LBA_Vysotska", "Ukraine_LIA", "Ukraine_Roman_BuhDniester")   # n<=2 or low coverage
CENTRAL <- c("Czechia_LBA", "Czechia_IA_Hallstatt", "Czechia_IA_LaTene",
             "Germany_Magdalenenberg_IA_Hallstatt", "Hungary_LBA_EIA_Kyjatice_Piliny")
BALTIC <- c("Latvia_LBA_FlatGraves", "Lithuania_LBA", "Estonia_LBA_StoneCistGraves", "Estonia_EIA_TarandGraves")
TRZ <- "Poland_BA_Trzciniec"

# ---- relabelled filesets: bed/bim symlinked, only the .fam differs ---------
fam0 <- read.table(paste0(base, ".fam"), stringsAsFactors = FALSE, colClasses = "character")
make_set <- function(name, split_gustorzyn) {
  pre <- file.path(out, paste0("set5_", target, "_", name))
  for (ext in c(".bed", ".bim")) {
    if (file.exists(paste0(pre, ext))) file.remove(paste0(pre, ext))
    file.symlink(normalizePath(paste0(base, ext)), paste0(pre, ext))
  }
  f <- fam0
  hit <- f[[1]] %in% names(POOL); f[[1]][hit] <- POOL[f[[1]][hit]]
  if (split_gustorzyn) {
    trz <- f[[1]] == TRZ
    f[[1]][trz & f[[2]] %in% GUSTORZYN] <- "Trz_Gustorzyn"
    f[[1]][trz & !(f[[2]] %in% GUSTORZYN)] <- "Trz_rest"
  }
  write.table(f, paste0(pre, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  pre
}
pooled <- make_set("pooled", FALSE)
split  <- make_set("split", TRUE)
counts <- table(read.table(paste0(pooled, ".fam"), colClasses = "character")[[1]])
gcount <- table(read.table(paste0(split, ".fam"), colClasses = "character")[[1]])
cat(sprintf("  Trzciniec: %d (Gustorzyn %d, rest %d); %s n=%s\n", counts[TRZ], gcount["Trz_Gustorzyn"],
            gcount["Trz_rest"], paste(SOUTHEAST, collapse = "/"),
            paste(as.integer(counts[SOUTHEAST]), collapse = "/")))
present <- function(p) p[p %in% names(counts)]

# ---- A: f4(Mbuti, X; Trzciniec, T) on per-statistic SNP sets ----------------
cat("\n=== A. f4(Mbuti, X; Trzciniec, ", target, "): z > 0 = X closer to ", target, " than to Trzciniec ===\n", sep = "")
xs <- setdiff(present(c(CANDIDATES, SOUTHEAST, MIGRATION, F4_ONLY, WHG, EHG, EEF, STEPPE)),
              c(TRZ, names(POOL), target, "Mbuti"))
fa <- f4(pooled, "Mbuti", xs, TRZ, target, allsnps = TRUE, verbose = FALSE)
fa <- data.frame(X = fa$pop2, n = as.integer(counts[fa$pop2]), f4 = fa$est, se = fa$se, z = fa$z,
                 stringsAsFactors = FALSE)
fa$reading <- ifelse(fa$z >= 3, paste0("X closer to ", target),
               ifelse(fa$z <= -3, "X closer to Trzciniec", ""))
fa <- fa[order(-fa$z), ]
print(transform(fa, f4 = signif(f4, 3), se = signif(se, 3), z = round(z, 2)), row.names = FALSE)
write.table(fa, file.path(out, paste0("30_f4_", target, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- B: Gustorzyn vs the rest of Trzciniec ----------------------------------
cat("\n=== B. Kuyavian (Gustorzyn) vs the rest of Trzciniec ===\n")
fb <- rbind(f4(split, "Mbuti", target, "Trz_rest", "Trz_Gustorzyn", allsnps = TRUE, verbose = FALSE),
            f4(split, "Mbuti", "Trz_Gustorzyn", "Trz_rest", target, allsnps = TRUE, verbose = FALSE),
            f4(split, "Mbuti", c(WHG, EHG, EEF, STEPPE), "Trz_rest", "Trz_Gustorzyn", allsnps = TRUE, verbose = FALSE))
fb <- data.frame(stat = sprintf("f4(%s, %s; %s, %s)", fb$pop1, fb$pop2, fb$pop3, fb$pop4),
                 f4 = fb$est, se = fb$se, z = fb$z)
print(transform(fb, f4 = signif(f4, 3), se = signif(se, 3), z = round(z, 2)), row.names = FALSE)
cat("  row 1: z > 0 = target shares more drift with Gustorzyn than with the rest of Trzciniec\n")
cat("  row 2: z > 0 = Gustorzyn closer to the target than to the rest of Trzciniec\n")
write.table(fb, file.path(out, paste0("30_gustorzyn_", target, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- C: qpAdm on one joint cache --------------------------------------------
pool_cand <- present(unique(c(setdiff(CANDIDATES, c(names(POOL), F4_ONLY)), SOUTHEAST, MIGRATION)))
pool_cand <- pool_cand[as.integer(counts[pool_cand]) >= 2]
right_fixed <- present(RIGHT_DISTAL)
f2dir <- file.path(out, paste0("f2cache30_", target))
extract_f2(pooled, outdir = f2dir, pops = c(target, right_fixed, pool_cand), maxmiss = 0,
           maxmem = 1000, overwrite = TRUE, format = "plink", verbose = FALSE)
f2b <- read_f2(f2dir)
cat(sprintf("\n=== C. qpAdm, joint cache: %d populations, %s SNPs ===\n", length(pool_cand) + length(right_fixed) + 1,
            format(sum(readRDS(file.path(f2dir, "block_lengths_f2.rds"))), big.mark = ",")))
MAX_SE_FIT <- 0.25
fit <- function(left, rotate = TRUE) {
  right <- c(right_fixed, if (rotate) setdiff(pool_cand, left))
  m <- tryCatch(suppressMessages(qpadm(f2b, left = left, right = right, target = target, verbose = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) return(data.frame(n_src = length(left), sources = paste(left, collapse = " + "), p = NA,
                                    weights = "", min_w = NA, max_se = NA, verdict = "ERROR"))
  p <- m$rankdrop$p[m$rankdrop$f4rank == length(left) - 1L]
  w <- m$weights[match(left, m$weights$left), ]
  verdict <- if (p < 0.05) "REJECTED" else if (any(w$weight <= 0) || any(w$se > MAX_SE_FIT)) "UNIDENTIFIABLE" else "FIT"
  data.frame(n_src = length(left), sources = paste(left, collapse = " + "), p = p,
             weights = paste(sprintf("%.2f±%.2f", w$weight, w$se), collapse = " | "),
             min_w = min(w$weight), max_se = max(w$se), verdict = verdict, stringsAsFactors = FALSE)
}
se_in <- intersect(c(SOUTHEAST, MIGRATION), pool_cand)
models <- c(lapply(c(TRZ, se_in), function(s) s),
            lapply(se_in, function(s) c(TRZ, s)),
            unlist(lapply(se_in, function(s) lapply(setdiff(pool_cand, c(s, TRZ)), function(o) c(s, o))), recursive = FALSE))
resC <- do.call(rbind, lapply(models, fit))
resC <- resC[order(resC$n_src, -resC$p), ]
print(transform(resC[, c("sources", "p", "weights", "verdict")], p = signif(p, 3)), row.names = FALSE)
write.table(resC, file.path(out, paste0("30_qpadm_", target, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- D: 3-source models -------------------------------------------------------
cat("\n=== D. 3-source: Trzciniec + {Central-EU or south-east} + {Baltic}; all w > 0, SE < 0.1 ===\n")
mid <- intersect(c(CENTRAL, se_in), pool_cand); bal <- intersect(BALTIC, pool_cand)
resD <- do.call(rbind, lapply(mid, function(m) do.call(rbind, lapply(bal, function(b) fit(c(TRZ, m, b))))))
resD$strict <- with(resD, verdict == "FIT" & min_w > 0 & max_se < 0.1)
resD <- resD[order(-resD$p), ]
print(transform(resD[, c("sources", "p", "weights", "verdict", "strict")], p = signif(p, 3)), row.names = FALSE)
cat(sprintf("  strict fits: %d of %d\n", sum(resD$strict, na.rm = TRUE), nrow(resD)))
write.table(resD, file.path(out, paste0("30_qpadm3_", target, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s/30_*_%s.tsv\n", out, target))
