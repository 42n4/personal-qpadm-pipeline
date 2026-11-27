# ==============================================================================
# set3/21_region_tests.R -- which medieval Polish region is the target closest to?
#
# Four tests on the set-3 f2 cache (regions with n >= MIN_REGION_N only):
#  A. Region distinctness: qpWave(left = A, B; right = deep + distal outgroups).
#     p < 0.05 = the two regions are genetically distinguishable at all. If most
#     pairs are NOT, no target test below can separate them either.
#  B. Pairwise f4(Mbuti, TARGET; A, B): z <= -3 target closer to A, z >= 3 closer
#     to B. Summed into a per-region score (wins - losses at |z| >= 2).
#  C. Outgroup f3(Mbuti; TARGET, X) for regions and reference periods.
#  D. Rotating qpAdm: pool = regions + reference periods; for each 1-source and
#     2-source model the unused pool members join the right set (as set 2/13).
#     Verdicts FIT / FIT-unphysical / UNIDENTIFIABLE (SE > MAX_SE) / REJECTED.
#
# Usage:  Rscript aadr/set3/21_region_tests.R <merged_prefix> <aadr_prefix> [target]
# Env:    MAX_SE (0.25), MAX_SRC (2)
# Output: <run_dir>/set3/21_{distinct,f4pairs,f3,rotating}_<target>.tsv
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 220)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "regions.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 21_region_tests.R <merged_prefix> <aadr_prefix> [target]")
target  <- if (length(args) >= 3) args[3] else "PW"
max_se  <- as.numeric(Sys.getenv("MAX_SE", "0.25"))
max_src <- as.integer(Sys.getenv("MAX_SRC", "2"))
d <- set3_dirs(args[1], target)
pops <- readLines(paste0(d$prefix, "_pops.txt"))
fam  <- read.table(paste0(d$prefix, ".fam"), stringsAsFactors = FALSE)
n_by <- table(fam[[1]])
regs <- names(REGIONS)[names(REGIONS) %in% pops & names(REGIONS) %in% names(n_by)]
small <- regs[n_by[regs] < MIN_REGION_N]
regs  <- setdiff(regs, small)
periods <- intersect(c(sapply(PERIOD_SPLITS, `[[`, "label"), REF_PERIODS), pops)
right_fixed <- c(RIGHT3, DISTAL3)
f2b <- read_f2(d$f2, pops = pops)
save <- function(x, name) write.table(x, file.path(d$out, sprintf("21_%s_%s.tsv", name, target)),
                                      sep = "\t", quote = FALSE, row.names = FALSE)

cat("================================================================================\n")
cat(sprintf("  SET 3 / 21: %s vs medieval Polish regions\n", target))
cat("================================================================================\n")
cat(sprintf("  regions used (n >= %d): %s\n", MIN_REGION_N,
            paste(sprintf("%s(%d)", regs, as.integer(n_by[regs])), collapse = ", ")))
if (length(small)) cat(sprintf("  too small, f3 only: %s\n", paste(small, collapse = ", ")))
cat(sprintf("  reference periods: %s\n\n", paste(periods, collapse = ", ")))
sp <- paste0(d$prefix, "_sparse.txt")
if (file.exists(sp) && length(readLines(sp)))
  cat(sprintf("  NOT TESTED (too sparse for the joint cache; rerun with ALLSNPS=1): %s\n\n",
              paste(readLines(sp), collapse = ", ")))

# ---- A. distinctness ---------------------------------------------------------
pairs <- combn(regs, 2, simplify = FALSE)
A <- do.call(rbind, lapply(pairs, function(pr) {
  w <- tryCatch(suppressMessages(qpwave(f2b, left = pr, right = right_fixed, verbose = FALSE)),
                error = function(e) NULL)
  p <- if (is.null(w)) NA else w$rankdrop$p[w$rankdrop$f4rank == 0]
  data.frame(A = pr[1], B = pr[2], p = p, distinct = !is.na(p) & p < 0.05)
}))
cat(sprintf("--- A. region pairs distinguishable (qpWave p < 0.05): %d of %d ---\n",
            sum(A$distinct, na.rm = TRUE), nrow(A)))
print(transform(A[order(A$p), ], p = signif(p, 3)), row.names = FALSE)
save(A, "distinct")

# ---- B. pairwise f4 ----------------------------------------------------------
B <- do.call(rbind, lapply(pairs, function(pr) {
  r <- f4(f2b, "Mbuti", target, pr[1], pr[2])
  data.frame(A = pr[1], B = pr[2], f4 = r$est, se = r$se, z = r$z,
             closer = ifelse(r$z <= -2, pr[1], ifelse(r$z >= 2, pr[2], "")))
}))
score <- sapply(regs, function(r) sum(B$closer == r) - sum(B$closer != "" & (B$A == r | B$B == r) & B$closer != r))
cat(sprintf("\n--- B. f4(Mbuti, %s; A, B): pairs with |z| >= 2: %d of %d ---\n", target,
            sum(abs(B$z) >= 2), nrow(B)))
print(transform(B[abs(B$z) >= 2, ], f4 = signif(f4, 3), se = signif(se, 3), z = round(z, 2)), row.names = FALSE)
cat("\n  per-region score (pairs closer - pairs farther, |z| >= 2):\n")
print(sort(score, decreasing = TRUE))
save(B, "f4pairs")

# ---- C. outgroup f3 ----------------------------------------------------------
xs <- setdiff(c(regs, small, periods), target)
C <- f3(f2b, "Mbuti", target, xs)
C <- C[order(-C$est), c("pop3", "est", "se")]
names(C) <- c("population", "f3", "se")
C$n <- as.integer(n_by[C$population])
C$gap_to_top_se <- round((C$f3[1] - C$f3) / sqrt(C$se[1]^2 + C$se^2), 2)
cat(sprintf("\n--- C. outgroup f3(Mbuti; %s, X) ---\n", target))
print(transform(C, f3 = round(f3, 5), se = round(se, 5)), row.names = FALSE)
save(C, "f3")

# ---- D. rotating qpAdm -------------------------------------------------------
pool <- c(regs, periods)
models <- c(lapply(pool, identity), if (max_src >= 2) combn(pool, 2, simplify = FALSE) else list())
rows <- lapply(models, function(left) {
  m <- tryCatch(suppressMessages(qpadm(f2b, left = left, right = c(right_fixed, setdiff(pool, left)),
                                       target = target, verbose = FALSE)), error = function(e) NULL)
  if (is.null(m)) return(data.frame(n_src = length(left), sources = paste(left, collapse = " + "),
                                    p = NA, weights = "", verdict = "ERROR"))
  p <- m$rankdrop$p[m$rankdrop$f4rank == length(left) - 1L]
  w <- m$weights[match(left, m$weights$left), ]
  phys <- all(w$weight >= -2 * w$se & w$weight <= 1 + 2 * w$se)
  data.frame(n_src = length(left), sources = paste(left, collapse = " + "), p = p,
             weights = paste(sprintf("%.2f±%.2f", w$weight, w$se), collapse = " | "),
             verdict = if (p < 0.05) "REJECTED" else if (any(w$se > max_se)) "UNIDENTIFIABLE"
                       else if (phys) "FIT" else "FIT-unphysical")
})
D <- do.call(rbind, rows)
D <- D[order(match(D$verdict, c("FIT", "FIT-unphysical", "UNIDENTIFIABLE", "REJECTED", "ERROR")), D$n_src, -D$p), ]
cat(sprintf("\n--- D. rotating qpAdm, pool = %d (regions + periods) ---\n", length(pool)))
cat("  1-source:\n")
print(transform(D[D$n_src == 1, c("sources", "p", "verdict")], p = signif(p, 3)), row.names = FALSE)
two <- D[D$n_src == 2 & D$verdict %in% c("FIT", "FIT-unphysical"), ]
cat(sprintf("\n  2-source FIT with identifiable weights: %d of %d (UNIDENTIFIABLE: %d)\n",
            nrow(two), sum(D$n_src == 2), sum(D$n_src == 2 & D$verdict == "UNIDENTIFIABLE")))
if (nrow(two)) print(transform(head(two, 25)[, c("sources", "p", "weights")], p = signif(p, 3)), row.names = FALSE)
save(D, "rotating")
cat(sprintf("\n  Saved: %s/21_*_%s.tsv\n", d$out, target))
