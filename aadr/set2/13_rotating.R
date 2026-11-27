# ==============================================================================
# set2/13_rotating.R -- rotating qpAdm: candidates as outgroups for each other
#
# For every 1-source and 2-source model drawn from the candidate pool:
#   left  = the source(s)
#   right = RIGHT_DISTAL + WHG/EHG/EEF/Steppe + every OTHER pool candidate
# (Harney et al. 2021 "rotating" design). A model now fails if the target
# shares drift with any unused candidate beyond what its sources explain, so
# close European sources compete directly -- the property set 1 lacks.
#
# Verdicts:
#   FIT              p >= 0.05, every SE <= MAX_SE and weight within [-2 SE, 1 + 2 SE]
#   FIT-unphysical   p >= 0.05, SEs fine, but a weight outside that range
#   UNIDENTIFIABLE   p >= 0.05 but some SE > MAX_SE (default 0.25): the two
#                    sources are collinear, any weights fit -- not evidence
#   REJECTED         p <  0.05
# A 2-source FIT whose 1-source parts are both REJECTED is real admixture
# evidence; a FIT that only adds a second source to a 1-source FIT is not.
#
# Usage:  Rscript aadr/set2/13_rotating.R <merged_prefix> <aadr_prefix> [target]
# Env:    MAX_SRC=1 to skip the 2-source models (pool of 16 = 120 models).
# Output: <run_dir>/set2/13_rotating_<target>.tsv
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 220)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
# POPS_FILE lets another set (e.g. set4/pops.R) reuse these scripts with its own
# population lists and output directory; default = this set's pops.R.
source(Sys.getenv("POPS_FILE", file.path(here, "pops.R")))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 13_rotating.R <merged_prefix> <aadr_prefix> [target]")
target  <- if (length(args) >= 3) args[3] else "PW"
max_src <- as.integer(Sys.getenv("MAX_SRC", "2"))
max_se  <- as.numeric(Sys.getenv("MAX_SE", "0.25"))
d <- set2_dirs(args[1], target)
pops <- readLines(paste0(d$prefix, "_pops.txt"))
pool <- intersect(CANDIDATES, pops)
fixed_right <- c(RIGHT_DISTAL, WHG, EHG, EEF, STEPPE)
f2b <- read_f2(d$f2, pops = c(target, fixed_right, pool))

models <- c(lapply(pool, identity),
            if (max_src >= 2) combn(pool, 2, simplify = FALSE) else list())
cat("================================================================================\n")
cat(sprintf("  SET 2 / 13: rotating qpAdm, pool %d candidates, %d models\n", length(pool), length(models)))
cat(sprintf("  fixed right (%d): %s\n", length(fixed_right), paste(fixed_right, collapse = ", ")))
cat("================================================================================\n\n")

rows <- vector("list", length(models))
for (i in seq_along(models)) {
  left  <- models[[i]]
  right <- c(fixed_right, setdiff(pool, left))
  m <- tryCatch(suppressMessages(qpadm(f2b, left = left, right = right, target = target, verbose = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) {
    rows[[i]] <- data.frame(n_src = length(left), sources = paste(left, collapse = " + "),
                            p = NA, weights = "", verdict = "ERROR", stringsAsFactors = FALSE)
    next
  }
  p <- m$rankdrop$p[m$rankdrop$f4rank == length(left) - 1L]
  w <- m$weights[match(left, m$weights$left), ]
  physical <- all(w$weight >= -2 * w$se & w$weight <= 1 + 2 * w$se)
  rows[[i]] <- data.frame(
    n_src = length(left), sources = paste(left, collapse = " + "), p = p,
    weights = paste(sprintf("%.2f±%.2f", w$weight, w$se), collapse = " | "),
    verdict = if (p < 0.05) "REJECTED" else if (any(w$se > max_se)) "UNIDENTIFIABLE"
              else if (physical) "FIT" else "FIT-unphysical",
    stringsAsFactors = FALSE)
  if (i %% 20 == 0) cat(sprintf("  %d/%d models\n", i, length(models)))
}
res <- do.call(rbind, rows)

# Flag 2-source fits that need both sources (both 1-source parts rejected).
one <- setNames(res$verdict[res$n_src == 1], res$sources[res$n_src == 1])
res$note <- ""
for (k in which(res$n_src == 2 & res$verdict == "FIT")) {
  parts <- strsplit(res$sources[k], " + ", fixed = TRUE)[[1]]
  if (all(one[parts] == "REJECTED")) res$note[k] <- "both parts rejected alone: ADMIXTURE"
  else res$note[k] <- "a part fits alone"
}
ord <- order(match(res$verdict, c("FIT", "FIT-unphysical", "UNIDENTIFIABLE", "REJECTED", "ERROR")), res$n_src, -res$p)
res <- res[ord, ]

cat("\n--- 1-source models ---\n")
print(transform(res[res$n_src == 1, c("sources", "p", "weights", "verdict")], p = signif(p, 3)), row.names = FALSE)
if (max_src >= 2) {
  two <- res[res$n_src == 2 & res$verdict %in% c("FIT", "FIT-unphysical"), ]
  cat(sprintf("\n--- 2-source models that fit with identifiable weights: %d of %d (%d more fit but are UNIDENTIFIABLE) ---\n",
              nrow(two), sum(res$n_src == 2), sum(res$n_src == 2 & res$verdict == "UNIDENTIFIABLE")))
  if (nrow(two) > 0) print(transform(two[, c("sources", "p", "weights", "verdict", "note")], p = signif(p, 3)), row.names = FALSE)
}
out <- file.path(d$out, paste0("13_rotating_", target, ".tsv"))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Tally: %s\n  Saved: %s\n", paste(names(table(res$verdict)), table(res$verdict), sep = "=", collapse = "  "), out))
