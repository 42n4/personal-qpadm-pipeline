# ==============================================================================
# set2/14_f4_scan.R -- what separates the target from its best sources?
#
# For each of the TOP best 1-source candidates S from script 13 (by p):
#   f4(Mbuti, X; S, TARGET)   for every other population X in the cache
#   z >= +3 : X shares MORE drift with TARGET than with S
#             -> TARGET carries ancestry related to X that S lacks
#   z <= -3 : X is closer to S than to TARGET
#             -> S carries X-related ancestry that TARGET lacks
# |z| < 3 everywhere = S and TARGET symmetric to the whole panel (a clade).
#
# Usage:  Rscript aadr/set2/14_f4_scan.R <merged_prefix> <aadr_prefix> [target]
# Env:    TOP=3 (number of sources), SOURCES="A,B" overrides the choice.
# Output: <run_dir>/set2/14_f4_scan_<target>.tsv
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
# POPS_FILE lets another set (e.g. set4/pops.R) reuse these scripts with its own
# population lists and output directory; default = this set's pops.R.
source(Sys.getenv("POPS_FILE", file.path(here, "pops.R")))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 14_f4_scan.R <merged_prefix> <aadr_prefix> [target]")
target <- if (length(args) >= 3) args[3] else "PW"
d <- set2_dirs(args[1], target)
pops <- readLines(paste0(d$prefix, "_pops.txt"))

srcs <- strsplit(Sys.getenv("SOURCES", ""), ",")[[1]]
if (length(srcs) == 0) {
  rot <- file.path(d$out, paste0("13_rotating_", target, ".tsv"))
  if (!file.exists(rot)) stop("run 13_rotating.R first, or set SOURCES=")
  r <- read.delim(rot, stringsAsFactors = FALSE)
  r <- r[r$n_src == 1 & !is.na(r$p), ]
  srcs <- head(r$sources[order(-r$p)], as.integer(Sys.getenv("TOP", "3")))
}
f2b <- read_f2(d$f2, pops = pops)

cat("================================================================================\n")
cat(sprintf("  SET 2 / 14: f4(%s, X; S, %s) scan, S = %s\n", OUTGROUP, target, paste(srcs, collapse = ", ")))
cat("================================================================================\n")
all_rows <- list()
for (s in srcs) {
  xs <- setdiff(pops, c(OUTGROUP, target, s))
  r <- f4(f2b, OUTGROUP, xs, s, target)
  r <- r[order(-r$z), c("pop2", "est", "se", "z")]
  names(r) <- c("X", "f4", "se", "z")
  r$source <- s
  r$reading <- ifelse(r$z >= 3, paste0("target has extra X-related ancestry vs ", s),
               ifelse(r$z <= -3, paste0(s, " has extra X-related ancestry vs target"), ""))
  sig <- r[abs(r$z) >= 3, ]
  cat(sprintf("\n--- S = %s: %d of %d X with |z| >= 3 ---\n", s, nrow(sig), nrow(r)))
  if (nrow(sig) > 0) print(transform(sig[, c("X", "f4", "se", "z")], f4 = signif(f4, 3), se = signif(se, 3), z = round(z, 2)), row.names = FALSE)
  else cat("  none -- S and target are symmetric to every population in the panel\n")
  all_rows[[s]] <- r
}
out <- file.path(d$out, paste0("14_f4_scan_", target, ".tsv"))
write.table(do.call(rbind, all_rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", out))
