# ==============================================================================
# set2/11_outgroup_f3.R -- who shares the most drift with the target?
#
# f3(Mbuti; TARGET, X) for every population in the set-2 cache. Larger f3 =
# more shared genetic drift with TARGET since the split from Mbuti. This is a
# similarity RANKING, not a model test: it needs no right set at all, so it is
# immune to the outgroup-choice problem of set 1.
#
# Read with care: outgroup f3 measures the branch TARGET and X share, so X's
# own later drift does not bias it, but small or low-coverage ancient groups
# get wide SEs. Treat populations within ~2 SE of each other as ties (column
# gap_to_top_se); modern groups usually top this list simply by being recent.
#
# Usage: Rscript aadr/set2/11_outgroup_f3.R <merged_prefix> <aadr_prefix> [target]
# Output: <run_dir>/set2/11_outgroup_f3_<target>.tsv
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 200)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
# POPS_FILE lets another set (e.g. set4/pops.R) reuse these scripts with its own
# population lists and output directory; default = this set's pops.R.
source(Sys.getenv("POPS_FILE", file.path(here, "pops.R")))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 11_outgroup_f3.R <merged_prefix> <aadr_prefix> [target]")
target <- if (length(args) >= 3) args[3] else "PW"
d <- set2_dirs(args[1], target)
pops <- readLines(paste0(d$prefix, "_pops.txt"))
others <- setdiff(pops, c(target, OUTGROUP))

f2b <- read_f2(d$f2, pops = pops)
r <- f3(f2b, OUTGROUP, target, others)
r <- r[order(-r$est), c("pop3", "est", "se", "z")]
names(r) <- c("population", "f3", "se", "z")
r$kind <- ifelse(r$population %in% CANDIDATES, "ancient candidate",
          ifelse(r$population %in% MODERNS, "modern",
          ifelse(r$population %in% c(WHG, EHG, EEF, STEPPE), "distal source", "outgroup")))
r$gap_to_top_se <- round((r$f3[1] - r$f3) / sqrt(r$se[1]^2 + r$se^2), 2)

cat("================================================================================\n")
cat(sprintf("  SET 2 / 11: outgroup f3(%s; %s, X) -- higher = closer\n", OUTGROUP, target))
cat("================================================================================\n\n")
print(transform(r, f3 = round(f3, 5), se = round(se, 5), z = round(z, 1)), row.names = FALSE)
out <- file.path(d$out, paste0("11_outgroup_f3_", target, ".tsv"))
write.table(r, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  gap_to_top_se < 2 = statistically tied with the top population.\n  Saved: %s\n", out))
