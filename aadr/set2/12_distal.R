# ==============================================================================
# set2/12_distal.R -- deep-ancestry proportions: target vs every candidate
#
# Fits the same distal qpAdm models to the target AND to each candidate /
# modern population, against the ancient right set RIGHT_DISTAL:
#   WHG+EHG+EEF     (Loschbour, Karelia HG, Turkey_N)
#   WHG+Steppe+EEF  (Loschbour, Samara Yamnaya, Turkey_N)
# then compares the target's proportions with each population's:
#   dz = (w_target - w_pop) / sqrt(se_target^2 + se_pop^2)   per component.
# A population whose |dz| < 2 on every component is "distally indistinguishable"
# from the target; |dz| >= 3 on any component rules it out as a sole source.
#
# Unlike set 1, the right set here contains lineages that are differentially
# related to WHG/EHG/EEF/steppe, so these proportions are identifiable.
#
# Usage:  Rscript aadr/set2/12_distal.R <merged_prefix> <aadr_prefix> [target]
# Output: <run_dir>/set2/12_distal_<target>.tsv        (fits, all populations)
#         <run_dir>/set2/12_distal_vs_<target>.tsv     (dz comparison table)
# ==============================================================================
suppressPackageStartupMessages(library(admixtools))
options(width = 220)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))))
# POPS_FILE lets another set (e.g. set4/pops.R) reuse these scripts with its own
# population lists and output directory; default = this set's pops.R.
source(Sys.getenv("POPS_FILE", file.path(here, "pops.R")))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 12_distal.R <merged_prefix> <aadr_prefix> [target]")
target <- if (length(args) >= 3) args[3] else "PW"
d <- set2_dirs(args[1], target)
pops <- readLines(paste0(d$prefix, "_pops.txt"))
f2b <- read_f2(d$f2, pops = pops)

subjects <- c(target, intersect(c(CANDIDATES, MODERNS), pops))
rows <- list()
for (mname in names(DISTAL_MODELS)) {
  left <- DISTAL_MODELS[[mname]]
  for (s in subjects) {
    m <- tryCatch(suppressMessages(qpadm(f2b, left = left, right = RIGHT_DISTAL,
                                         target = s, verbose = FALSE)),
                  error = function(e) { cat(sprintf("  ERROR %s / %s: %s\n", mname, s, conditionMessage(e))); NULL })
    if (is.null(m)) next
    p <- m$rankdrop$p[m$rankdrop$f4rank == length(left) - 1L]
    w <- m$weights[match(left, m$weights$left), ]   # component order = model order
    rows[[length(rows) + 1L]] <- data.frame(
      model = mname, population = s, p = p,
      c1 = w$weight[1], c1_se = w$se[1], c2 = w$weight[2], c2_se = w$se[2],
      c3 = w$weight[3], c3_se = w$se[3],
      verdict = if (p < 0.05) "REJECTED" else if (any(w$weight < -2 * w$se)) "OK-negative" else "OK",
      stringsAsFactors = FALSE)
  }
}
fits <- do.call(rbind, rows)
comp_names <- lapply(DISTAL_MODELS, function(x) c("WHG", if (STEPPE %in% x) "Steppe" else "EHG", "EEF"))

cat("================================================================================\n")
cat(sprintf("  SET 2 / 12: distal proportions, right = %d ancient/deep outgroups\n", length(RIGHT_DISTAL)))
cat("================================================================================\n")
cmp_rows <- list()
for (mname in names(DISTAL_MODELS)) {
  cn <- comp_names[[mname]]
  f <- fits[fits$model == mname, ]
  tg <- f[f$population == target, ]
  cat(sprintf("\n--- %s  (components: %s) ---\n", mname, paste(cn, collapse = ", ")))
  if (nrow(tg) == 0) { cat("  target fit failed\n"); next }
  cat(sprintf("  %s: p=%.3f  %s=%.3f±%.3f  %s=%.3f±%.3f  %s=%.3f±%.3f  [%s]\n\n", target, tg$p,
              cn[1], tg$c1, tg$c1_se, cn[2], tg$c2, tg$c2_se, cn[3], tg$c3, tg$c3_se, tg$verdict))
  o <- f[f$population != target, ]
  dz <- function(a, ase, b, bse) (a - b) / sqrt(ase^2 + bse^2)
  o$dz1 <- dz(tg$c1, tg$c1_se, o$c1, o$c1_se)
  o$dz2 <- dz(tg$c2, tg$c2_se, o$c2, o$c2_se)
  o$dz3 <- dz(tg$c3, tg$c3_se, o$c3, o$c3_se)
  o$max_abs_dz <- pmax(abs(o$dz1), abs(o$dz2), abs(o$dz3))
  o$vs_target <- ifelse(o$max_abs_dz < 2, "indistinguishable",
                 ifelse(o$max_abs_dz < 3, "borderline", "different"))
  o <- o[order(o$max_abs_dz), ]
  show <- data.frame(population = o$population, p = round(o$p, 3),
                     x1 = sprintf("%.3f", o$c1), x2 = sprintf("%.3f", o$c2), x3 = sprintf("%.3f", o$c3),
                     dz1 = round(o$dz1, 1), dz2 = round(o$dz2, 1), dz3 = round(o$dz3, 1),
                     fit = o$verdict, vs_target = o$vs_target)
  names(show)[3:8] <- c(cn, paste0("dz_", cn))
  print(show, row.names = FALSE)
  o$model <- mname
  cmp_rows[[mname]] <- o
}

out1 <- file.path(d$out, paste0("12_distal_", target, ".tsv"))
out2 <- file.path(d$out, paste0("12_distal_vs_", target, ".tsv"))
write.table(fits, out1, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, cmp_rows), out2, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  c1..c3 = components in model order (%s / %s).\n",
            paste(comp_names[[1]], collapse = ","), paste(comp_names[[2]], collapse = ",")))
cat(sprintf("  dz = (target - population) / combined SE; |dz| >= 3 on any component rules a sole source out.\n"))
cat(sprintf("  Saved: %s\n         %s\n", out1, out2))
