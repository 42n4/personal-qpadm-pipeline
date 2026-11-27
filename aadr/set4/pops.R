# ==============================================================================
# set4/pops.R -- PRE-MEDIEVAL comparison: no medieval or modern population at all
#
# Sets 2/3 let medieval groups compete; PW then simply matches medieval Poland.
# Set 4 removes every medieval and modern group from BOTH the candidate pool
# and the right set, so the question becomes: which Bronze Age / Iron Age /
# Roman-period population -- or mix of two -- explains the target?
# Runs the set-2 scripts with POPS_FILE=aadr/set4/pops.R (output in set4/).
#
# AADR v66 has NO Lusatian-culture genomes from Poland, and none of
# Zarubintsy / Kyiv culture / Pomeranian face-urn / Milograd (checked
# 2026-10-22); the pool is the nearest LBA/IA/Roman ring around Poland.
# All labels present in HO and 1240K.
# ==============================================================================

SET_NAME <- "set4"
OUTGROUP <- "Mbuti"

# Deep outgroups (set-2 list without Russia_Malta_UP: n=1, ~1x, it alone cuts
# most SNPs under maxmiss=0).
RIGHT_DISTAL <- c("Mbuti", "Russia_UstIshim_IUP", "Russia_Kostenki_UP", "Han", "Papuan",
                  "Karitiana", "Israel_Natufian", "Iran_GanjDareh_N",
                  "Georgia_KotiasKlde_Mesolithic", "Morocco_Iberomaurusian")

WHG    <- "Luxembourg_Loschbour_Mesolithic"
EHG    <- "Russia_Karelia_Mesolithic_HG"
EEF    <- "Turkey_N"
STEPPE <- "Russia_Samara_EBA_Yamnaya"
DISTAL_MODELS <- list(
  "WHG+EHG+EEF"    = c(WHG, EHG, EEF),
  "WHG+Steppe+EEF" = c(WHG, STEPPE, EEF)
)

CANDIDATES <- c(
  # Poland
  "Poland_CordedWare",                    # 2900-2300 BCE
  "Poland_BA",                            # mostly Early BA 2300-1700 BCE
  "Poland_BA_Trzciniec",                  # 1500-1300 BCE
  "Poland_IA_Wielbark",                   # 200 BCE - 200 CE
  "Poland_Roman_Wielbark",
  "Poland_IA_Przeworsk",                  # n=2
  # Baltic
  "Lithuania_LBA",                        # Turlojiškė, 1250-500 BCE
  "Latvia_LBA_FlatGraves",
  "Estonia_LBA_StoneCistGraves",
  "Estonia_EIA_TarandGraves",
  # East
  "Ukraine_IA_Lusatian",                  # Rovantsi, Volhynia, 1000-700 BCE
  "Ukraine_IA_Hallstatt_Thracian",
  "Ukraine_LIA_Scythian",
  "Russia_Roman_RyazanOka",
  # South
  "Czechia_LBA",                          # Bohemia 1200-800 BCE
  "Czechia_IA_Hallstatt",
  "Czechia_IA_LaTene",
  "Slovakia_EIA_Vekerzug",
  "Hungary_LBA_EIA_Kyjatice_Piliny",      # n=2
  # North / West
  "Sweden_IA",
  "Denmark_Simonsborg_EarlyRoman",
  "Germany_Magdalenenberg_IA_Hallstatt"
)

MODERNS <- character(0)

MIN_N <- 2L

# ---- helpers -----------------------------------------------------------------

# Read .ind, return available populations and counts.
panel_counts <- function(aadr_prefix) {
  ind <- read.table(paste0(aadr_prefix, ".ind"), stringsAsFactors = FALSE)
  table(ind[[3]])
}

# Keep labels present with n >= MIN_N; warn about the rest.
available <- function(pops, counts, what) {
  n <- ifelse(pops %in% names(counts), as.integer(counts[pops]), 0L)
  keep <- pops[n >= MIN_N | pops %in% c("Russia_UstIshim_IUP", "Russia_Malta_UP")]
  gone <- setdiff(pops, keep)
  if (length(gone) > 0)
    cat(sprintf("  WARNING (%s): not in panel or n<%d, skipped: %s\n",
                what, MIN_N, paste(gone, collapse = ", ")))
  keep
}

set2_dirs <- function(merged_prefix, target) {
  out <- file.path(dirname(merged_prefix), SET_NAME)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  list(out = out,
       prefix = file.path(out, paste0("set2_", target)),
       f2 = file.path(out, paste0("f2cache_", target)))
}
