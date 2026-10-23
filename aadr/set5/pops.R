# ==============================================================================
# set5/pops.R -- set 4 (pre-medieval candidates only) + SOUTH-EASTERN candidates
#
# Set 4 had only Baltic groups as "eastern" sources. Set 5 adds the Ukrainian
# forest-steppe Iron Age and the Chernyakhiv-era (Roman period) groups of
# Saag et al. 2025 (Sci Adv), still with no medieval or modern population
# anywhere. Used by:
#   - the set-2 builder: POPS_FILE=aadr/set5/pops.R Rscript aadr/set2/10_build_f2.R ...
#   - aadr/set5/30_trzciniec_tests.R (Trzciniec admixture tests, pooling of
#     the small forest-steppe labels, Gustorzyn split)
# Do NOT run 13_rotating.R on this list: MIN_N is 1 here, so singletons are
# extracted for pooling and would rotate as single-individual sources.
# ==============================================================================

SET_NAME <- "set5"
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
  "Germany_Magdalenenberg_IA_Hallstatt",
  # --- set5: south-east (forest-steppe / Chernyakhiv era), Saag et al. 2025 ---
  "Ukraine_LBA_Vysotska",                             # Petrykiv, ~1150 BCE, n=2, LOW coverage
  "Ukraine_EIA1_ForestSteppe_Scythian",               # pooled in 30_ as Ukraine_EIA_ForestSteppe
  "Ukraine_EIA1_ForestSteppe_Scythian-1",
  "Ukraine_EIA2_ForestSteppe_Scythian",
  "Ukraine_EIA2_BuhDniester_ForestSteppe_Scythian",
  "Ukraine_Kyiv_EIA1_ForestSteppe_Scythian-1",
  "Ukraine_MIA_Scythian",                             # Bilsk hillfort, forest-steppe
  "Ukraine_LIA",                                      # Petrykiv (Ternopil), ~240 BCE, n=1
  "Ukraine_Classical",                                # Shushaky, 235-350 CE, n=10
  "Ukraine_Shushaku_Classical_Chernyakhov",           # n=1, pooled with Ukraine_Classical? kept apart
  "Ukraine_Roman_BuhDniester",                        # n=1
  "Moldova_Scythian"
)

MODERNS <- character(0)

# set5 extracts singletons too: 30_trzciniec_tests.R pools them into
# multi-site groups (forest-steppe EIA) and never rotates them alone.
MIN_N <- 1L

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
