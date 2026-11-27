# ==============================================================================
# set2/pops.R -- population lists shared by every set-2 script (sourced, not run)
#
# Set 1 (steps 04-09) compares PW with candidates against Mbuti/Yoruba/Han only.
# Those outgroups cannot tell one European source from another (34/38 single
# sources fit PW on WGS x HO). Set 2 uses ancient, differentially related
# outgroups and lets the candidates act as outgroups for each other.
#
# Labels verified present in AADR v66.p1 HO and 1240K (2026-10-21). Scripts
# skip a candidate the panel lacks (with a warning) but stop if a distal
# reference or a fixed outgroup is missing.
# ==============================================================================

# Fixed outgroup for f3/f4.
OUTGROUP <- "Mbuti"

# Distal right set: deep lineages that are differentially related to the
# European ancestral components (Lazaridis/Olalde-style).
RIGHT_DISTAL <- c(
  "Mbuti",                          # African root
  "Russia_UstIshim_IUP",            # 45 ka, basal non-African Eurasian
  "Russia_Kostenki_UP",             # 37 ka, early European
  "Russia_Malta_UP",                # MA1, Ancient North Eurasian
  "Han",                            # East Asian
  "Papuan",                         # Oceanian / Denisovan-admixed
  "Karitiana",                      # Native American (ANE-related)
  "Israel_Natufian",                # Levant HG (Basal Eurasian)
  "Iran_GanjDareh_N",               # Iran Neolithic
  "Georgia_KotiasKlde_Mesolithic",  # Caucasus HG (CHG)
  "Morocco_Iberomaurusian"          # Taforalt, North African UP
)

# Distal sources.
WHG    <- "Luxembourg_Loschbour_Mesolithic"
EHG    <- "Russia_Karelia_Mesolithic_HG"
EEF    <- "Turkey_N"
STEPPE <- "Russia_Samara_EBA_Yamnaya"

DISTAL_MODELS <- list(
  "WHG+EHG+EEF"    = c(WHG, EHG, EEF),
  "WHG+Steppe+EEF" = c(WHG, STEPPE, EEF)
)

# Proximal candidate pool: rotated (each model's unused candidates join the
# right set). Ancient only -- a modern population in the pool would simply
# win as the nearest relative and say nothing about the past.
CANDIDATES <- c(
  "Ukraine_IA_Lusatian",              # set-1 leader
  "Poland_IA",                        # label says IA, but 128/130 are 900-1200 CE Piast-era
  "Poland_BA",                        # n=25, mostly EARLY Bronze Age 2300-1700 BCE (Chyleński
                                      # 2023); only 3 are 1600-1100 BCE -- set 3 splits them
  "Poland_BA_Trzciniec",
  "Poland_CordedWare",
  "Poland_IA_Wielbark",
  "Poland_IA_Przeworsk",              # n=2
  "Poland_Roman_Wielbark",
  "Poland_EarlySlav",
  "Poland_EarlyMedieval_Slav",
  "Czechia_EarlyMedieval_EarlySlav",
  "Lithuania_LBA",
  "Sweden_IA",
  "Slovakia_MigrationPeriod",
  "Croatia_EarlyMedieval_Avar",
  "Russia_Ivanovo_EarlyMedieval",
  "Hungary_Conqueror_Commoner"
)

# Present-day groups: compared in f3 and distal tables, never rotated.
MODERNS <- c("Belarusian", "Ukrainian_North", "Ukrainian", "Lithuanian",
             "Russian", "Czech", "Hungarian", "Romanian")

SET_NAME <- "set2"   # output subdirectory (set4/pops.R overrides)

MIN_N <- 2L   # Przeworsk has n=2; everything else in the lists has more

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
