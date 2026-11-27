# ==============================================================================
# set3/regions.R -- medieval Polish burials regrouped by historical region
# (sourced by set-3 scripts, not run)
#
# AADR labels follow publications, not geography: Poland_IA (n=130) is in fact
# 128 Piast-era burials 900-1200 CE from 19 sites, Poland_EarlySlav is Ląd and
# Łęgowo, Poland_EarlyMedieval_Slav is Santok. Set 3 reassigns each medieval
# individual (Poland, 400-1450 CE in the .anno, "-o" outliers excluded) to a
# region by SITE NAME -- not by coordinates, which are wrong for at least
# "Poznań Sołacz" (given as 54.49 N, i.e. Słupsk).
#
# Region -> sites (from v66.p1 .anno, checked 2026-10-22):
# ==============================================================================

REGIONS <- list(
  Reg_GreaterPoland_Gniezno = c("Ostrów Lednicki", "Dziekanowice-22", "Gniezno", "Gołuń",
                                "Łęgowo", "Zielonka", "Oblaczkowo"),
  Reg_GreaterPoland_Poznan  = c("Poznań", "Poznań-Śródka", "Poznań Sołacz", "Sowinki"),
  Reg_GreaterPoland_East    = c("Ląd", "Rumin"),
  Reg_Kuyavia               = c("Markowice", "Balczewo"),
  Reg_Silesia               = c("Niemcza", "Milicz", "Groszowice"),
  Reg_Masovia               = c("Płońsk", "Czersk"),
  Reg_LesserPoland          = c("Końskie", "Krakow", "Sandomierz"),
  Reg_Lubusz_Santok         = c("Santok", "Cedynia"),
  Reg_Bodzia_elite          = c("Bodzia")      # elite cemetery with Scandinavian/eastern links
)
# Poland_BA (n=25) split by date: AADR has almost no Lusatian-culture samples
# from Poland. 22 are Early Bronze Age (2300-1700 BCE); 3 at Świniary Stare /
# context-dated 1600-1100 BCE are the closest thing to the Lusatian window.
PERIOD_SPLITS <- list(
  list(group = "Poland_BA", label = "Per_Poland_EBA",          min_bp = 3501, max_bp = 99999),
  list(group = "Poland_BA", label = "Per_Poland_BA_late_n3",   min_bp = 0,    max_bp = 3500)
)

MIN_REGION_N <- 5L      # smaller regions are reported but not used as qpAdm sources

# Non-regional references (labels as in AADR).
REF_PERIODS <- c(
  "Poland_BA_Trzciniec",             # 1500-1300 BCE
  "Ukraine_IA_Lusatian",             # Volhynia, 1000-700 BC
  "Poland_IA_Wielbark",
  "Poland_IA_Przeworsk",
  "Czechia_EarlyMedieval_EarlySlav",
  "Lithuania_LBA",
  "Sweden_IA"
)

# Deep outgroups: set-2 right set minus Russia_Malta_UP (n=1, ~1x; it alone
# cuts most SNPs under maxmiss=0).
RIGHT3 <- c("Mbuti", "Russia_UstIshim_IUP", "Russia_Kostenki_UP", "Han", "Papuan",
            "Karitiana", "Israel_Natufian", "Iran_GanjDareh_N",
            "Georgia_KotiasKlde_Mesolithic", "Morocco_Iberomaurusian")
DISTAL3 <- c("Luxembourg_Loschbour_Mesolithic", "Russia_Karelia_Mesolithic_HG",
             "Turkey_N", "Russia_Samara_EBA_Yamnaya")

# ---- helpers -----------------------------------------------------------------

anno_path <- function(aadr_prefix) {
  # v66.p1_HO.aadr.patch.PUB -> v66.p1_HO.aadr.PUB.anno (symlinked next to the geno)
  p <- paste0(sub("\\.patch\\.PUB$", ".PUB", aadr_prefix), ".anno")
  if (!file.exists(p)) stop("annotation file not found: ", p)
  p
}

# Individuals (IID, region, site, date, HO SNPs) for every region.
region_members <- function(aadr_prefix) {
  a <- read.delim(anno_path(aadr_prefix), stringsAsFactors = FALSE, quote = "",
                  check.names = FALSE, comment.char = "")
  id <- a[[1]]; grp <- a[[15]]; site <- gsub('"', "", a[[16]]); ctry <- a[[17]]
  bp <- suppressWarnings(as.numeric(a[[11]]))
  site <- sub(" \\(.*", "", site)
  keep <- ctry == "Poland" & !is.na(bp) & bp >= 500 & bp <= 1550 &   # 400-1450 CE
          !grepl("-o", grp, fixed = TRUE)
  out <- data.frame(IID = id[keep], group = grp[keep], site = site[keep],
                    date = a[[13]][keep], stringsAsFactors = FALSE)
  out$region <- NA_character_
  for (r in names(REGIONS)) out$region[out$site %in% REGIONS[[r]]] <- r
  un <- unique(out$site[is.na(out$region)])
  if (length(un) > 0) cat(sprintf("  WARNING: medieval Polish sites with no region (skipped): %s\n",
                                  paste(un, collapse = ", ")))
  out[!is.na(out$region), ]
}

# Individuals of PERIOD_SPLITS with their new labels.
period_members <- function(aadr_prefix) {
  a <- read.delim(anno_path(aadr_prefix), stringsAsFactors = FALSE, quote = "",
                  check.names = FALSE, comment.char = "")
  bp <- suppressWarnings(as.numeric(a[[11]]))
  do.call(rbind, lapply(PERIOD_SPLITS, function(sp) {
    k <- a[[15]] == sp$group & !is.na(bp) & bp >= sp$min_bp & bp <= sp$max_bp
    data.frame(IID = a[[1]][k], label = rep(sp$label, sum(k)), stringsAsFactors = FALSE)
  }))
}

set3_dirs <- function(merged_prefix, target) {
  out <- file.path(dirname(merged_prefix), if (Sys.getenv("ALLSNPS", "0") == "1") "set3_allsnps" else "set3")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  list(out = out,
       prefix = file.path(out, paste0("set3_", target)),
       f2 = file.path(out, paste0("f2cache_", target)))
}
