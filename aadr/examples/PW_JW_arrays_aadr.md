# PW and JW consumer arrays vs AADR v66.p1 HO — set 1 (steps 04–09)

Rerun 2026-10-22 on **correctly decoded** AADR genotypes (see "TGENO" in
`../../CLAUDE.md`). Replaces the June 2026 reports `PW_FF_slavic_aadr.md` and
`JW_MH_slavic_aadr.md`, whose AADR genotypes were mis-decoded; their conclusions
are void and the files were removed (they remain in git history only).

| Target | Data | Panel SNPs typed | SNPs in the f2 caches |
|---|---|---|---|
| PW | FamilyTreeDNA Family Finder array | ~57K of 584K HO | 57K (04/05), 48K (06), 20K (07), 45K (08), 20K (09) |
| JW | MyHeritage array | ~57K of 584K HO | 57K (04/05), 48K (06), 20K (07), 46K (08), 20K (09) |

Set 1 uses only modern outgroups (Mbuti, Yoruba, Han; step 05 rotates six sets).
These outgroups see Europe as almost one point, so set 1 can reject sources with
East-Asian, Siberian, Near-Eastern or early-farmer excess, but it cannot rank one
European source above another. Read every `OK` below as "not rejected", never as
"best". For discriminating results see `PW_WGS_aadr.md` (set 2–4, WGS).

## Steps 04–05: Germany / Corded Ware component

Model A = Poland_EarlySlav + Germany_Esperstedt_CordedWare.

| Right set | PW p | PW Germany weight (z) | JW p | JW Germany weight (z) |
|---|---|---|---|---|
| Mbuti+Yoruba+Han | 0.80 | −0.79 (−0.12) | 0.10 | 2.71 (0.97) |
| Mbuti+Yoruba+Papuan | 0.82 | −0.67 (−0.24) | 0.11 | 4.92 (0.39) |
| Mbuti+Han+Biaka | 0.43 | −2.45 (−0.38) | 0.18 | 2.19 (1.03) |
| Mbuti+Yoruba+Nganasan | 0.45 | 1.08 (0.54) | 0.09 | 1.88 (1.19) |
| Mbuti+Yoruba | 0.57 | −1.19 (−0.36) | 0.05 | 3.05 (0.82) |
| Mbuti+Han | 0.96 | 0.09 (0.04) | 0.11 | 2.27 (0.89) |

The Germany/CWC weight is never significant (|z| ≤ 1.2) and its SE is larger
than the weight: no evidence for a separate Germany_IA/CWC component. This
June conclusion survives on correct data, and on WGS too (6/6 HO, 5/5 1240K).

## Step 06: Polish reference populations

| Model | PW p | JW p |
|---|---|---|
| Poland_IA_Wielbark | 0.889 | 0.055 |
| Poland_Roman_Wielbark | 0.830 | 0.100 |
| Poland_EarlySlav | 0.814 | 0.133 |
| Poland_EarlyMedieval_Slav | 0.619 | 0.413 |
| Poland_CordedWare | 0.499 | 0.089 |
| Poland_IA (medieval 900–1200 CE despite the label) | 0.482 | 0.276 |
| Poland_BA_Trzciniec | 0.461 | 0.293 |
| Russia_Samara_EBA_Yamnaya | 0.136 | 0.026 (rejected) |
| Poland_GlobularAmphora | 0.006 (rejected) | <1e-4 (rejected) |
| Slovakia_N_LBK | 0.002 (rejected) | <1e-4 (rejected) |
| Yamnaya + LBK, Yamnaya + GAC, Yamnaya + LBK + GAC | rejected | rejected |

Every Polish source from Corded Ware to the Piast era fits both targets; even
Yamnaya alone is not rejected for PW. Only pure early farmers are excluded. All
2-source Polish models are `OK-unphysical` (one weight negative): the sources
are collinear against these outgroups.

## Step 07: 38 single sources outside Poland_*

Rejected for both targets (p < 0.05): Jew_Ashkenazi, Germany_Medieval_Jewish,
Austria_Avar, Hungary_EarlyAvar-oLowEastAsia, Russia_Tatarstan_EarlyMedieval_EarlyBulgar.
The other 33 are not rejected for either target (p 0.05–0.998), including
Avars from Croatia/Hungary, Huns, Sarmatians, Hungarian Conquerors, Vikings,
Wielbark, every Baltic and Slavic group, and modern Northern/Eastern Europeans.

| Top by p | PW p | JW p |
|---|---|---|
| Lithuania_LBA | 0.961 | 0.969 |
| Denmark_IA | 0.946 | 0.885 |
| Estonian | 0.892 | 0.646 |
| Lithuanian | 0.859 | 0.998 |
| Poland_Roman_Wielbark | 0.857 | 0.678 |
| Belarusian | 0.684 | 0.895 |
| Ukraine_IA_Lusatian | 0.596 | 0.329 |

The June "Ukraine_IA_Lusatian p ≈ 0.9–1.0, best single source" result is void.
On correct data it is one of 33 non-rejected sources; on WGS with ancient
outgroups (set 2) it is rejected.

## Step 08: rotating 2-source scan (bases Ukraine_IA_Lusatian, Poland_EarlySlav; 35 second sources each)

| Verdict | PW | JW |
|---|---|---|
| SIGNIFICANT-2src | 0 | 0 |
| OK-2src-nonsignificant / OK-second-dominant / OK-base-sufficient | 16 | 8 |
| UNIDENTIFIABLE | 11 | 19 |
| REJECTED | 43 | 43 |

No hidden second component is detectable. The rejected second sources are the
ones with Siberian, Central-Asian, Caucasian, Near-Eastern or Ashkenazi ancestry
(Altaian, Buryatia Xiongnu, Kazakh, Nganasan, Mansi, Armenian, Georgian, Iranian,
Druze, Jew_Ashkenazi, Hungary_Conqueror_Commoner ...). With a European second
source, models fit but the weight is never significant (max |z_second| < 2).

## Step 09: Przeworsk vs Wielbark

| Test | PW | JW |
|---|---|---|
| 1-source Poland_IA_Przeworsk (n=2) | p 0.774 | p 0.359 |
| 1-source Poland_IA_Wielbark | p 0.663 | p 0.261 |
| 1-source Ukraine_IA_Lusatian | p 0.190 | p 0.866 |
| B1 f4(Mbuti, Sweden_IA; Przeworsk, Wielbark) | z 1.16 (n.s.) | z 1.06 (n.s.) |
| B2 f4(Mbuti, T; Wielbark, Przeworsk) | z −1.10 (n.s.) | z −0.36 (n.s.) |
| B4 same with Han | z −0.76 (n.s.) | z −0.03 (n.s.) |

Przeworsk and Wielbark cannot be told apart on ~20K SNPs, and neither target
leans to one of them. On WGS with ancient outgroups (set 2) both are rejected as
sole sources for PW; in set 4 Przeworsk returns only as part of a 2-source model
(Przeworsk ~0.6 + Latvia_LBA_FlatGraves ~0.4).

## PW vs JW

Rankings differ between the two targets, but every difference lies inside the
noise of a test that cannot rank European sources. The June claim "JW drifts
SOUTH" is void. The direct test on shared SNPs (`../pair_compare.R`, in
`PW_WGS_aadr.md`) finds PW and JW indistinguishable (qpWave p 0.48–0.86; 0 of
63 f4 statistics with |z| ≥ 3).

## What a consumer array can and cannot answer here

- Can: exclude gross components (East-Asian/Siberian steppe, Near-Eastern,
  Ashkenazi, pure early farmers).
- Cannot: say which European population, period or region a target descends
  from. Neither more candidates nor more 2-source models fix that. It needs ancient
  outgroups (set 2) and more SNPs (WGS, `01b_prepare_wgs.sh`). Even set 2 on
  the JW array (~10K SNPs in the joint cache) fits almost everything.
