# PW — 30x WGS vs AADR v66.p1 (2026-10-21)

First results on **correctly decoded** AADR genotypes. All June 2026 AADR results
were computed on TGENO files that admixtools mis-decoded and are superseded — see
"TGENO" in `../../CLAUDE.md`; the consumer-array reruns are in `PW_JW_arrays_aadr.md`.

Target: consumer 30x WGS (CRAM, GRCh37), genotyped at panel sites by
`01b_prepare_wgs.sh` (578,836 HO / 1,147,026 1240K autosomal SNPs). Run dirs:
`./PW_WGS_HO/`, `./PW_WGS_1240K/` (not published).

> **Label warning:** `Poland_IA` is not Iron Age. 128/130 individuals are
> 900-1200 CE Piast-state burials (Stolarek et al. 2023). See `../../CLAUDE.md`.

## Set 2 (ancient outgroups) — the results that discriminate

SNPs after `maxmiss=0` across 40 populations: HO 46,927; 1240K 82,172.

### Rotating qpAdm, 1 source (right = 11 deep outgroups + WHG/EHG/EEF/Steppe + other candidates)

| Source | p HO | p 1240K | Verdict |
|---|---|---|---|
| Poland_EarlySlav | 0.605 | 0.146 | fits |
| Lithuania_LBA | 0.475 | 0.123 | fits |
| Poland_IA (despite the name: Piast-era 900-1200 CE) | 0.374 | 0.535 | fits |
| Poland_EarlyMedieval_Slav | 0.140 | 0.441 | fits |
| Ukraine_IA_Lusatian | 0.030 | 0.0065 | rejected |
| Czechia_EarlyMedieval_EarlySlav | 0.028 | 0.030 | rejected |
| Poland_BA_Trzciniec | 0.016 | 6e-4 | rejected |
| Russia_Ivanovo_EarlyMedieval | 2e-3 | 2e-4 | rejected |
| Poland_Roman_Wielbark | 8e-5 | 1e-6 | rejected |
| Sweden_IA | 6e-5 | 5e-5 | rejected |
| Poland_IA_Wielbark | 1e-5 | 1e-7 | rejected |
| Poland_IA_Przeworsk (n=2) | 5e-6 | 2e-9 | rejected |
| Slovakia_MigrationPeriod, Hungary_Conqueror_Commoner, Poland_CordedWare, Croatia_EarlyMedieval_Avar | < 1e-9 | < 1e-14 | rejected |

The same four sources fit on both panels and the same twelve are rejected.

### Rotating qpAdm, 2 sources

No 2-source model is needed: every identifiable 2-source fit with p above 0.13
contains one of the four 1-source fits. The "both parts rejected alone" fits are
borderline (p 0.05–0.13; the 0.127 one on 1240K is Czechia_EarlyMedieval_EarlySlav
with a −0.11 ± 0.12 second weight, i.e. Czechia alone). One pattern repeats on
both panels:
Ukraine_IA_Lusatian ~0.75 + Croatia_EarlyMedieval_Avar or Slovakia_MigrationPeriod
~0.25 — i.e. the Volhynian Lusatians alone are slightly too far from PW, and a
Central/South-European component closes the gap. Weak evidence, not a finding.
Models whose weight SE exceeds 0.25 are reported as UNIDENTIFIABLE (HO 44, 1240K 31).

### f4 scan

f4(Mbuti, X; S, PW) for S = the best three sources: **0 of 31–37 populations
with |z| ≥ 3** on either panel. PW and Poland_IA / Poland_EarlySlav /
Poland_EarlyMedieval_Slav / Lithuania_LBA are symmetric to everything tested.

### Distal (right = 11 deep outgroups)

| Model | Panel | p | WHG | Steppe | EEF |
|---|---|---|---|---|---|
| WHG+Steppe+EEF | HO | 0.299 | 0.04 ± 0.08 | 0.59 ± 0.08 | 0.37 ± 0.06 |
| WHG+Steppe+EEF | 1240K | 0.894 | 0.05 ± 0.06 | 0.62 ± 0.06 | 0.33 ± 0.05 |

WHG+EHG+EEF is marginal on HO (p 0.06) and rejected on 1240K. SEs of 0.05–0.08
make almost every candidate distally indistinguishable from PW; this table places
PW, it does not rank sources.

### Outgroup f3(Mbuti; PW, X)

Top 15 within 0.8 SE of each other (Ukraine_IA_Lusatian nominally first) — a tie.

## Set 1 (Mbuti/Yoruba/Han outgroups) on correct data

- Step 04/05: Germany/CWC weight not significant, 6/6 outgroup sets (HO), 5/5 (1240K).
- Step 07 HO: 12/38 rejected — every group with East-Asian/Siberian ancestry
  (Avars, Sarmatians, Early Bulgars, Hungarian Conquerors, Russians, Finns,
  Estonians) plus Jew_Ashkenazi and Germany_Medieval_Jewish. Among European
  sources it cannot choose (Denmark_IA p 0.96 tops a list set 2 rejects).
- Step 09: f4(Mbuti, PW; Wielbark, Przeworsk) z = −0.61 (n.s.); set 2 rejects
  both as sole sources anyway.

## What changed against the (invalid) array-era conclusions

- "Ukraine_IA_Lusatian best single source" — **not supported**: rejected on both
  panels once candidates compete. What fits is medieval Poland (Poland_IA = Piast-era
  cemeteries 900-1200 CE, Poland_EarlySlav, Poland_EarlyMedieval_Slav) and
  Lithuania_LBA -- i.e. PW is continuous with the early medieval Polish population.
- Wielbark / Przeworsk as sole sources — **rejected**.
- No Avar / steppe-nomad / Ashkenazi ancestry — now supported by real data (set 1
  step 07 rejections; no set-2 model needs them).

Caveats: Poland_IA_Przeworsk n=2; set-2 SNP counts are limited by low-coverage
ancient references (Malta, Kostenki); p-values across sources are not a ranking.

## Set 3 — medieval Polish regions (strict mode; added 2026-10-22)

Regions built by site from the AADR .anno (see `aadr/set3/`). Rotating qpAdm,
1 source, PW (HO / 1240K):

| Region / group | HO p | 1240K p | |
|---|---|---|---|
| Greater Poland East (Ląd, Rumin) | 0.353 | 0.174 | fits |
| Lubusz (Santok, Cedynia) | 0.352 | 0.107 | fits |
| Masovia (Płońsk, Czersk) | 0.331 | 0.179 | fits |
| Czechia_EarlyMedieval_EarlySlav | 0.325 | 0.141 | fits |
| Greater Poland — Poznań | 0.141 | 0.352 | fits |
| Greater Poland — Gniezno | 0.084 | 0.214 | fits |
| Silesia (Niemcza, Milicz, Groszowice) | 0.003 | 0.003 | rejected |
| Kuyavia (Markowice, Balczewo) | 0.003 | 5e-5 | rejected |
| Bodzia elite | 0.006 | 0.028 | rejected |

Regions are barely distinguishable from each other against deep outgroups
(qpWave: 0/28 pairs on HO, 6/28 on 1240K, Kuyavia in 5 of those 6). Not tested
(too sparse): Lesser Poland, Przeworsk, Lithuania_LBA, the 3 late-BA individuals.
`ALLSNPS=1` mode was tried and is unreliable (see `../set3/README.md`).

## Set 4 — pre-medieval sources only (added 2026-10-22)

22 LBA/IA/Roman groups around Poland, no medieval or modern population anywhere.
SNPs: HO 66,844; 1240K 115,898 (JW array: 6,991 — too few to discriminate).

- **Every single pre-medieval source is rejected for PW on both panels**
  (best: Trzciniec p 0.02 / 3e-4; Poland_BA 0.003 / 1e-4; Ukraine_IA_Lusatian 3e-4 / 2e-6).
- 2-source models that FIT with identifiable weights on BOTH panels:

| Model | HO weights, p | 1240K weights, p |
|---|---|---|
| Poland_IA_Przeworsk + Latvia_LBA_FlatGraves | 0.61/0.39, 0.287 | 0.56/0.44, 0.062 |
| Lithuania_LBA + Sweden_IA | 0.37/0.63, 0.163 | 0.32/0.68, 0.211 |
| Estonia_LBA_StoneCistGraves + Sweden_IA | 0.34/0.66, 0.156 | 0.27/0.73, 0.138 |
| Ukraine_IA_Lusatian + Sweden_IA | 0.39/0.61, 0.134 | 0.36/0.64, 0.153 |
| Poland_BA_Trzciniec + Sweden_IA | 0.61/0.39, 0.132 | 0.50/0.50, 0.085 |
| Latvia_LBA_FlatGraves + Sweden_IA | 0.31/0.69, 0.115 | 0.35/0.65, 0.254 |
| Poland_Roman_Wielbark + Lithuania_LBA | 0.58/0.42, 0.114 | 0.63/0.37, 0.114 |
| Poland_IA_Wielbark + Lithuania_LBA | 0.56/0.44, 0.092 | 0.63/0.37, 0.078 |
| Poland_Roman_Wielbark + Latvia_LBA_FlatGraves | 0.65/0.35, 0.084 | 0.63/0.37, 0.131 |
| Estonia_EIA_TarandGraves + Germany_Magdalenenberg_IA_Hallstatt | 0.63/0.37, 0.073 | 0.66/0.34, 0.065 |
| Estonia_EIA_TarandGraves + Czechia_LBA | 0.57/0.43, 0.059 | 0.60/0.40, 0.056 |
| Estonia_EIA_TarandGraves + Czechia_IA_LaTene | 0.57/0.43, 0.056 | 0.59/0.41, 0.061 |

Common structure: **an East-Baltic LBA/EIA group (~30-45 %) + a Central or
Northern European Iron Age group (~55-70 %)**. Przeworsk + Latvia_LBA is the
best-supported pair and the only one that also fits JW with the same weights
(0.62/0.38, p 0.72; JW otherwise uninformative). Caveats: ~230 two-source models
were tested, so a few fits at p 0.05-0.15 are expected by chance; Przeworsk is
n=2; this says which ANCESTRY ring reproduces PW, not which people migrated.

## PW vs JW (added 2026-10-22, `../pair_compare.R`)

Shared SNPs: 57,238 (JW's MyHeritage array ∩ HO panel), PW WGS + JW array merged.

| Test | Medieval-regions panel (set 3 groups) | Pre-medieval panel (set 4 groups) |
|---|---|---|
| qpWave PW,JW clade, right = 10 deep outgroups | p 0.864 | p 0.607 |
| qpWave PW,JW clade, right = all other groups | p 0.482 (28 groups) | p 0.652 (37 groups) |
| f4(Mbuti, X; PW, JW), \|z\| >= 3 | 0 of 27 | 0 of 36 |
| f4 sign | 23/27 positive, mean z +0.64 | 35/36 positive, mean z +1.39 |

PW and JW are not distinguishable. The near-uniform positive f4 (every X
"slightly closer to JW") is a target-level offset, most likely WGS-vs-array
genotyping on the same SNPs, not ancestry toward any particular group. Distal
WHG/Steppe/EEF shares overlap within SE and swap order between the two panels.
The earlier "JW drifts SOUTH" finding came from mis-decoded AADR data and is void.

## Own PCA + an external nMonte runner (added 2026-10-22, `../pca/`)

`30_build_pca.py`: 25-PC West-Eurasian PCA from converted AADR HO (1,439 modern
basis individuals, 426,281 SNPs), lsqproject of 17,285 individuals + PW (WGS) +
JW (array). Population distances vs Davidski G25: r = 0.95 (68 modern pops),
0.993 (15 European pops). `32_to_g25_centroids.py` rebuilds the runner's
catalogue ids (Davidski/G25 cluster names of a separate nMonte runner, not
included here) from AADR sites: 118/184 available. NOT available (not in AADR):
Germany_IA (Saxony-Anhalt clusters), Hungary_Mig_Himod, Latvia_Neo_Kivutkalns,
Gródek, Ukrainian medieval sites -- so the canonical G25 models cannot be
reproduced like-for-like; the catalogue also has no Baltic LBA source.

| Stage | PW best model (dist) | JW best model (dist) |
|---|---|---|
| 01D mixed free | Płońsk 69% + Santok 16% + Germany_BA/CW Esperstedt 15% (0.0207) | Ląd 51% + Płońsk 41% + Markowice 4% + Anderten 5% (0.0299) |
| 04D weighted | Płońsk 79% + Germany_CW 12% + Fröjel 6% + Sintashta 3% (0.0205) | Płońsk 74% + Drantum 20% + ... (0.0316) |
| 03B deep only | Gustorzyn (Trzciniec) 68% + Trzciniec 18% + Krusza 8% + Sintashta 6% (0.0271) | – |

Medieval Masovia (Płońsk) dominates both; agrees with qpAdm set 3 (Masovia,
Greater Poland, Lubusz fit). Distances are in this PCA's units (not G25's).
