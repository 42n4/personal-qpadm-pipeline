# PW — 30x WGS vs AADR v66.p1 + Gretzinger et al. 2025 (2026-10-23)

Extends `PW_WGS_aadr.md` (sets 1–4, PW vs JW, own PCA) with the south-eastern
candidates AADR already had but set 4 never used (**set 5**), the Migration- and
Slavic-period genomes of Gretzinger et al. 2025 (**set 6**), and a PCA with those
590 individuals projected into the same space.

Target: consumer 30x WGS (CRAM, GRCh37), genotyped at panel sites by
`01b_prepare_wgs.sh` (578,836 HO / 1,147,026 1240K autosomal SNPs). Every result
below was computed on **both** panels; nothing is quoted that holds on only one.

> **Label warnings.** `Poland_IA` is not Iron Age — it is 900–1200 CE Piast-era
> burials. `MP` = Migration Period, `SP` = Slavic Period: chronological labels, not
> ethnic verdicts (MP medians 275–500 CE, SP medians 775–1250 CE).

## The question

`PW_WGS_aadr.md` left set 4 with: no single pre-medieval source fits PW; the pairs
that fit are an East-Baltic LBA group plus a Central/Northern European Iron Age one;
and PW has ~3× less hunter-gatherer ancestry than Trzciniec at an unchanged steppe
share. So something richer in early-farmer ancestry diluted the Bronze Age Polish
profile between ~1200 BCE and the Middle Ages. Sets 5 and 6 ask **which direction**
that came from, and whether it is visible before or after the Migration Period.

## New reference data

**Gretzinger et al. 2025**, *Nature*, doi:10.1038/s41586-025-09437-6 — genotypes
published as the Poseidon package `2025_Gretzinger_SlavicPeriod` v1.0.0 (PLINK,
snpSet 1240K, CC-BY). All 1,233,013 SNP ids are **identical to the AADR 1240K panel**
and HO is a strict subset of it, so the package merges with no liftover and no id
translation; all three MD5 sums match `POSEIDON.yml`.

| Group | n used | Median date | Sites |
|---|---|---|---|
| Poland_MP | 14 | 334 CE | Gródek, Koziarnia, Żarska, Boczne |
| Poland_SP | 15 | 1074 CE | Gródek |
| NWUkraine_SP | 25 | 1250 CE | Korolivka, Pidhirtsi (Volhynia/Galicia) |
| EGermany_MP | 117 | 500 CE | Saxony-Anhalt |
| EGermany_SP | 210 | 1100 CE | Saxony-Anhalt |
| NWBalkan_MP / SP | 21 / 70 | 275 / 775 CE | Croatia |
| Latvia_MA | 4 | 1027 CE | Latvia |

Individuals below 50,000 SNPs were dropped (114 of 590). Joint cache after
`maxmiss=0`: HO 55,180 / 1240K 97,618 SNPs.

## Set 5 — south-eastern candidates (`set5/`)

Set 4's only "eastern" sources were Baltic. Set 5 adds the Ukrainian forest-steppe
Iron Age and Chernyakhiv-era groups (Saag et al. 2025, in AADR as
`SaagThomasSciAdv2025`) plus `Moldova_Scythian` (Krzewińska et al. 2018), still with
no medieval or modern population anywhere.

### Nothing is closer to PW than Trzciniec

f4(Mbuti, X; Trzciniec, PW) — a positive z would mean X is closer to PW:

| X | z HO | Reading |
|---|---|---|
| Ukraine_Chernyakhiv_era | +0.02 | equidistant |
| Turkey_N (EEF) | +0.02 | equidistant |
| Moldova_Scythian | −0.48 | — |
| Latvia_LBA_FlatGraves | −2.24 | nearer Trzciniec |
| Yamnaya / Poland_BA / Loschbour | −2.9 … −3.0 | nearer Trzciniec |

**No candidate gives a significantly positive z on either panel.** The south-eastern
groups sit equidistant; the Baltic and older Polish groups sit on Trzciniec's side.

### Kuyavian Trzciniec is not a closer relative

Trzciniec in AADR is >80 % Lesser Poland; Gustorzyn (Kuyavia, n=6) is the only
northern series. Splitting it out:

| Statistic | z HO | z 1240K |
|---|---|---|
| f4(Mbuti, PW; Trz_rest, Trz_Gustorzyn) | −0.60 | −0.46 |
| f4(Mbuti, Trz_Gustorzyn; Trz_rest, PW) | **−3.11** | **−3.86** |

PW does not lean toward the Kuyavian series, and Gustorzyn is closer to the rest of
Trzciniec than to PW. The regional split changes nothing.

### Direction remains undecided

f4(PW, Trzciniec; Central-European, south-eastern) over 9 pairs × 2 panels:
**all |z| < 1.8**. The weak, non-significant lean is toward `Ukraine_Chernyakhiv_era`.

A methodological catch worth recording: in set 4 the Central-European two-source
models were *rejected* on 1240K (p 0.012–0.032); in set 5, with the same target and
the same panels, they *pass* (p 0.088–0.25). The only change is which groups sit in
the rotating right set. **A p in the 0.05–0.2 band is a property of the design, not
of the past** — see "What this does not show".

## Set 6 — Migration- and Slavic-period genomes (`set6/`)

`set6/40_add_external.R` merges an external PLINK package into a set-N fileset
(coverage filter, allele-conflict retry, target ∩ panel re-filter), then the set-5
tests run unchanged with `SET=set6`. The merge produced **zero allele conflicts**.

### Every single source is still rejected

| Source | p HO | p 1240K |
|---|---|---|
| Poland_BA_Trzciniec | 0.046 | 0.010 |
| Poland_SP | 1.4e-3 | 1.3e-4 |
| NWUkraine_SP | 3.8e-6 | 3.4e-11 |
| **Poland_MP** | **6.2e-7** | **1.9e-9** |
| EGermany_MP | 8.1e-11 | 3.1e-12 |

**Poland_MP is the most strongly rejected Polish source ever tested for PW.** The
Migration-Period population of south-eastern Poland is not PW's ancestor on its own.
Poland_SP — the Slavic-period layer at the same site — is rejected 3–5 orders of
magnitude less hard.

`Latvia_MA` appears to fit alone (p 0.74 / 0.46). It does not: with n=4 it also fits
with **any** second source at ~0.9 weight and p 0.42–0.79, including Scythians,
Corded Ware and Ryazan. That is the weak-anchor signature behind `aadr_eu`'s
`MIN_BASE_N >= 15` rule. Excluded from every conclusion here.

### Gródek has continuity that excludes PW

| Statistic | z HO | z 1240K |
|---|---|---|
| f4(Mbuti, **PW**; Poland_MP, Poland_SP) | −0.62 | −0.14 |
| f4(Mbuti, **Poland_MP**; Poland_SP, **PW**) | **−3.72** | **−4.45** |
| f4(Mbuti, EGermany_MP; Poland_MP, Poland_SP) | −6.46 | −6.89 |
| f4(Mbuti, Trzciniec; Poland_MP, Poland_SP) | −2.76 | −2.42 |

Read together:

1. PW is **equidistant** from the two phases at Gródek — it does not prefer the
   Slavic-period layer.
2. The two phases **share drift that PW lacks**, significant on both panels: there is
   local continuity across the Migration→Slavic boundary at that site, and PW sits
   outside that line.
3. Poland_MP is strongly Germanic-shifted, consistent with its I1/I2/J2b/E1b
   Y-haplogroups and Gródek's position on the eastern edge of the Wielbark range.
   Poland_SP carries R1a-Z280/CTS1211 instead.
4. Trzciniec is nearer the Migration-Period phase than the Slavic-Period one: the
   earlier layer kept more of the Bronze Age substrate.

### Trzciniec + X: seven second sources still fit

| Second source | HO weights, p | 1240K weights, p |
|---|---|---|
| **Poland_SP** | 0.59 / 0.41, p 0.32 | 0.58 / 0.42, p 0.32 |
| NWUkraine_SP | 0.64 / 0.36, p 0.32 | 0.68 / 0.32, p 0.16 |
| Poland_MP | 0.73 / 0.27, p 0.16 | 0.70 / 0.30, p 0.22 |
| EGermany_SP | 0.69 / 0.31, p 0.19 | 0.71 / 0.29, p 0.14 |
| NWBalkan_SP | 0.73 / 0.27, p 0.17 | 0.70 / 0.30, p 0.22 |
| Moldova_Scythian | 0.75 / 0.25, p 0.23 | 0.80 / 0.20, p 0.059 |
| EGermany_MP | 0.78 / 0.22, p 0.090 | 0.70 / 0.30, p 0.090 |

Trzciniec + Poland_SP is the most consistent (identical weights on both panels), but
real Migration- and Slavic-period samples **did not narrow the direction** that the
Scythian proxies could not. No three-source model meets the strict criterion (all
weights > 0, every SE < 0.1).

## PCA with Gretzinger projected (`pca/`, `--extra-plink`)

`30_build_pca.py` gained `--extra-plink PREFIX`: every individual of a PLINK set is
projected through the existing lsqproject target path, FID as the label. Each SNP's
allele pair is re-checked against the panel and the genotype flipped (`2-g`) when
PLINK's A1 differs — on this package **0 flips, 0 unmatched** across 584,131 HO SNPs.
Basis: 1,439 modern West-Eurasians in 147 populations, 426,281 SNPs.

### Distances to PW (25 PCs, group centroids)

| Group | n | to PW | to JW |
|---|---|---|---|
| EGermany_SP | 260 | **0.0245** | 0.0289 |
| Poland_SP | 15 | 0.0287 | 0.0296 |
| NWUkraine_SP | 25 | 0.0323 | 0.0332 |
| Latvia_MA | 5 | 0.0341 | 0.0462 |
| **Poland_MP** | 16 | **0.0456** | 0.0458 |
| NWBalkan_SP | 80 | 0.0457 | 0.0451 |
| EGermany_MP | 154 | 0.0625 | 0.0638 |
| NWBalkan_MP | 35 | 0.1114 | 0.1063 |

**Every Slavic-period group is closer to PW than every Migration-period group.**

### Do not rank the top three from that table

`EGermany_SP` heads the list partly because it has **n=260**. A centroid estimated
from 260 individuals has standard error ~0.0027; `Poland_SP` (n=15) has ~0.0065, so
its centroid is more than twice as noisy and sits further away for free. Using the
median distance to **individual members**, which sample size does not affect, the
order reverses and the three are a tie:

| Group | n | to centroid | median to members | centroid SE |
|---|---|---|---|---|
| Poland_SP | 15 | 0.0287 | **0.0366** | 0.0065 |
| NWUkraine_SP | 25 | 0.0323 | 0.0374 | 0.0063 |
| EGermany_SP | 260 | **0.0245** | 0.0378 | 0.0027 |

f4 settles it, being immune to both effects — f4(Mbuti, PW; A, B), z < 0 = PW nearer A:

| Comparison | z HO | z 1240K | Verdict |
|---|---|---|---|
| EGermany_SP vs Poland_SP | −1.95 | −1.34 | **n.s.** — indistinguishable |
| **EGermany_SP vs NWUkraine_SP** | **−3.24** | **−3.94** | **PW nearer EGermany_SP** |
| Poland_SP vs NWUkraine_SP | −0.55 | −1.38 | n.s. |

The directional reading does **not** survive a closer look, for two reasons.

**NWUkraine_SP is not one group.** It is two sites that differ in both place and date:

| Site | n | Dates | Coordinates |
|---|---|---|---|
| Pidhirtsi | 12 | 1050–1150 | 50.78 N, 23.78 E — **11 km from Gródek** |
| Korolivka | 13 | 1250 | 48.74 N, 25.99 E — 277 km south-east |

Splitting it and repeating f4(Mbuti, PW; A, B):

| Comparison | z HO | z 1240K |
|---|---|---|
| EGermany_SP vs **Korolivka** (1250, genuinely eastern) | −1.70 | −1.93 — **n.s.** |
| EGermany_SP vs **Pidhirtsi** (1100, beside Gródek) | −2.79 | **−4.04** |
| Pidhirtsi vs Korolivka | +1.15 | +2.08 — n.s. |

The combined signal was driven by **Pidhirtsi — the site that is not eastern at all**.
Against Korolivka, the only genuinely eastern and later series, PW is **not**
significantly nearer EGermany_SP on either panel. So this cannot be an east–west axis.

**The contrast is also confounded with date.** NWUkraine_SP's median is 1250 CE against
1074 (Poland_SP) and 1100 (EGermany_SP) — 150–175 years, on post-Mongol-invasion
Volhynia. "West versus east" is inseparable here from "11th century versus 13th".

### Why qpWave does not settle this question

A first run (rank 0, right = 10 deep outgroups; the 45-population cache, 55,180 HO /
97,618 1240K SNPs) appeared to separate PW from EGermany_SP (p 0.195 / 0.039) while
keeping it a clade with Poland_SP (0.404 / 0.276). That reading is withdrawn.

**A qpWave p never means "clade" — only "not distinguishable at this power."** Cutting
the cache to ~12 populations leaves 314,000–659,000 SNPs, because `maxmiss=0` keeps far
more sites when fewer groups compete. At that power PW separates from **every** medieval
group, Gródek included (EGermany blocks 14/14 and 7/7 reject; PW + Gródek area 0.0127 /
0.0009). At 55K it separated from almost none. That is the test behaving correctly, not
failing — but it means **a qpWave p is incomparable between runs with different SNP
counts**. Report n and the SNP count on the same row, and compare only within one cache.

**At high power, most of what separates PW from these groups is technology, not
descent.** PW is diploid 30x WGS; the references are capture data, largely
pseudo-haploid, with their own error profile, reference bias and deamination. Those
differences are *systematic*, so they grow into detectable "structure" as SNPs
accumulate. PW is also separated from all of them by 800–1000 years of its own drift.
The same artefact is visible elsewhere in this project: in `PW_WGS_aadr.md` the PW-vs-JW
f4 scan gave a near-uniform positive offset, attributed to WGS-vs-array genotyping on
shared SNPs. So "7/7 blocks and Gródek rejected" is **not** evidence that PW descends
from no medieval population — read it as the method detecting the data type.

**The one contrast with both confounders controlled.** EGermany_SP was cut into disjoint
blocks matched to the comparison group's size, in a single cache:

| Pair | n each | SNPs | p HO | p 1240K |
|---|---|---|---|---|
| PW + Poland_SP | 15 | 314K / 612K | 0.113 | 0.023 |
| PW + EGermany blocks | 15 | 314K / 612K | median 0.008 | median 0.001 |
| PW + Gródek area | 27 | 339K / 659K | 0.0127 | 0.0009 |
| PW + EGermany blocks | 27 | 339K / 659K | median 0.0013 | median 0.0000 |

At matched n *and* matched SNPs, PW is **consistently about an order of magnitude less
distinguishable from the Slavic-period population around Gródek than from the East
German one**. Both rejections are significant, so this is a difference of degree, not of
kind, and it does not license a ranking. It is, however, the only comparison in this
whole series where both main confounders were held fixed at once — weak, but controlled.

Set beside the f4 result it is coherent rather than contradictory: PW is **excluded from
the local MP→SP continuity at Gródek** (z −3.72 / −4.45) yet is **least distinguishable
from the Slavic-period phase there**. Those two say the same thing as the MP→SP
convergence table — PW shares the component that arrived in the Slavic period, but not
the local pre-Slavic substrate.

### A finding about the dataset, not about PW

**`Poland_SP` and `Pidhirtsi` are one population under two country labels.** Clade at
p 0.281 (HO) / 0.287 (1240K), consistent on both panels, and f4(Mbuti, PW; Poland_SP,
Pidhirtsi) is n.s. (z −0.77 / −2.00). They are recorded 11 km apart with overlapping
dates, one as Poland and one as Ukraine. Merge them (n=27) under a geographically
neutral name — `BugRegion_SP` — rather than either country label.

It follows that **`NWUkraine_SP` as shipped mixes two populations** separated by 277 km
and 150–200 years: Pidhirtsi (n=12, 1050–1150 CE, beside Gródek) and Korolivka (n=13,
1250 CE, far south-east). Anyone using `2025_Gretzinger_SlavicPeriod` should split it
before treating it as a geographic unit; a Poland_SP / NWUkraine_SP contrast is not an
east–west axis.

### Register of retractions

Three independent routes each produced an apparent ranking among the Slavic-period
groups, and each dissolved under control:

| Route | Apparent result | What dissolved it |
|---|---|---|
| PCA centroid distance | EGermany_SP nearest PW | centroid SE scales with n (0.0027 at n=260 vs 0.0065 at n=15); median distance to members reverses the order into a tie |
| f4 vs NWUkraine_SP | PW nearer western groups, z −3.2/−3.9 | the signal came from Pidhirtsi, which is 11 km from Gródek, not eastern; vs Korolivka it is n.s. on both panels |
| qpWave | PW separable from EGermany_SP, clade with Poland_SP | p tracks SNP count; at high power everything separates, largely on WGS-vs-capture differences |
| Trzciniec + X scan | seven second sources fit at 20–40 % | target n=1; against the medieval regions (n=18–48) the same models pass 0/10 |

That three different methods agree the direction is unresolvable is a stronger result
than any one of them alone.

**Net result: no directional claim is supportable.** PW sits with the Slavic-period
cluster as a whole; which part of it cannot be resolved with this data.

## Closing the qpAdm line: retargeting onto the medieval regions

PW is a single genome. Every model comparison above therefore runs at n=1, where qpAdm
has little power to reject anything. Since set 3 found PW statistically indistinguishable
from the early medieval populations of Greater Poland, Masovia and Lubusz, those
populations can be used as the target instead — n=18–48 rather than n=1, for no new data,
and a conclusion about them is a conclusion about PW.

### First pass: Santok alone

The cheapest version of the idea needs no new fileset. `Poland_EarlyMedieval_Slav`
(Santok, **n=46**) was already in the set-6 fileset, so the `Trzciniec + X` scan could be
re-run against it with the same 14 candidates, the same panels and the same right set as
for PW:

| Target | n | HO passing | 1240K passing |
|---|---|---|---|
| **PW** | 1 | **13 of 14** | 3 of 14 |
| **Santok** | 46 | **0 of 14** | 1 of 14 |

The tie collapses. The same sources that fit a single genome 13 times over fit a
46-individual population from the same region **zero** times.

The single survivor on 1240K was `Trzciniec 0.64 ± 0.06 + EGermany_MP 0.36 ± 0.06`
(p 0.052), which on HO was also the top non-degenerate model (p 0.040, just below
threshold). Everything else was rejected outright, `Poland_SP` (Gródek) included at
p 0.0028 / 0.00016 — consistent with the Gródek-excludes-PW result above. That survivor is
the lead recorded at the end of this section; it is not a finding.

### Second pass: all regions — and a warning about the first

Two near-identical targets give opposite answers. `Poland_EarlyMedieval_Slav` (n=46,
Santok) and `Reg_Lubusz_Santok` (n=48, Santok + Cedynia) **share 46 of their individuals**
— they differ by two burials — yet:

| Target | n | shared | HO passing |
|---|---|---|---|
| `Poland_EarlyMedieval_Slav` (set 6, 14 candidates) | 46 | 46/48 | **0/14** |
| `Reg_Lubusz_Santok` (set 7, 10 candidates) | 48 | 46/48 | **9/10** |

What changed besides those two individuals: the candidate pool (14 vs 10), hence the
rotating right set, hence the SNP count under `maxmiss=0` (167,807 vs 187,793 on HO).
This is the set-4-versus-set-5 lesson again, now with almost the same people on both
sides: **the verdict tracks the design, not the past.** Neither number should be read as a
property of Santok. The regional scan below is reported for the same reason — as an
illustration of how the method behaves, not as a result about these populations.

`set7` = the set-3 regional fileset (medieval Polish individuals relabelled by site) with
the Gretzinger groups merged in by `set6/40_add_external.R`, pointing `FROM_SET` at set3.
Scan: `Trzciniec + X` for ten candidates, right = 10 deep outgroups plus the unused
candidates, a model counted as passing only at p > 0.05 with all weights > 0 and every
SE < 0.25.

| Target | n | HO SNPs | HO pass | 1240K SNPs | 1240K pass |
|---|---|---|---|---|---|
| **PW** | 1 | 201,955 | **7/10** | 366,501 | **3/10** |
| Reg_Masovia | 18 | 162,922 | 7/10 | 294,910 | 2/10 |
| Reg_GreaterPoland_East | 26 | 190,414 | 1/10 | 345,599 | 0/10 |
| Reg_GreaterPoland_Poznan | 30 | 130,743 | 1/10 | 236,844 | 0/10 |
| Reg_Silesia | 34 | 199,902 | 0/10 | 362,509 | 0/10 |
| Reg_GreaterPoland_Gniezno | 34 | 201,175 | 0/10 | 364,887 | 0/10 |
| Reg_Lubusz_Santok | 48 | 187,793 | 9/10 | 340,486 | 0/10 |

**The SNP count is part of the power, not a constant.** Under `maxmiss=0` each target gets
its own cache, so coverage of the target group sets how many sites survive — the same
trap that invalidated the qpWave table. On 1240K the pass rate is monotonic in n × SNPs
(PW 0.4 M → 3 passes; Masovia 5.3 M → 2; everything at 7.1 M or above → 0).

**Masovia's survival is what the power account predicts, not a signal.** It is the
smallest region tested (n=18) and therefore the weakest at rejecting. Its agreement with
the PCA/nMonte Płońsk result is tempting and must be resisted: it is the same shape as
`Poland_SP` looking "closest" in qpWave and `Latvia_MA` fitting everything in set 6 —
the least powerful group looking the most compatible.

**A subsampling control splits the verdict.** Taking the two full-rejecting regions down
to Masovia's n=18:

| Region | n | HO SNPs | HO pass | 1240K pass |
|---|---|---|---|---|
| Silesia | 34 → 18 | 200K → 139K | 0/10 → **6/10** | 0/10 → 0/10 |
| Gniezno | 34 → 18 | 201K → 187K | 0/10 → **0/10** | 0/10 → 0/10 |
| Santok | 48 → 18 | 188K → 138K | 9/10 → 9/10 | 0/10 → 0/10 |

Silesia behaves exactly as the power account predicts. **Gniezno does not** — at n=18 with
its SNP count barely moved (−7 %, against Silesia's −31 %) it still rejects all ten. So
power explains some of the regional pattern and not all of it, and the control is itself
unequally powered. On 1240K nothing passes at any n.

**The obvious alternative — that Gniezno is internally heterogeneous — was tested and does
not hold.** qpAdm models the target as one population, so a mixed group (plausible for a
capital with elite burials) would reject every simple model regardless of n. Dispersion in
the existing 25-PC space, mean distance of members to their own region centroid:

| Region | n | mean | sd |
|---|---|---|---|
| Poznań | 17 | 0.0526 | 0.0119 |
| Santok | 30 | 0.0511 | 0.0138 |
| Silesia | 26 | 0.0445 | 0.0124 |
| Greater Poland East | 24 | 0.0398 | 0.0111 |
| **Gniezno** | 24 | **0.0362** | 0.0122 |
| **Masovia** | 9 | **0.0295** | 0.0060 |

Gniezno is among the *tighter* regions, and the decisive comparison runs the wrong way:
Silesia is more dispersed yet starts passing when subsampled, while the tighter Gniezno
does not.

**What this dispersion measure does not isolate:** distance to a centroid in 25 PCs
carries projection noise, which grows as an individual's coverage falls. The pair the
verdict rests on is exactly where coverage differs most — subsampling cost Silesia 31 % of
its SNPs (its total was carried by a few well-covered samples) against Gniezno's 7 %
(evenly well covered). So part of Silesia's larger spread may be projection noise from its
weaker individuals rather than real heterogeneity. This does not reverse the result —
Gniezno is still not dispersed, and the hypothesis still fails — but it weakens the one
contrast that decided it. Masovia is indeed the most homogeneous real region, which fits it being the
easiest to model — but that does not account for Gniezno. (PCA n differs from qpAdm n:
the PCA applies its own coverage floor.)

**Gniezno's residual is therefore recorded as open, not interpreted.** It is one small
unexplained inequality in a section whose conclusion — that f-statistics have stopped
discriminating here — does not depend on it.

### The result worth carrying outside this project

**The populations PW is statistically indistinguishable from reject the models that fit
PW.** Seven of ten candidates pass with PW as target on HO; the same ten against Gniezno,
Silesia, Poznań or Greater Poland East pass none. This is a *measured* demonstration —
not an assumption — that a qpAdm target of n=1 does not discriminate between models. It
generalises well beyond this dataset, and it is a standard move in amateur and
professional aDNA analysis alike.

### What "no model fits" does and does not mean

qpAdm tests the null that the target is *exactly* a mixture of the named sources. For any
real population that null is false, so with enough power every simplification is rejected.
The useful regime was always the low-power one, where qpAdm discriminates *between*
models; at high n and high SNP counts it only reports "not exactly this", which is
trivially true.

So the correct conclusion is **not** that the history needs three or more sources. It is
that **f-statistics have stopped being a discriminating instrument for this question**.
Pushing further — more candidates, more source combinations — cannot recover resolution
that the method no longer has. The qpAdm line is closed here.

What would move the question is a different class of evidence: haplotype sharing. Of the
590 Gretzinger individuals, 266 carry ≥600,000 1240K sites and 391 carry ≥200,000, which
is within the usual working range for imputation followed by IBD detection. The data-type
asymmetry that was a liability above becomes an asset there — PW's diploid 30x WGS is an
ideal IBD target, where in f4 it was the source of a systematic offset. IBD also asks the
question f-statistics could not pose: not what PW's admixture proportions are, but whether
PW shares long segments with specific individuals, and how long — because segment length
dates the common ancestor. The cost is the ENA reads (`PRJEB81250`) plus an imputation
pipeline; the Poseidon package ships genotypes, not reads.

## What this adds

- The Migration-Period population of south-eastern Poland is **excluded** as PW's
  ancestor, more strongly than any other Polish group, by qpAdm and by PCA.
- Gródek shows **local continuity** across the Migration→Slavic boundary that PW is
  outside of — eastern Poland has its own line, and it does not lead to PW.
- The Slavic-period shift is **real and multi-regional**, visible as three independent
  populations converging on PW's position.
- Within the Slavic period, **no sub-direction is resolvable**: the apparent
  west-versus-Volhynia contrast collapses once NWUkraine_SP is split by site, and
  qpWave puts PW closest to Poland_SP rather than EGermany_SP.

## What this does not show

- **Not a donor.** Seven second sources fit Trzciniec equally well at 20–40 % — and
  retargeting onto the medieval regions (n=18–48) shows that tie was an artefact of PW's
  n=1: the same models pass 0/10 against Gniezno, Silesia, Poznań and Greater Poland
  East. The honest statement is that no two-source model is identifiable here at all.
- **Register of retractions, continued.** A fourth apparent result — "PW ≈ Trzciniec +
  one source, seven candidates" — dissolved the same way as the other three, under the
  same cause: a comparison run at a power the data could not support.
- **Not a date.** qpAdm cannot distinguish an admixture event from an unsampled
  source population that already had the mixture.
- **p is design-dependent** in the 0.05–0.2 band: the same models flip between set 4
  and set 5 purely from the rotating right set. Only |z| ≥ 3 f4 results and
  rejections by orders of magnitude are treated as findings here.
- **Not a west–east axis.** The SP groups differ in date (1074–1100 vs 1250), in
  sample size (15 vs 260) and, for NWUkraine_SP, in internal composition. Any ranking
  among them mixes those three things with ancestry. Only the MP/SP contrast, which is
  replicated in three independent regions, is treated as a finding here.
- **`Trzciniec + Poland_SP` is circular, not a donor model.** Poland_SP is itself an
  already-admixed population of 1074 CE, i.e. younger than the profile it is meant to
  explain (set 3 places PW's matches at 900–1200 CE). Identical weights on both panels
  show consistency, not validity. It is useful only as a diagnostic — "was the
  EEF-richer component already present in an 11th-century Slavic population" — and is
  disqualified as a source for exactly the reason Great Moravia was.
- **The regional gap persists.** Every Polish MP/SP sample in this package is from
  **south-eastern** Poland (Gródek, on the Bug). Greater Poland and Masovia before
  ~900 CE — PW's own region, where set 3 places his medieval matches — remain
  unsampled. The Polish Iron Age was cremating its dead; that is a gap in the record,
  not in the analysis.

## Caveats on the new data

- Poland_MP n=14, Poland_SP n=15, NWUkraine_SP n=25.
- The SP groups are **site-concentrated, not regionally spread**: EGermany_SP's four
  sites lie within 23 km of each other in Saxony-Anhalt; Poland_SP is a single site.
  No regional sub-analysis is possible from this package.
- `Pidhirtsi (Podhorce)` is recorded in the `.janno` at 50.78 N, 23.78 E — 11 km from
  Gródek and on the Polish border — while the well-known Pidhirtsi lies ~150 km away
  near Lviv. Homonym or coordinate error, unverified against the paper. This is not a
  footnote: splitting NWUkraine_SP on it is what dissolved the apparent west–east
  signal, so the site's true location decides how that section should read. If the
  coordinates are wrong and Pidhirtsi is in fact near Lviv, the split must be redone.
- **Target and references differ in data type.** PW is diploid 30x WGS; the ancient
  groups are capture data, largely pseudo-haploid, with deamination and reference bias.
  That difference is systematic and becomes detectable as apparent structure at high SNP
  counts — see "Why qpWave does not settle this question". Any statistic comparing PW
  directly with an ancient group carries it; only contrasts *between* two ancient groups
  against PW (f4 of the form `f4(Mbuti, PW; A, B)`) cancel most of it.
- `NWUkraine_SP` and `Poland_SP` are not independent groups — see "A finding about the
  dataset, not about PW".
- ~400 two-source models were tested in one pass: a p of 0.05–0.2 is noise.

## Sources

- Gretzinger J. et al. 2025. Ancient DNA connects large-scale migration with the
  spread of Slavs. *Nature*, doi:10.1038/s41586-025-09437-6. Genotypes: Poseidon
  Community Archive, `2025_Gretzinger_SlavicPeriod` v1.0.0.
- Saag L. et al. 2025. North Pontic crossroads. *Science Advances* (AADR key
  `SaagThomasSciAdv2025`).
- Krzewińska M. et al. 2018. *Science Advances* (`Moldova_Scythian`, Hlinaia).
- Chyleński M. et al. 2023. *Nature Communications*,
  doi:10.1038/s41467-023-40072-9 (Trzciniec).
- Stolarek I. et al. 2023. *Genome Biology*, doi:10.1186/s13059-023-03013-9.
- Mallick S. et al. 2024. The AADR. *Scientific Data* 11:182.

Poseidon packages: list at `https://server.poseidon-adna.org/packages`; genotype
files via Git LFS at
`https://media.githubusercontent.com/media/poseidon-framework/community-archive/master/<pkg>/<file>`
(plain raw.githubusercontent returns LFS pointers).
