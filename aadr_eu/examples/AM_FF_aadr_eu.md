# AM — public PGP-HMS FamilyTreeDNA array vs AADR v66.p1 HO (aadr_eu)

Rerun 2026-10-22 on **correctly decoded** AADR genotypes (see "TGENO" in
`../../CLAUDE.md`). The June 2026 AM results (continental f3 z ≈ +50 on every axis,
"25 → 0 SIGNIFICANT-2src", the French/Orcadian rejections) were computed on
mis-decoded data and are void.

Target: AM, a US participant of the Harvard Personal Genome Project who published
an autosomal FamilyTreeDNA file (`site.txt`). Right set for steps 07–08: Mbuti,
Yoruba, Han, Papuan, Karitiana. One SNP (rs114922848, A/C on the array vs T/C in
AADR) was an array miscall and is excluded at the merge.

| Step | SNPs after `maxmiss=0` | Runtime, peak RAM |
|---|---|---|
| 07 (58 single sources) | 20,118 | 9.6 min, 53 GB |
| 07b (continental f3) | same cache | seconds |
| 08 (8 bases × 50 second sources) | 43,445 | 7.3 min, 48 GB |

## 07b: is AM admixed with another continent?

Admixture f3(AM; E, X) for E = English, French, Spanish, Italian_North, Orcadian,
Icelandic and X = Yoruba (Africa), Han (East Asia), Karitiana (Americas),
Papuan (Oceania): **24 of 24 not significant**, z from −0.27 (English, Han) to
+1.32. No detectable non-European admixture, so a Europe-only left panel is adequate.
f3 is near zero with SE ~0.0015, so a few-percent component could still hide.

## 07: single sources

Rejected (p < 0.05): Austria_Avar, Hungary_EarlyAvar (both subsets),
Russia_Tatarstan_EarlyMedieval_EarlyBulgar, Jew_Ashkenazi,
Germany_Medieval_Jewish, Italy_Lazio_ImperialRoman_Roman (0.017), Finnish (0.020),
Ukraine_IA_Lusatian (0.0499).

The other 49 are not rejected, from England_Saxon and French through Spain_IA to
Belarusian and Sarmatians (highest p: Hungary_Madaras Sarmatian 0.82,
Germany_AltInden_EarlyMedieval 0.81, Belarusian 0.79). On ~20K SNPs with modern
outgroups this step excludes East-Asian-shifted, Near-Eastern-shifted and
Mediterranean-Roman sources. It does not rank European ones; same limit as PW and
JW (`../../aadr/examples/PW_JW_arrays_aadr.md`).

## 08: rotating 2-source scan

With fewer populations in the cache, step 08 keeps 43K SNPs instead of 20K, and the
bases alone separate more:

| Base alone | p |
|---|---|
| Germany_AltInden_EarlyMedieval (Lower Rhine, early medieval) | 0.525 |
| Italian_North | 0.028 |
| France_GrandEst_IA2 | 0.019 |
| France_Bucy_le_Long_IA_LaTene | 0.013 |
| England_IA | 0.011 |
| England_Saxon | 0.008 |
| Poland_EarlySlav | 0.007 |
| Spanish | 0.006 |

| Verdict (401 models) | n |
|---|---|
| SIGNIFICANT-2src | 0 |
| OK-base-sufficient | 1 (AltInden + Lezgin, second weight negative) |
| UNIDENTIFIABLE | 154 |
| REJECTED | 246 |

No second source improves on Germany_AltInden_EarlyMedieval alone. The 2-source
fits on the rejected bases are all UNIDENTIFIABLE, with weights such as −2.9/+3.9:
collinear sources, not admixture.

## Reading

- AM is European, with no detectable admixture from another continent.
- On 43K SNPs a continental North-West European early-medieval source
  (Germany_AltInden) fits alone. England IA/Saxon, French IA, Spanish,
  North Italian and Early Slavic bases are rejected alone at p 0.006–0.028.
  That points to a Rhine/North-Sea-continental profile rather than British,
  Iberian or Slavic. It is one test with modern outgroups only, so it is a lead, not
  a result.
- To separate North-West European sources properly, AM would need the set-2
  design (ancient outgroups, rotating qpAdm) with a West-European candidate list,
  and preferably WGS. Set 2 currently carries Polish/Baltic candidates only.
