# Set 3: the target vs medieval Polish REGIONS

AADR labels follow publications, not geography (`Poland_IA` = 130 Piast-era
burials 900-1200 CE from 19 sites; `Poland_EarlySlav` = Ląd + Łęgowo;
`Poland_EarlyMedieval_Slav` = Santok). Set 3 regroups every medieval Polish
individual (400-1450 CE, outliers excluded) by site into historical regions
(`regions.R`), splits `Poland_BA` by date, and tests the target against them.

| Script | What |
|---|---|
| `regions.R` | region -> site map, period splits, reference periods, outgroups |
| `20_build_f2.R` | reads the individuals, relabels, merges, measures per-group coverage, builds the f2 cache |
| `21_region_tests.R` | A: are regions distinguishable at all (qpWave); B: pairwise f4(Mbuti, T; A, B); C: outgroup f3; D: rotating qpAdm over regions + periods |

```bash
cd /path/to/personal-qpadm-pipeline   # run dirs are created here
M=./PW_WGS_HO/merged_pw_aadr_final_ready
A=$AADR   # converted-GENO prefix, see ../../README.md §2
Rscript aadr/set3/20_build_f2.R $M $A PW && Rscript aadr/set3/21_region_tests.R $M $A PW
```

## Strict mode only (default)

Groups covering < `SPARSE_FRAC` (0.6) of the target's SNPs are dropped from the
joint cache and listed as NOT TESTED (Lesser Poland, the 3 late-BA individuals,
Przeworsk, Lithuania_LBA on PW WGS). The rest must all be present at every SNP.

**Do not use `ALLSNPS=1` (pairwise overlap, maxmiss=1) for conclusions.** Tested
2026-10-22 on PW WGS: it rejected almost every source with absurd p (Lithuania_LBA
1e-169, Czechia_EarlyMedieval_EarlySlav 1e-97) and flipped verdicts between HO and
1240K (Gniezno 0.44 vs 0.047, Silesia 0.13 vs 0.002). Different SNP sets per pair
break qpAdm's f4 covariance for a diploid WGS target against sparse pseudo-haploid
groups. The mode is kept only to reproduce that diagnosis.

Array targets (consumer chips, ~57K panel SNPs) keep ~10K SNPs here and almost
everything fits -- too little power for regional questions.
