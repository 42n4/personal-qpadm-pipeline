# personal-qpadm-pipeline — Claude Code project context

qpAdm / f-statistics analysis of personal genomes (consumer arrays or WGS)
against the Allen Ancient DNA Resource (AADR v66.p1). Two variants share one
code base and this file:

- `aadr/` — Poland and its neighbours (Slavic, Baltic, Lusatian, Wielbark,
  Przeworsk, steppe sources). The full pipeline: set 1 (steps 00–09), WGS input
  (01b), set 2/3/4, PW-vs-JW pair test, own PCA.
- `aadr_eu/` — West/Central-European fork of steps 00–03, 07, 08 for
  admixed or non-Polish targets. It retains every `aadr/` population and adds British/Anglo-Saxon,
  continental Germanic, Gaulish/French, Italic/Roman, Iberian and modern
  English/French/Italian/Spanish/Orcadian/Icelandic anchors. Steps 04–06, 09 are
  Polish-target probes: run them from `aadr/`.

Targets are named by short labels only: **PW** (FamilyTreeDNA array and 30x WGS),
**JW** (MyHeritage array), **AM** (public PGP-HMS FamilyTreeDNA file). No personal
genotypes, file names or genealogy belong in this repository.

**Results — read before any new analysis, do not re-derive:**
`aadr/examples/PW_WGS_aadr.md` (WGS, sets 1–4, pair test, PCA),
`aadr/examples/PW_WGS_aadr_gretzinger.md` (sets 5-6 + Gretzinger 2025 PCA),
`aadr/examples/PW_JW_arrays_aadr.md` (arrays, set 1),
`aadr_eu/examples/AM_FF_aadr_eu.md` (AM, aadr_eu).

## TGENO: every result before 2026-10-21 was computed on garbage

AADR v66.p1 `.geno` files are **TGENO** (transposed: a fixed 48-byte header, then
one record of ceil(nsnp/4) bytes per INDIVIDUAL). admixtools 2.0.10
`read_packedancestrymap()` reads only SNP-major GENO and decodes TGENO silently into
garbage: HO missingness 38.3 % instead of 3.4 %, Han heterozygosity 0.077 against
a Hardy-Weinberg expectation of 0.374, Mbuti closer to Han than to Yoruba.

- Convert once: `python3 -I aadr/tgeno2geno.py <TGENO.geno> <out.geno>` (HO ~2.5 min;
  output records are max(48, ceil(nind/4)) bytes per SNP), then symlink the
  `.snp/.ind/.anno` next to it and always use that prefix.
- Check every new AADR download: `head -c 6 X.geno` must print `GENO`, never `TGENO`.
- Garbage signature in an existing f2 cache: f2(Mbuti, Han) < f2(Mbuti, Yoruba).

The June 2026 "closed findings" (Ukraine_IA_Lusatian best source p≈0.9–1.0, JW drifts
south, Avar/G25 artefacts "falsified by step 8") were all computed on mis-decoded data.
The table below replaces them.

## Findings on correct data (2026-10-22)

| Question | Answer | Where |
|---|---|---|
| Can set 1 (Mbuti/Yoruba/Han outgroups) rank European sources? | **No.** On arrays 33/38 sources in step 07 and every Polish source in step 06 fit; on PW WGS 26/38. It only rejects East-Asian/Siberian, Near-Eastern, Ashkenazi and pure early-farmer sources | arrays + WGS, steps 06–07 |
| Germany/CWC component (Germany_IA proxy) | **Not supported**: weight never significant, 6/6 outgroup sets (arrays, WGS HO), 5/5 (WGS 1240K) | steps 04–05 |
| Hidden second source (steppe, Siberian, Near-Eastern, Ashkenazi, Avar) | **None**: 0 SIGNIFICANT-2src on either array; such second sources are rejected outright | step 08 |
| Best single sources for PW (WGS, ancient outgroups) | Poland_IA (= Piast-era 900–1200 CE), Poland_EarlySlav, Poland_EarlyMedieval_Slav, Lithuania_LBA fit on HO and 1240K. Ukraine_IA_Lusatian, Wielbark, Przeworsk, Trzciniec, Corded Ware, Sweden_IA, Avars, Hungarian Conquerors are rejected | set 2 |
| Medieval Polish regions for PW | Greater Poland East (Ląd), Gniezno, Poznań, Masovia (Płońsk), Lubusz (Santok), Czechia EarlySlav fit. Silesia, Kuyavia, Bodzia are rejected. Regions are barely distinguishable from each other | set 3 |
| Pre-medieval sources only | No single source fits PW. Pairs that fit = East-Baltic LBA/EIA ~30–45 % + Central/Northern European IA ~55–70 %. Best: Przeworsk 0.6 + Latvia_LBA_FlatGraves 0.4, which also fits JW | set 4 |
| Distal PW | WHG ~5 %, Steppe ~60 %, EEF ~35 % (SE 0.05–0.08) | set 2, 12 |
| PW vs JW | **Indistinguishable**: qpWave p 0.48–0.86, 0/63 f4 with \|z\| ≥ 3 | `pair_compare.R` |
| Przeworsk vs Wielbark channel | Indistinguishable on arrays; both rejected as sole sources on WGS | step 09, set 2 |
| Own PCA + nMonte (G25-style) | PW ≈ medieval Płońsk 69 % + Santok 16 % + CW/BA 15 %; JW ≈ Ląd 51 % + Płońsk 41 %. Agrees with set 3 | `aadr/pca/` |
| AM (aadr_eu) | European, no continental admixture (f3 24/24 n.s.). On 43K SNPs Germany_AltInden_EarlyMedieval fits alone (p 0.53); England IA/Saxon, French IA, Spanish, Italian_North, Poland_EarlySlav bases are rejected alone (p 0.006–0.028); 0 SIGNIFICANT-2src of 401 | aadr_eu 07, 07b, 08 |

Read p-values across sources as fit / no fit, never as a ranking.

## Repository layout

```
CLAUDE.md, README.md, LICENSE, .gitignore   # shared by both variants; .gitignore is deny-by-default
aadr/
  00_merge_arrays.sh         union-merge two arrays of one person before step 01
  01_prepare_files.sh        consumer CSV (FTDNA / MyHeritage / Ancestry) -> PLINK
  01b_prepare_wgs.sh         WGS CRAM/BAM -> PLINK genotyped at the AADR panel sites
  tgeno2geno.py              AADR TGENO -> GENO converter (run once per panel)
  02_build_aadr_subset.R     AADR population subset -> PLINK (memory-safe pops= filter)
  03_merge_and_test.sh       merge target + subset -> merged_pw_aadr_final_ready
  04_run_qpadm_test.R        set 1: Germany/CWC component test (model A/B + popdrop)
  05_robustness_outgroups.R  set 1: the same over 6 outgroup sets
  06_ancestry_models.R       set 1: 18 Poland_* models (Poland_IA = MEDIEVAL despite the label)
  07_extended_models.R       set 1: 38 single sources outside Poland_*
  08_two_source.R            set 1: rotating 2-source scan, SIGNIFICANT-2src detector
  09_przeworsk_channel.R     set 1: Przeworsk vs Wielbark (1-source, f4 symmetry, 2-source)
  set2/                      ancient outgroups: build f2, outgroup f3, distal, rotating qpAdm, f4 scan
  set3/                      medieval Polish regions by site (strict coverage mode)
  set4/                      pre-medieval candidates only (pops.R for the set-2 scripts)
  set5/                      set 4 + Ukrainian forest-steppe / Chernyakhiv-era candidates;
                             30_trzciniec_tests.R (f4 scan, Gustorzyn split, qpAdm, 3-source)
  set6/                      40_add_external.R: merge an external PLINK package (Poseidon)
                             into a set-N fileset, then SET=set6 runs the set-5 tests
  pair_compare.R             two targets on shared SNPs (qpWave clade, f4 scan, distal)
  pca/                       own 25-PC PCA (30_build_pca.py; --extra-plink projects every
                             individual of a PLINK set) + mapping to G25 catalogue ids (32)
  examples/                  result reports + README_commands.md (full command history)
aadr_eu/
  00-03                      same as aadr/
  07_extended_models.R       + West/Central-EU batch, right = Mbuti Yoruba Han Papuan Karitiana
  07b_continental_f3.R       admixture f3(T; West European, Yoruba/Han/Karitiana/Papuan) on the 07 cache
  08_two_source.R            + EU batch; hardened SIGNIFICANT-2src rule (below)
  examples/                  AM report, pointer to the PGP-HMS public genomes
```

## Invocation

Run from the repository root. Every run directory (`./PW_FF`, `./PW_WGS_HO`, ...)
is created under the current directory and is git-ignored; AADR stays read-only.
`$AADR` = a converted-GENO prefix, e.g. `/data/aadr/v66.p1/geno/v66.p1_HO.aadr.patch.PUB`.

```bash
SAMPLE_ID=PW ./aadr/01_prepare_files.sh /path/to/raw.csv PW_FF   # or 01b for WGS
Rscript aadr/02_build_aadr_subset.R $AADR ./PW_FF
./aadr/03_merge_and_test.sh ./PW_FF PW                           # arg 2 = SAMPLE_ID
Rscript aadr/0N_*.R ./PW_FF/merged_pw_aadr_final_ready $AADR PW  # 05, 07, 08, 09
Rscript aadr/06_ancestry_models.R ./PW_FF/merged_pw_aadr_extended $AADR PW
Rscript aadr/04_run_qpadm_test.R ./PW_FF/merged_pw_aadr_final_ready PW
```

Sets 2–4, WGS, pair test and PCA: `aadr/examples/README_commands.md` and the
README in each `set*/` directory.

- Step 03 skips work if `merged_pw_aadr_final_ready.*` exists.
- f2 caches land in the run directory (`dirname(prefix)`), never in CWD.
- **WGS: use HO for set 1** (1240K v66.p1 lacks 14 modern labels the scripts use:
  Ukrainian, Belarusian, Lithuanian, Romanian, Jew_Ashkenazi, Mordovian, Chuvash,
  Khanty, Nganasan, Selkup, Mongol ...). Set 2–4 work on both panels; read only
  results that hold on both.

## Hard-won pitfalls — pre-fixed in the scripts, do not regress

Numbers 1–16 are cited by number in script comments; keep them stable.

1. plink2 needs `--pedmap`, not bare `--file`. (01, 3a)
2. plink2 cannot read packedancestrymap `.geno/.snp/.ind`. (02)
3. `packedancestrymap_to_plink()` loads the whole panel (~129 GB) — OOM. Use
   `read_packedancestrymap(pops=...)`. (02)
4. admixtools has no exported `write_plink()` — PED/MAP written by hand. (02)
5. `read_packedancestrymap()$ind` columns are `X1/X2/X3`; access positionally. (02)
6. SNPs 100 % missing in a subset get `A1=A2="0"` — excluded before merge. (3b)
7. plink2 `--pmerge` cannot do non-concatenating merges — use `p-link --bmerge`
   (plink 1.9, binary name `p-link`). (3e)
8. plink 1.9 `--sort-vars` cannot combine with `--make-bed` — go through pgen. (3c)
9. "3+ alleles" from one-sided missingness or an array miscall (A/C vs T/C) is not
   a strand problem; `--flip` does not fix it. Exclude the SNP IDs from both sides
   and merge again. Steps 03, 07, 08 (both variants), set 2 and `pair_compare.R` do this.
10. [resolved] plink 1.9 v1.90b2p segfaulted at ~99 % of `--bmerge`; use ≥ v1.90b7.
11. `extract_f2(auto_only=TRUE)` cannot parse X/Y/MT strings — `plink2 --chr 1-22`
    first (3f, and on every `--make-bed` in aadr_eu 07/08).
12. After merge the `.fam` FID is the sample ID; step 3g rewrites it to the
    population (`aadr_subset_pop_lookup.tsv`), which `extract_f2`/`qpadm` match on.
13. f2 caches land in the run directory (`dirname(prefix)`), never in CWD; never run
    from a root-owned directory.
14. f2 caches of steps 07–09 are target-specific (`..._<TARGET>_f2cache/`); never
    share one between targets.
15. `options(width = 200)` at the top of every R step, or wide tables wrap.
16. **Never edit an R or shell step while it runs.** Rscript reads the file
    incrementally; an edit mid-run ends in `unexpected end of input`.
17. **WGS (01b):** genotype only at panel sites (all WGS positions break the
    f-statistic ascertainment). Do not use the provider's variant-only VCF (no rows
    at hom-ref sites). `bcftools call -m` emits no GQ at invariant sites — 01b uses
    3 × ref reads. plink2 writes "." as A1 at homozygous sites — 01b fills it from
    the panel pair, otherwise pitfall 9 drops ~77 % of SNPs. REF must be the FASTA
    the CRAM was aligned to (contig lengths are checked; CRAMs often lack M5 tags).
18. **Array targets in set 2:** keep only target ∩ panel SNPs (10_build_f2.R does);
    otherwise the array's off-panel SNPs enter the joint fileset.
19. **Set 3 strict mode:** groups covering < 60 % of the target's SNPs are dropped
    from the joint cache (a sparse late-BA group collapsed it to a few thousand
    SNPs). `ALLSNPS=1` (pairwise overlap) gives absurd p-values and flips between
    panels — diagnosis only.
20. **qpAdm weights:** set 2 reports a FIT whose weight SE > 0.25 as UNIDENTIFIABLE.
21. **aadr_eu memory (07/08):** SNP pre-filter to the target's genotyped set
    (`*_dense`), `--chr 1-22`, and `extract_f2(maxmem = 1000)`, which must be below its own
    no-split estimate or it does not split. The PED matrix of the extra
    populations dominates RAM (aadr 07: ~36 GB, aadr_eu 07: ~52 GB peak on HO).
22. **aadr_eu `extra_needed` must include `right_pops`**, or new outgroups (Papuan,
    Karitiana) are never merged ("Populations missing in indfile"). A new right pop
    forces an ext7/2src rebuild (delete `merged_pw_aadr_{ext7,2src}*`).
23. **aadr_eu 08 SIGNIFICANT-2src rule:** base n ≥ 15 (`MIN_BASE_N`; low-n bases are
    demoted to second sources), p > p_base, both weights in (0.1, 0.9), |z| > 3 and
    Benjamini-Hochberg `q_second < 0.05`. The lax rule produced 25 false hits, all on
    the n=3 Ukraine_IA_Lusatian base. `aadr/08` still has the lax rule (0 hits on
    the arrays anyway).

## SNP-count expectations

- Consumer array: ~57K of 584K HO sites typed. After `maxmiss=0`, ~57K remain in
  steps 04/05, ~48K in 06, ~45K in 08, ~20K in 07/09, and ~7–10K in set 2–4
  (sparse ancient references). A union of two arrays of one person adds nothing.
- WGS 30x: 578,836 / 579,720 autosomal HO sites typed (1,147,026 on 1240K). After
  `maxmiss=0`: set 2 HO 46,927 / 1240K 82,172; set 4 HO 66,844 / 1240K 115,898.
  The limit is now ancient missingness, not the target.

## AADR label traps

- **`Poland_IA` is not Iron Age**: all 130 individuals are dated 900–1200 CE
  (early Piast-state cemeteries, Stolarek et al. 2023, Genome Biology). Two Markowice
  samples say "850-1050 BCE" in the text field but 1000 ± 58 BP in the numeric one
  (a BCE/CE typo).
- No Polish Lusatian-culture genomes, and no Zarubintsy, Kyiv, Pomeranian or
  Milograd groups. `Ukraine_IA_Lusatian` = 3 individuals from Volhynia.
- No `Germany_IA` matching Davidski's Saxony-Anhalt cluster; `Germany_Esperstedt_CordedWare`
  is the proxy (tests the CWC substrate, not IA Germanic ancestry).
- `Poland_IA_Przeworsk` n = 2 (both female, Gąski, 100–300 CE); steps 06–07 skip
  n < 3, step 09 and sets 2/4 admit it.
- Labels follow publications, not geography: set 3 regroups medieval Poland by site.
- Y-haplogroups (`.anno`, YFull auto-calls): every ancient male tends to get his
  own terminal branch, so n = 1 per exact branch is normal.

## Tools

- R `admixtools` (uqrmaie1, `remotes::install_github`), `plink2`, `p-link`
  (plink 1.9 ≥ v1.90b7), `bcftools`/`htslib` for 01b (tested with 1.24), Python 3 + numpy for
  `tgeno2geno.py` and `pca/`.
- AADR v66.p1 from Harvard Dataverse `doi:10.7910/DVN/FFIDCW` (file API; no R
  download helper). Dataverse md5sum files carry a leading line number: check with
  `awk '{print $2"  "$3}' | md5sum -c`.

## What not to do

- Do not run anything on an unconverted (TGENO) AADR prefix.
- Do not quote set-1 p-values as a ranking of European sources.
- Do not reintroduce `Ukraine_IA_Lusatian` as "best source" or "JW drifts south";
  both were artefacts of mis-decoded data.
- Do not quote weighted or mixed G25/nMonte components (Hungary_Mig_Himod,
  Latvia_Neo_Kivutkalns) as ancestry on their own: tail-PC amplification. But set 4
  finds the same shape (East-Baltic ~35–45 % + Central European), so the Baltic
  part is not plain noise.
- Do not add low-n populations as 2-source bases (aadr_eu 08 demotes them).
- Do not commit run directories, genotypes, raw files or anything that names a
  person — `.gitignore` is deny-by-default; re-include new code paths by name.
