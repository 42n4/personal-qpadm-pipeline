# personal-qpadm-pipeline — qpAdm ancestry analysis of a personal genome against AADR

This pipeline analyses your personal DNA (a consumer array file or a WGS CRAM)
against thousands of ancient genomes from the Allen Ancient DNA Resource (AADR)
using **qpAdm** and other f-statistics. The result tells you which ancient
populations are statistically compatible with your ancestry, and which are not.

This is a different, more rigorous approach than G25/nMonte (geometric distance).
qpAdm works on allele frequencies and can tell genuine shared ancestry from
geometric coincidence — provided the outgroups can see the difference (§3, set 1 vs set 2).

**Two variants, one code base:**

| Directory | For | Steps |
|---|---|---|
| `aadr/` | Poland and neighbours (Slavic, Baltic, Lusatian, Wielbark, Przeworsk, steppe sources) | 00–09 (set 1), 01b (WGS), set 2/3/4, pair test, own PCA |
| `aadr_eu/` | Western/Central European or admixed (e.g. US) targets | 00–03, 07, 08 with an added West/Central-European batch |

> **AADR v66.p1 warning.** Its `.geno` files are in the transposed TGENO format,
> which admixtools 2.0.10 decodes silently into garbage. Convert them first (§2).
> Every result this repository published before 2026-10-22 was computed on
> mis-decoded data; the reports in `aadr/examples/` and `aadr_eu/examples/` are
> the corrected reruns.

---

## 0. Scientific purpose

This repository contains an experimental research pipeline for analysing
consumer autosomal DNA using qpAdm and the Allen Ancient DNA Resource (AADR).

The project was developed as an independent scientific and educational study
to evaluate the applicability of formal f-statistics to personal ancestry
analysis. It is intended for reproducible research, methodological evaluation,
and learning.

No personal genomic data are included in this repository. Users must analyse
their own legally obtained genotype files (e.g. FamilyTreeDNA, MyHeritage or
AncestryDNA) and must obtain the AADR reference dataset separately according
to its licence terms.

The pipeline is provided for research purposes only. The results should not be
interpreted as medical, forensic, or genealogical proof without independent
validation.

---

## 1. Buy and download your raw DNA data

You need a **raw DNA file** (CSV format) from one of these providers:

### Family Tree DNA — Family Finder
1. Order at **familytreedna.com** → "Family Finder" autosomal test
2. After results arrive: **myFTDNA** → top-right menu → **"Download Raw Data"**
   → "Chromosome Browser Raw Data" → download the zip
3. Unzip: `unzip <filename>.zip` → you get a `.csv` file

### MyHeritage DNA
1. Order at **myheritage.com/dna**
2. After results arrive: **DNA** tab → **"Manage DNA kits"** → your kit →
   **"Download"** → "Download raw DNA data" → confirm by email link
3. Unzip: `unzip MyHeritage_raw_dna_data.zip` → you get a `.csv` file

### Ancestry DNA  *(partial support — step 01 auto-detects)*
Download via: **DNA** → **"Settings"** → **"Download DNA Raw Data"**

> **Tip:** If you have results from multiple providers for the same person,
> run step 01 separately for each, then merge with step 00 before continuing.
> See `aadr/examples/README_commands.md` for the merge workflow.

---

## 2. Prerequisites — software and data

```bash
# R package
Rscript -e 'remotes::install_github("uqrmaie1/admixtools")'

# Tools (Gentoo, overlay pwr; other distros: install plink2, plink 1.9, bcftools)
sudo emerge sci-biology/plink2-bin   # binary: plink2
sudo emerge sci-biology/plink        # binary: p-link  (plink 1.9, v1.90b7+)
sudo emerge sci-biology/bcftools     # only for WGS input (step 01b)

# AADR v66.p1 HO (and optionally 1240K) — Harvard Dataverse doi:10.7910/DVN/FFIDCW
# The .geno files are TGENO. Convert once (HO: ~2.5 min), symlink the rest:
mkdir -p /data/aadr/v66.p1/geno && cd /data/aadr/v66.p1/geno
python3 -I /path/to/personal-qpadm-pipeline/aadr/tgeno2geno.py \
    ../v66.p1_HO.aadr.patch.PUB.geno v66.p1_HO.aadr.patch.PUB.geno
ln -s ../v66.p1_HO.aadr.patch.PUB.snp ../v66.p1_HO.aadr.patch.PUB.ind ../v66.p1_HO.aadr.PUB.anno .
head -c 6 v66.p1_HO.aadr.patch.PUB.geno     # must print GENO, never TGENO

export AADR=/data/aadr/v66.p1/geno/v66.p1_HO.aadr.patch.PUB   # used below
```

---

## 3. Run the pipeline (`aadr/`)

Replace `NN` with your own short label (2–4 letters). Run everything from the
repository root: run directories (`./NN_FF`, `./NN_WGS_HO` ...) are created there
and are git-ignored.

### Step 01 — raw CSV to PLINK

```bash
SAMPLE_ID=NN ./aadr/01_prepare_files.sh /path/to/your_ftdna_rawdata.csv NN_FF   # FTDNA
SAMPLE_ID=NN ./aadr/01_prepare_files.sh /path/to/MyHeritage_raw_dna_data.csv NN_MH
```

Two arrays of the same person can be union-merged with `00_merge_arrays.sh`
(see `aadr/examples/README_commands.md`); it does not add usable SNPs.

### Step 01b — WGS instead of a CSV *(CRAM/BAM, GRCh37)*

```bash
REF=/path/GRCh37.primary_assembly.genome.fa.bgz JOBS=16 SAMPLE_ID=NN \
  ./aadr/01b_prepare_wgs.sh /path/NN.cram $AADR NN_WGS_HO
```

- Genotypes the reads at the panel's autosomal SNPs only (bcftools mpileup + call),
  so steps 02–09 run unchanged. Give step 02 the same AADR prefix.
- `REF` must be the FASTA the CRAM was aligned to; contig lengths are checked.
- Filters: `MIN_DP=8 MAX_DP=100 MIN_GQ=20 MIN_MQ=20 MIN_BQ=20` (env overrides).
- 30x WGS: ~579K of 580K HO sites typed in ~7 min with 16 jobs.
- Do not use the provider's variant-only VCF: it cannot tell hom-ref from missing.

### Steps 02–03 — AADR subset and merge

```bash
Rscript aadr/02_build_aadr_subset.R $AADR ./NN_FF
./aadr/03_merge_and_test.sh ./NN_FF NN          # 2nd argument = your SAMPLE_ID
```

### Set 1 — steps 04–09 (outgroups Mbuti, Yoruba, Han)

```bash
M=./NN_FF/merged_pw_aadr_final_ready
Rscript aadr/04_run_qpadm_test.R      $M NN          # Germany/Corded Ware component test
Rscript aadr/05_robustness_outgroups.R $M $AADR NN   # the same over 6 outgroup sets
Rscript aadr/06_ancestry_models.R ./NN_FF/merged_pw_aadr_extended $AADR NN   # Poland_* models
Rscript aadr/07_extended_models.R     $M $AADR NN    # 38 single sources outside Poland_*
Rscript aadr/08_two_source.R          $M $AADR NN    # rotating 2-source scan
Rscript aadr/09_przeworsk_channel.R   $M $AADR NN    # Przeworsk vs Wielbark (optional)
```

Set 1 rejects sources with gross non-European or early-farmer excess. It
**cannot rank European sources**: the three modern outgroups see Europe as
nearly one point, so dozens of European sources "fit" at once.

### Set 2 — ancient outgroups (the comparisons that discriminate)

```bash
for s in 10_build_f2 11_outgroup_f3 12_distal 13_rotating 14_f4_scan; do
  Rscript aadr/set2/$s.R $M $AADR NN || break; done
```

Outgroup f3, distal WHG/Steppe/EEF, rotating qpAdm (unused candidates become
outgroups) and an f4 scan. Population lists: `aadr/set2/pops.R`. Details:
`aadr/set2/README.md`.

### Set 3 — medieval Polish regions; Set 4 — pre-medieval sources only

```bash
Rscript aadr/set3/20_build_f2.R $M $AADR NN && Rscript aadr/set3/21_region_tests.R $M $AADR NN
export POPS_FILE=$PWD/aadr/set4/pops.R     # set 4 = the set-2 scripts with this list
for s in 10_build_f2 11_outgroup_f3 12_distal 13_rotating 14_f4_scan; do
  Rscript aadr/set2/$s.R $M $AADR NN || break; done
unset POPS_FILE
```

See `aadr/set3/README.md` and `aadr/set4/README.md`.

### Two targets on shared SNPs; own PCA

```bash
Rscript aadr/pair_compare.R <set3_or_set4_prefix> NN /path/other/pw_plink XX $AADR ./NN_vs_XX
python3 -I aadr/pca/30_build_pca.py --prefix $AADR --anno ${AADR%.patch.PUB}.PUB.anno \
    --target NN=./NN_WGS_HO/pw:wgs --relabel aadr/pca/relabel_regions.tsv --out ./NN_WGS_HO/pca
```

`30_build_pca.py` builds a 25-PC West-Eurasian PCA from AADR moderns (smartpca-style
normalisation) and projects ancients and targets by least squares (lsqproject).
Its distances agree with Davidski's G25 at r = 0.95. `32_to_g25_centroids.py`
turns it into centroids for a G25-style nMonte runner.

**WGS: run set 1 on HO** (1240K lacks 14 of the modern labels the scripts use); run
sets 2–4 on HO and 1240K and trust only what holds on both.

---

## 4. The `aadr_eu/` variant

Same steps 00–03; then:

```bash
Rscript aadr_eu/07_extended_models.R ./NN_FF/merged_pw_aadr_final_ready $AADR NN
Rscript aadr_eu/08_two_source.R     ./NN_FF/merged_pw_aadr_final_ready $AADR NN
```

- Adds England IA/Saxon/Medieval, Ireland, Netherlands, Germany (Alt-Inden),
  Hungary Lombards, France IA, Italy Roman/Piceni, Spain IA/EBA, and modern
  English, French, Italian_North, Spanish, Orcadian and Icelandic, to every `aadr/` population.
- Right set = Mbuti, Yoruba, Han, Papuan, Karitiana.
- Step 08 uses a strict `SIGNIFICANT-2src` rule: base n ≥ 15, improvement over the
  base, both weights in 0.1–0.9, |z| > 3 and a Benjamini-Hochberg `q_second < 0.05`.
- For an admixed target, check the continents:
  `Rscript aadr_eu/07b_continental_f3.R ./NN_FF/merged_pw_aadr_final_ready NN` (after 07)
  computes f3(NN; West European, X) for X = Yoruba, Han, Karitiana, Papuan.
  A negative estimate with z < −3 means real admixture with X; then the Europe-only
  left panel cannot model the target.
- Public genomes to try: PGP-HMS (`aadr_eu/examples/site.txt`).
- Peak RAM: ~52 GB on HO (the PED matrix of the extra populations).

---

## 5. Reading the results

| File | What it shows |
|------|---------------|
| `NN_FF/ancestry_models_results.tsv` | step 06, Poland_* models |
| `NN_FF/extended_models_results.tsv` | step 07, single sources |
| `NN_FF/two_source_results_NN.tsv` | step 08, second components |
| `NN_FF/przeworsk_channel_*_NN.tsv` | step 09 |
| `NN_FF/set2/*.tsv`, `set3/`, `set4/` | sets 2–4 |

- A model is compatible with the data when `p > 0.05`. A higher p is **not** a
  better ancestry match; compare sources only inside a design that can reject them.
- `UNIDENTIFIABLE` (step 08, set 2): the model fits but a weight SE is large —
  the sources are collinear and any weights fit. Not a finding.
- `OK-unphysical`: one weight is negative — collinear sources.
- `SIGNIFICANT-2src`: a real second component (step 08).

---

## 6. Example results: PW, JW (Polish), AM (US public genome)

Short labels only; no personal data is in this repository.

- `aadr/examples/PW_WGS_aadr.md` — PW 30x WGS: sets 1–4, PW vs JW, own PCA.
- `aadr/examples/PW_JW_arrays_aadr.md` — PW (FTDNA) and JW (MyHeritage) arrays, set 1.
- `aadr_eu/examples/AM_FF_aadr_eu.md` — AM (PGP-HMS FTDNA file), aadr_eu.

In brief:
- On arrays, set 1 rejects only East-Asian/Siberian, Near-Eastern, Ashkenazi and
  pure early-farmer sources. It finds no hidden second component and no separate
  Germany/Corded Ware component. It cannot choose among European sources.
- On WGS with ancient outgroups, PW fits early medieval Poland: Poland_IA
  (= Piast-era 900–1200 CE despite the label), Poland_EarlySlav,
  Poland_EarlyMedieval_Slav, and also Lithuania_LBA. It rejects the Volhynian
  Lusatians, Wielbark, Przeworsk, Trzciniec and Corded Ware as sole sources.
  Regionally, Greater Poland, Masovia and Lubusz fit, while Silesia and Kuyavia do
  not. With pre-medieval sources only, PW needs two: East-Baltic LBA ~30–45 % plus
  Central/Northern European IA.
- PW and JW are indistinguishable on their shared SNPs.
- AM (US, PGP-HMS) is European with no detectable admixture from another continent.
  A continental North-West European early-medieval source (Germany_AltInden) fits
  alone, while British, French IA, Iberian and Slavic bases are rejected alone. That is a lead
  from modern outgroups only.

---

## 7. Runtimes, memory and SNP counts (HO panel, 16 cores, 64 GB RAM)

| Step | Array target | 30x WGS target |
|------|---------|---------|
| 01 / 01b | < 1 min | ~7 min (JOBS=16) |
| 02 | ~1.5 min, 8 GB | same |
| 03–06 | ≤ 2 min each, ≤ 12 GB | a few min each |
| 07 | ~4 min, 36 GB (aadr_eu ~9 min, 52 GB) | similar |
| 08 | ~4 min, 32 GB (aadr_eu ~7 min, 48 GB) | similar |
| 09 | < 1 min | a few min |

- Consumer arrays type ~57K of the 584K HO sites. After `maxmiss=0` across the
  ancient groups, 20–57K SNPs remain, depending on the step. This is expected, not a bug.
- 30x WGS types ~579K HO sites; after `maxmiss=0`, set 2 keeps ~47K (HO) / ~82K
  (1240K). The limit is now the sparse ancient genomes.

---

## 8. Developer reference

- `CLAUDE.md` — pitfalls (all pre-fixed in the scripts), label traps, findings,
  conventions. Shared by both variants.
- `aadr/examples/README_commands.md` — full command history, including the
  array-merge workflow.

## 9. Acknowledgements

Parts of this project were developed with the assistance of Claude Code
(Anthropic) as an AI programming assistant. All methodological decisions,
implementation, testing, and validation were performed by the repository
author.
