# ══════════════════════════════════════════════════════════════════════════════
# OPTIONAL: merge RAW data from multiple platforms before step 01
# Script: aadr/00_merge_arrays.sh <plink_A> <plink_B> <outdir>
# Input:  two pw_plink filesets after step 01 (same person, different platforms)
# Output: <outdir>/pw_plink.{bed,bim,fam} — SNP union, then normal pipeline
#
# Supported CSV formats (step 01):
#   FTDNA Family Finder : RSID,CHROMOSOME,POSITION,RESULT  (no ## header)
#   MyHeritage          : ## reference=build37 + RSID,CHROMOSOME,POSITION,RESULT
#   Ancestry            : rsid,chromosome,position,allele1,allele2  (two allele columns!)
#                         → step 01 auto-detects format from header
#
# Two arrays of the same person → merge at step 00:
#   ./aadr/01_prepare_files.sh /path/to/FTDNA_PW.csv     PERSON_FTDNA
#   ./aadr/01_prepare_files.sh /path/to/MyHeritage_PW.csv PERSON_MH
#   SAMPLE_ID=PW ./aadr/00_merge_arrays.sh \
#       PERSON_FTDNA/pw_plink PERSON_MH/pw_plink PERSON_MERGED
#   # continue pipeline from step 02 with ./PERSON_MERGED as outdir
#
# Three arrays (FTDNA + MyHeritage + Ancestry) — iteratively:
#   SAMPLE_ID=PW ./aadr/00_merge_arrays.sh \
#       PERSON_FTDNA/pw_plink PERSON_MH/pw_plink PERSON_MH_FF
#   SAMPLE_ID=PW ./aadr/00_merge_arrays.sh \
#       PERSON_MH_FF/pw_plink PERSON_ANC/pw_plink PERSON_MERGED
#
# Practical notes (from experience merging MH_PW + PW_FF):
#   - Overlap MH ∩ FTDNA: ~556K / ~600K SNPs (92%) — arrays overlap heavily
#   - Union MH+FTDNA: ~612K SNPs, but after extract_f2(maxmiss=0) FEWER than either alone
#     (~49,218 vs ~50,114 / ~49,883) — edge-of-overlap missing pairs reduce yield
#   - Ancestry DNA (Illumina OmniExpress v2.5, ~730K SNPs) has different content than FTDNA/MH;
#     three-platform union should yield more unique SNPs overlapping AADR
#   - Ancestry format: two allele columns (allele1, allele2) — step 01 must concatenate them
#     into a single RESULT string before conversion; TODO: add support to 01_prepare_files.sh
#   - Allele conflicts (missnp): typically 1–10 SNPs — removed automatically by step 00
#   - "same-position variants" warnings (rs8192284 / rs2228145 etc.) are normal —
#     different rsIDs at the same position; do not block the merge
# ══════════════════════════════════════════════════════════════════════════════

# ── Paths used below ─────────────────────────────────────────────────────────
# Run from the repository root; run directories (./PW_FF, ./JW_MH, ./PW_WGS_HO ...)
# are created there and are git-ignored.
AADR=/data/aadr/v66.p1/geno/v66.p1_HO.aadr.patch.PUB       # CONVERTED GENO prefix (tgeno2geno.py)
AADR1240=/data/aadr/v66.p1/geno/v66.p1_1240K.aadr.patch.PUB
REF=/data/ref/GRCh37.primary_assembly.genome.fa.bgz         # FASTA the CRAM was aligned to

# ── PW (FTDNA Family Finder array) ───────────────────────────────────────────
./aadr/01_prepare_files.sh /path/to/PW_familyfinder.csv PW_FF
Rscript aadr/02_build_aadr_subset.R $AADR ./PW_FF
./aadr/03_merge_and_test.sh ./PW_FF PW
Rscript aadr/04_run_qpadm_test.R      ./PW_FF/merged_pw_aadr_final_ready PW
Rscript aadr/05_robustness_outgroups.R ./PW_FF/merged_pw_aadr_final_ready $AADR PW
Rscript aadr/06_ancestry_models.R     ./PW_FF/merged_pw_aadr_extended    $AADR PW
Rscript aadr/07_extended_models.R     ./PW_FF/merged_pw_aadr_final_ready $AADR PW
Rscript aadr/08_two_source.R          ./PW_FF/merged_pw_aadr_final_ready $AADR PW
Rscript aadr/09_przeworsk_channel.R   ./PW_FF/merged_pw_aadr_final_ready $AADR PW

# ── JW (MyHeritage array) — 5 differences vs PW ─────────────────────────────
# 1. SAMPLE_ID=JW   (env var for step 1; arg 2 for steps 3-9)
# 2. different CSV
# 3. different outdir (JW_MH instead of PW_FF)
# 4. step 3 requires 2nd argument JW  (default "PW" → grep fails in step 3g)
# 5. steps 4-9 use JW as target_id
SAMPLE_ID=JW ./aadr/01_prepare_files.sh /path/to/JW_myheritage.csv JW_MH
Rscript aadr/02_build_aadr_subset.R $AADR ./JW_MH
./aadr/03_merge_and_test.sh ./JW_MH JW
for s in 05_robustness_outgroups 07_extended_models 08_two_source 09_przeworsk_channel; do
  Rscript aadr/$s.R ./JW_MH/merged_pw_aadr_final_ready $AADR JW || break; done
Rscript aadr/06_ancestry_models.R ./JW_MH/merged_pw_aadr_extended $AADR JW
# Results: PW_JW_arrays_aadr.md

# ── PW (30x WGS, CRAM) — HO, then 1240K as a cross-check ─────────────────────
REF=$REF JOBS=16 SAMPLE_ID=PW ./aadr/01b_prepare_wgs.sh /path/to/PW.cram $AADR ./PW_WGS_HO
Rscript aadr/02_build_aadr_subset.R $AADR ./PW_WGS_HO
./aadr/03_merge_and_test.sh ./PW_WGS_HO PW
M=./PW_WGS_HO/merged_pw_aadr_final_ready
# set 1 = steps 04-09 as above, with $M
# set 2 (ancient outgroups, rotating qpAdm, f4 scan):
for s in 10_build_f2 11_outgroup_f3 12_distal 13_rotating 14_f4_scan; do Rscript aadr/set2/$s.R $M $AADR PW || break; done
# set 3 (medieval Polish regions by site):
Rscript aadr/set3/20_build_f2.R $M $AADR PW && Rscript aadr/set3/21_region_tests.R $M $AADR PW
# set 4 (pre-medieval sources only) = set-2 scripts with another population file:
POPS_FILE=$PWD/aadr/set4/pops.R bash -c 'for s in 10_build_f2 11_outgroup_f3 12_distal 13_rotating 14_f4_scan; do Rscript aadr/set2/$s.R '"$M $AADR"' PW || break; done'
# 1240K: the same with $AADR1240 and ./PW_WGS_1240K (step 05 skips outgroups missing from 1240K)
# Results: PW_WGS_aadr.md

# ── PW vs JW on shared SNPs ───────────────────────────────────────────────────
Rscript aadr/pair_compare.R PW_WGS_HO/set3/set3_PW PW ./JW_MH/pw_plink JW $AADR PW_WGS_HO/pair_JW_set3
Rscript aadr/pair_compare.R PW_WGS_HO/set4/set4_PW PW ./JW_MH/pw_plink JW $AADR PW_WGS_HO/pair_JW_set4

# ── Own transparent PCA (25 PCs, West-Eurasian modern basis, lsqproject) ────
python3 -I aadr/pca/30_build_pca.py --prefix $AADR --anno ${AADR%.patch.PUB}.PUB.anno \
    --target PW=./PW_WGS_HO/pw:wgs --target JW=./JW_MH/pw \
    --relabel aadr/pca/relabel_regions.tsv --out ./PW_WGS_HO/pca
# optional: feed an external G25-style nMonte runner (its GroupDef catalogue):
python3 -I aadr/pca/32_to_g25_centroids.py --pca ./PW_WGS_HO/pca --anno ${AADR%.patch.PUB}.PUB.anno \
    --runner /path/to/g25_nmonte_runner.py --out ./PW_WGS_HO/pca/nmonte

# ── AM (public PGP-HMS FamilyTreeDNA file) — aadr_eu variant ────────────────
SAMPLE_ID=AM ./aadr_eu/01_prepare_files.sh /path/to/familytreedna_am_autosomal.csv AM_FF
Rscript aadr_eu/02_build_aadr_subset.R $AADR ./AM_FF
./aadr_eu/03_merge_and_test.sh ./AM_FF AM
Rscript aadr_eu/07_extended_models.R ./AM_FF/merged_pw_aadr_final_ready $AADR AM
Rscript aadr_eu/07b_continental_f3.R ./AM_FF/merged_pw_aadr_final_ready AM
Rscript aadr_eu/08_two_source.R     ./AM_FF/merged_pw_aadr_final_ready $AADR AM
# Results: ../../aadr_eu/examples/AM_FF_aadr_eu.md
