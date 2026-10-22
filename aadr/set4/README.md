# Set 4: pre-medieval sources only

Same scripts as set 2, run with `POPS_FILE=aadr/set4/pops.R` (output `<run_dir>/set4/`).
No medieval or modern population appears anywhere -- not as a candidate, not in
the right set -- so the question is which Bronze Age / Iron Age / Roman-period
population, or pair of populations, explains the target. AADR v66 has no
Lusatian-culture genomes from Poland and no Zarubintsy, Kyiv culture, Pomeranian
or Milograd samples; the pool is the LBA/IA/Roman ring around Poland (22 groups,
see `pops.R`).

```bash
cd /path/to/personal-qpadm-pipeline; export POPS_FILE=$PWD/aadr/set4/pops.R
M=./PW_WGS_HO/merged_pw_aadr_final_ready; A=$AADR   # converted-GENO prefix, see ../../README.md §2
for s in 10_build_f2 11_outgroup_f3 12_distal 13_rotating 14_f4_scan; do Rscript aadr/set2/$s.R $M $A PW || break; done
```

Read only models that fit on BOTH panels (HO and 1240K) with identifiable weights.
