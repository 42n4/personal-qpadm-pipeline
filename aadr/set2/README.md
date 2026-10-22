# Set 2: comparisons that can tell European sources apart

Set 1 (steps 04-09) tests the target against candidates with only Mbuti, Yoruba
and Han as outgroups. Those cannot separate one European source from another:
on PW WGS x HO, 34 of 38 single sources fit, and the p-values do not rank them.
Set 2 runs on the same step-03 merge but asks different questions.

| Script | Question | Needs a right set? |
|---|---|---|
| `10_build_f2.R` | adds the populations in `pops.R`, builds one f2 cache | - |
| `11_outgroup_f3.R` | which population shares the most drift with the target | no |
| `12_distal.R` | WHG / EHG or steppe / EEF shares of the target vs every candidate | ancient deep outgroups |
| `13_rotating.R` | rotating qpAdm: unused candidates become outgroups | ancient + rotating pool |
| `14_f4_scan.R` | which populations pull the target away from its best sources | no |

All population lists live in `pops.R`; edit them there.

## Run

```bash
cd /path/to/personal-qpadm-pipeline   # run dirs are created here
M=./PW_WGS_HO/merged_pw_aadr_final_ready
A=$AADR   # converted-GENO prefix, see ../../README.md §2
for s in 10_build_f2 11_outgroup_f3 12_distal 13_rotating 14_f4_scan; do
  Rscript aadr/set2/$s.R $M $A PW || break
done
```

Outputs land in `<run_dir>/set2/`. Script 10 is idempotent on the fileset and
always rebuilds the f2 cache. `MAXMISS` (script 10), `MAX_SRC` (13) and
`TOP` / `SOURCES` (14) are environment overrides.

## Reading the results

- **11**: treat populations within 2 SE of the top (`gap_to_top_se < 2`) as tied.
- **12**: `|dz| >= 3` on any component rules a population out as a sole source;
  `|dz| < 2` everywhere means distally indistinguishable from the target.
- **13**: only `FIT` rows matter. `UNIDENTIFIABLE` = fits, but a weight SE exceeds
  `MAX_SE` (0.25): collinear sources, any weights fit, not evidence.
  A 2-source `FIT` is evidence of admixture only when both sources are rejected
  on their own (`note` column).
- **14**: `z >= 3` means the target has ancestry related to X that the source
  lacks; `|z| < 3` everywhere means source and target form a clade.

Works on HO and 1240K: every reference label is present in both panels
(checked 2026-10-21). Candidates missing from a panel are skipped with a warning.
