#!/usr/bin/env bash
# ==============================================================================
# 01b_prepare_wgs.sh -- WGS (CRAM/BAM) -> PLINK at AADR panel positions
#
# WGS counterpart of 01_prepare_files.sh. Instead of parsing a consumer-array
# CSV it genotypes the aligned reads at exactly the SNPs of an AADR panel
# (.snp file), so the output pw_plink.{bed,bim,fam} carries AADR SNP IDs and
# drops straight into steps 02-09 unchanged.
#
# WHY GENOTYPE AT PANEL SITES (not the provider's variant VCF):
#   - a variant-only VCF has no record at hom-ref sites, so "absent" cannot be
#     told apart from "no coverage"; calling from the reads gives an explicit
#     0/0 with depth at every covered panel site.
#   - restricting to panel sites keeps the f-statistic ascertainment of the
#     panel (CLAUDE.md: never use all WGS positions).
#
# REFERENCE: the FASTA must be the one the CRAM was aligned to. A typical
# consumer-WGS CRAM's @RG DS names GRCh37.primary_assembly (GENCODE naming: chr1..chr22,
# chrM = rCRS, GL000191.1 ...). Its @SQ lines carry no M5 checksums, so htslib
# cannot detect a wrong reference -- the script compares contig lengths of
# chr1-22 between the CRAM header and the FASTA index instead. The FASTA may
# be plain or bgzip-compressed; htslib builds .fai/.gzi on first use if the
# directory is writable.
#
# OUTPUT CHOICES:
#   - diploid calls (bcftools mpileup + call -m), like the consumer arrays;
#     AADR ancients are pseudo-haploid, which qpAdm handles per population.
#   - a call is kept only if MIN_DP <= DP <= MAX_DP and GQ >= MIN_GQ, and only
#     if both called alleles are among the panel's two alleles (a third allele
#     is a sequencing/mapping error at a biallelic panel SNP).
#   - bcftools call -m emits no GQ at invariant (0/0, ALT ".") sites; there
#     GQ is taken as 3 x ref-read count (capped 99), so MIN_DP 8 ~ GQ 24.
#   - autosomes 1-22 only (steps 3f/4-9 use autosomes only anyway).
#   - PED/MAP written the same way as step 01 (letters, hom sites as "A A"),
#     so step 03's merge sees the same shape of input as for the arrays.
#
# Usage:
#   REF=/path/GRCh37.primary_assembly.genome.fa.gz \
#     ./aadr/01b_prepare_wgs.sh <sample.cram|bam> <aadr_prefix> [outdir]
#
#   <aadr_prefix>  e.g. $AADR_GENO/v66.p1_1240K.aadr.patch.PUB (converted GENO)
#                  (.snp is read; the SAME prefix must be given to step 02)
#   [outdir]       default: <cram stem>_<panel> under CWD, e.g. ./PW_WGS_1240K
#
# Environment (defaults):
#   REF        required for CRAM; reference FASTA (plain or bgzipped)
#   SAMPLE_ID  PW
#   JOBS       8      chromosomes genotyped in parallel
#   MIN_MQ 20  MIN_BQ 20  MIN_DP 8  MAX_DP 100  MIN_GQ 20
#
# Output:
#   <outdir>/pw_plink.{bed,bim,fam}  -- your sample at AADR panel SNPs, for step 02
#   <outdir>/pw.{map,ped}            -- intermediate PED/MAP
#   <outdir>/wgs_calls/chrN.tsv      -- raw per-site calls (CHROM POS REF ALT GT DP GQ AD)
#   <outdir>/wgs_summary.txt         -- site counts per filter stage
# ==============================================================================
set -euo pipefail

help() { sed -n '3,/^# =====/p' "$0" | sed '$d; s/^# \{0,1\}//'; }
die()  { printf 'ERROR: %b\n' "$*" >&2; exit 1; }

[[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]] && { help; exit 0; }

ALN="$(realpath "${1:?}")"
AADR_PREFIX="${2:?Usage: $0 <sample.cram|bam> <aadr_prefix> [outdir]}"
SAMPLE_ID="${SAMPLE_ID:-PW}"
JOBS="${JOBS:-8}"
MIN_MQ="${MIN_MQ:-20}"; MIN_BQ="${MIN_BQ:-20}"
MIN_DP="${MIN_DP:-8}";  MAX_DP="${MAX_DP:-100}"; MIN_GQ="${MIN_GQ:-20}"
REF="${REF:-}"

for t in bcftools htsfile plink2 python3 xargs; do
  command -v "$t" >/dev/null || die "$t not found in PATH"
done
[[ -f "$ALN" ]] || die "alignment not found: $ALN"
[[ -f "${AADR_PREFIX}.snp" ]] || die "panel not found: ${AADR_PREFIX}.snp"
case "$ALN" in
  *.cram) [[ -f "$ALN.crai" || -f "${ALN%.cram}.crai" ]] || die "CRAM index (.crai) missing"
          [[ -n "$REF" && -f "$REF" ]] || die "CRAM needs REF=<reference fasta>" ;;
  *.bam)  [[ -f "$ALN.bai" || -f "${ALN%.bam}.bai" ]] || die "BAM index (.bai) missing" ;;
  *)      die "expected a .cram or .bam file: $ALN" ;;
esac

panel="$(basename "$AADR_PREFIX")"; panel="${panel#*_}"; panel="${panel%%.*}"   # HO / 1240K / 2M
stem="$(basename "$ALN")"; stem="${stem%%.*}"
OUTDIR="$(realpath -m "${3:-$(pwd)/${stem}_${panel}}")"
mkdir -p "$OUTDIR/wgs_calls"
cd "$OUTDIR"

echo "================================================================================"
echo "  STEP 1b: WGS reads -> PLINK at AADR ${panel} panel sites"
echo "================================================================================"
echo "  Alignment : $ALN"
echo "  Panel     : ${AADR_PREFIX}.snp"
echo "  Reference : ${REF:-<from BAM, none given>}"
echo "  Output dir: $OUTDIR"
echo "  Filters   : MQ>=$MIN_MQ BQ>=$MIN_BQ  $MIN_DP<=DP<=$MAX_DP  GQ>=$MIN_GQ   jobs=$JOBS"
echo

# ---- Contig naming + reference sanity check --------------------------------
# No M5 tags in the CRAM header means htslib will decode against ANY fasta
# with matching names -- a wrong build silently yields wrong bases. Compare
# chr1-22 lengths header vs fasta before spending an hour on pileups.
htsfile -h "$ALN" | awk -F'\t' '$1=="@SQ"{n="";l="";for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)}print n"\t"l}' > aln_contigs.tsv
[[ -s aln_contigs.tsv ]] || die "could not read @SQ lines from $ALN"
if grep -q $'^chr1\t' aln_contigs.tsv; then CHR_PREFIX=chr
elif grep -q $'^1\t' aln_contigs.tsv; then CHR_PREFIX=
else die "neither chr1 nor 1 in the alignment header"; fi

if [[ -n "$REF" ]]; then
  if [[ ! -f "$REF.fai" ]]; then
    echo "  Indexing reference (one-off; needs a bgzip or plain FASTA)..."
    # A plain-gzip FASTA cannot be indexed -- htslib needs bgzip.
    if [[ "$REF" == *.gz ]] && ! htsfile "$REF" | grep -q BGZF; then
      die "REF is plain gzip; recompress: zcat '$REF' | bgzip -@8 > '${REF%.gz}.bgz.gz'"
    fi
    bcftools mpileup -f "$REF" -r "${CHR_PREFIX}22:1-2" "$ALN" >/dev/null 2>&1 \
      || die "could not index/read REF=$REF"
  fi
  [[ -f "$REF.fai" ]] || die "no $REF.fai after indexing attempt"
  bad=$(awk -F'\t' -v p="$CHR_PREFIX" 'NR==FNR{len[$1]=$2;next}
        {c=$1; for(i=1;i<=22;i++) if(c==p i){ if(!(c in len)) print c" missing in fasta";
          else if(len[c]!=$2) print c" length "$2" (aln) vs "len[c]" (fasta)" }}' \
        "$REF.fai" aln_contigs.tsv)
  [[ -z "$bad" ]] || die "reference does not match the alignment:\n$bad"
  echo "  Reference matches the alignment header on chr1-22 (names + lengths)."
fi

# ---- Panel sites -----------------------------------------------------------
# EIGENSTRAT .snp: ID chrom genpos physpos allele1 allele2; chrom 23/24/90 =
# X/Y/MT. Keep autosomal biallelic ACGT SNPs; one row per position (the panel
# has a handful of duplicate positions under different IDs -- keep the first).
echo "[1b-1] Extracting autosomal panel sites..."
python3 -I - "${AADR_PREFIX}.snp" "$CHR_PREFIX" <<'PYEOF'
import sys
snp, pre = sys.argv[1], sys.argv[2]
seen = set(); n = kept = dup = skip = 0
with open(snp) as f, open("panel_sites.tsv", "w") as out:
    for line in f:
        p = line.split()
        if len(p) < 6: continue
        n += 1
        sid, chrom, _, pos, a1, a2 = p[:6]
        if chrom not in {str(i) for i in range(1, 23)} or a1 not in "ACGT" or a2 not in "ACGT" \
           or len(a1) != 1 or len(a2) != 1 or a1 == a2:
            skip += 1; continue
        key = (chrom, pos)
        if key in seen: dup += 1; continue
        seen.add(key); kept += 1
        out.write(f"{pre}{chrom}\t{pos}\t{sid}\t{a1}\t{a2}\n")
print(f"  panel SNPs {n}, autosomal biallelic kept {kept}, skipped {skip}, duplicate positions {dup}")
PYEOF
sort -k1,1V -k2,2n panel_sites.tsv -o panel_sites.tsv
cut -f1,2 panel_sites.tsv > panel_targets.tsv

# ---- Per-chromosome genotyping (parallel) -----------------------------------
# -T (targets, streaming) not -R (regions, index jump per site): with ~10^6
# sites the streaming pass over each chromosome is far cheaper. -r confines
# each job to one chromosome so the jobs split the CRAM between them.
echo
echo "[1b-2] Genotyping panel sites, $JOBS chromosomes at a time (this is the long step)..."
export ALN REF MIN_MQ MIN_BQ CHR_PREFIX
call_chr() {
  local c="$1" out="wgs_calls/$1.tsv"
  [[ -f "$out.ok" ]] && { echo "  $c cached"; return 0; }
  rm -f "$out" "$out.ok"
  # Pipeline status: any stage failing must fail the job (pipefail per job).
  set -o pipefail
  bcftools mpileup ${REF:+-f "$REF"} -r "$c" -T panel_targets.tsv \
      -q "$MIN_MQ" -Q "$MIN_BQ" -a FORMAT/AD,FORMAT/DP --max-depth 500 -Ou "$ALN" 2>"wgs_calls/$c.mpileup.log" \
    | bcftools call -m -a GQ -Ou 2>"wgs_calls/$c.call.log" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t%AD]\n' 2>"wgs_calls/$c.query.log" > "$out" \
    && touch "$out.ok" && echo "  $c done ($(wc -l < "$out") sites)"
}
export -f call_chr
printf "${CHR_PREFIX}%s\n" $(seq 1 22) \
  | xargs -P "$JOBS" -I{} bash -c 'call_chr "$@"' _ {} || true
missing=()
for i in $(seq 1 22); do [[ -f "wgs_calls/${CHR_PREFIX}$i.tsv.ok" ]] || missing+=("${CHR_PREFIX}$i"); done
(( ${#missing[@]} == 0 )) || die "genotyping failed for: ${missing[*]} (see wgs_calls/<chr>.*.log); re-run to retry only those"

# ---- Filter + PED/MAP --------------------------------------------------------
echo
echo "[1b-3] Filtering calls and writing PED/MAP..."
python3 -I - "$SAMPLE_ID" "$MIN_DP" "$MAX_DP" "$MIN_GQ" <<'PYEOF'
import sys, glob
sample, min_dp, max_dp, min_gq = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
panel = {}
order = []
with open("panel_sites.tsv") as f:
    for line in f:
        c, p, sid, a1, a2 = line.rstrip("\n").split("\t")
        panel[(c, p)] = (sid, a1, a2); order.append((c, p))
st = dict(panel=len(panel), covered=0, low_dp=0, high_dp=0, low_gq=0, nocall=0,
          third_allele=0, dup_record=0, kept=0, hom_ref=0, het=0, hom_alt=0)
calls = {}
for fn in glob.glob("wgs_calls/*.tsv"):
    with open(fn) as f:
        for line in f:
            c, p, ref, alt, gt, dp, gq, ad = line.rstrip("\n").split("\t")
            key = (c, p)
            if key not in panel: continue
            # bcftools occasionally emits two SNV records for one position;
            # keep the first so counts match the written .map.
            if key in calls: st["dup_record"] += 1; continue
            # mpileup can add an indel record at a SNP position; only the SNV
            # record describes the panel site.
            if len(ref) != 1 or any(len(x) != 1 for x in alt.split(",") if x != "."): continue
            st["covered"] += 1
            dp = int(dp) if dp not in (".", "") else 0
            if gq not in (".", ""):
                gq = int(float(gq))
            elif gt == "0/0":
                # call -m leaves GQ empty at invariant sites (ALT "."). Each
                # ref read halves the odds of a missed het: GQ ~ 10*log10(2^n).
                nref = int(ad.split(",")[0]) if ad not in (".", "") else 0
                gq = min(99, round(3.0103 * nref))
            else:
                gq = 0
            if dp < min_dp: st["low_dp"] += 1; continue
            if dp > max_dp: st["high_dp"] += 1; continue
            if gq < min_gq: st["low_gq"] += 1; continue
            idx = gt.replace("|", "/").split("/")
            if len(idx) != 2 or "." in idx: st["nocall"] += 1; continue
            alleles = [ref] + (alt.split(",") if alt != "." else [])
            try: g = [alleles[int(i)] for i in idx]
            except (IndexError, ValueError): st["nocall"] += 1; continue
            _, a1, a2 = panel[key]
            if any(x not in (a1, a2) for x in g): st["third_allele"] += 1; continue
            if idx == ["0", "0"]: st["hom_ref"] += 1
            elif idx[0] != idx[1]: st["het"] += 1
            else: st["hom_alt"] += 1
            calls[key] = g; st["kept"] += 1
with open("pw.map", "w") as fm:
    for key in order:
        if key in calls:
            fm.write(f"{key[0].removeprefix('chr')}\t{panel[key][0]}\t0\t{key[1]}\n")
with open("pw.ped", "w") as fp:
    fp.write(f"{sample}\t{sample}\t0\t0\t0\t-9\t" +
             "\t".join(f"{calls[k][0]}\t{calls[k][1]}" for k in order if k in calls) + "\n")
with open("wgs_summary.txt", "w") as fs:
    for k, v in st.items(): fs.write(f"{k}\t{v}\n")
for k, v in st.items(): print(f"  {k:<13}{v:>10}")
if st["kept"] == 0:
    sys.exit("no calls survived the filters")
het_rate = st["het"] / st["kept"]
print(f"  het fraction among kept calls: {het_rate:.3f}")
PYEOF

# ---- PED/MAP -> binary PLINK (same as step 01) --------------------------------
echo
plink2 --pedmap pw --make-bed --out pw_plink > pw_plink.convert.log 2>&1 \
  || { tail -20 pw_plink.convert.log >&2; die "plink2 --pedmap failed"; }
# plink2 writes "." as A1 where the sample is homozygous (~75% of WGS panel
# sites). p-link --bmerge can flag "." vs a real allele as "3+ alleles"
# (CLAUDE.md pitfall 9) and step 3d would then drop those SNPs. Fill the
# unobserved allele from the panel pair: no genotype uses that code, so the
# .bed is unchanged.
python3 -I - <<'PYEOF'
pair = {}
with open("panel_sites.tsv") as f:
    for line in f:
        _, _, sid, a1, a2 = line.rstrip("\n").split("\t")
        pair[sid] = (a1, a2)
out = []; filled = 0
with open("pw_plink.bim") as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if c[4] in (".", "0"):
            a1, a2 = pair[c[1]]
            c[4] = a2 if c[5] == a1 else a1
            filled += 1
        if {c[4], c[5]} != set(pair[c[1]]):
            raise SystemExit(f"allele mismatch vs panel at {c[1]}: {c[4]}/{c[5]}")
        out.append("\t".join(c))
with open("pw_plink.bim", "w") as f:
    f.write("\n".join(out) + "\n")
print(f"  Filled the unobserved allele at {filled} homozygous sites; all {len(out)} match the panel pair.")
PYEOF

plink2 --bfile pw_plink --missing --out pw_plink_selfcheck > /dev/null 2>&1 \
  || die "plink2 --missing failed"
echo "  Self-check (F_MISS should be 0 -- only kept calls were written):"
cat pw_plink_selfcheck.smiss

echo
echo "================================================================================"
echo "  STEP 1b DONE. Output: $OUTDIR/pw_plink.{bed,bim,fam}"
echo "  Next (same AADR prefix!):"
echo "    Rscript aadr/02_build_aadr_subset.R $AADR_PREFIX $OUTDIR"
echo "    ./aadr/03_merge_and_test.sh $OUTDIR $SAMPLE_ID"
echo "================================================================================"
