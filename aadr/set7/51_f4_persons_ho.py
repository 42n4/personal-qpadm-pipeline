#!/usr/bin/env python3
"""
51_f4_persons_ho.py -- where single present-day people sit on the Central European - Baltic cline,
on alleles (Human Origins panel, present-day diploid data on every side).

For each individual i:   position_i = f4(Czech, Lithuanian; Czech, i) / f2(Czech, Lithuanian)
  0 = the Czech mean, 1 = the Lithuanian mean. The same number is computed for every present-day
  individual of the neighbouring populations (leave-one-out when i belongs to Czech or Lithuanian), so a
  target is read against the spread of real people on the same SNPs.
Also the plain form f4(Mbuti, Lithuanian; X, target) with a block jackknife.

Targets are PLINK ped/map stems (one person each); an array target drops strand-ambiguous SNPs.
Usage: python3 -I 51_f4_persons_ho.py OUTDIR LABEL=stem[:wgs] ...
"""
import os
import sys

import numpy as np

AADR = "/mnt/db1/aadr/v66.p1"
GENO = f"{AADR}/geno/v66.p1_HO.aadr.patch.PUB.geno"
IND = f"{AADR}/v66.p1_HO.aadr.patch.PUB.ind"
SNP = f"{AADR}/v66.p1_HO.aadr.patch.PUB.snp"
POPS = ["Mbuti", "Czech", "Lithuanian", "Polish", "Belarusian", "Ukrainian", "Ukrainian_North", "Russian",
        "Estonian", "Hungarian", "Croatian"]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_target(stem, wgs, sid_index, a1, a2, nsnp):
    g = np.full(nsnp, 3, np.int8)
    al = open(stem + ".ped").read().split()[6:]
    for k, line in enumerate(open(stem + ".map")):
        j = sid_index.get(line.split()[1])
        if j is None:
            continue
        x, y = al[2 * k], al[2 * k + 1]
        if "0" in (x, y):
            continue
        pair = {a1[j], a2[j]}
        if pair in ({"A", "T"}, {"C", "G"}) and not wgs:
            continue
        if not {x, y} <= pair:
            x, y = COMP.get(x, "?"), COMP.get(y, "?")
            if not {x, y} <= pair:
                continue
        g[j] = (x == a1[j]) + (y == a1[j])
    return g


def jack_ratio(num, den, ok, block):
    ub, inv = np.unique(block[ok], return_inverse=True)
    sn = np.bincount(inv, weights=num[ok]); sd = np.bincount(inv, weights=den[ok])
    tn, td = sn.sum(), sd.sum()
    est = tn / td
    loo = (tn - sn) / (td - sd)
    g = len(ub)
    return est, float(np.sqrt((g - 1) / g * ((loo - loo.mean()) ** 2).sum())), int(ok.sum())


def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    snp = [l.split() for l in open(SNP)]
    nsnp = len(snp)
    sid_index = {s[0]: i for i, s in enumerate(snp)}
    chrom = np.array([int(s[1]) if s[1].isdigit() else 99 for s in snp])
    gpos = np.array([float(s[2]) for s in snp])
    a1 = [s[4] for s in snp]; a2 = [s[5] for s in snp]
    auto = chrom <= 22
    block = chrom * 10000 + np.floor(gpos / 0.05).astype(int)
    ids = [l.split()[0] for l in open(IND)]; lab = [l.split()[2] for l in open(IND)]
    cols = [i for i, l in enumerate(lab) if l in POPS]
    rl = max(48, -(-len(ids) // 4))
    mm = np.memmap(GENO, dtype=np.uint8, mode="r", offset=rl, shape=(nsnp, rl))
    byte = np.array([c // 4 for c in cols]); shift = np.array([6 - 2 * (c % 4) for c in cols], dtype=np.uint8)
    G = np.empty((nsnp, len(cols)), np.uint8)
    for a in range(0, nsnp, 100000):
        G[a:a + 100000] = (np.asarray(mm[a:a + 100000][:, byte]) >> shift) & 3
    plab = [lab[c] for c in cols]; pid = [ids[c] for c in cols]
    people = {f"{plab[k]}:{pid[k]}": G[:, k].astype(np.int8) for k in range(len(cols)) if plab[k] != "Mbuti"}
    targets = {}
    for spec in sys.argv[2:]:
        L, stem = spec.split("=", 1)
        wgs = stem.endswith(":wgs"); stem = stem[:-4] if wgs else stem
        targets[L] = load_target(stem, wgs, sid_index, a1, a2, nsnp)
        print(f"target {L}: {int((targets[L] != 3).sum())} HO SNPs")

    def cnt(pop, drop=None):
        # .DG genomes are left out of every reference frequency: S_Czech-2.DG is the same person as
        # NA15728.HO, and a duplicate inside the reference inflates its sampling variance and shifts
        # every position upward (first run: all Czechs at +0.6 instead of 0)
        k = [i for i in range(len(cols)) if plab[i] == pop and pid[i] != drop and (pop in ("Mbuti", "Polish") or not pid[i].endswith(".DG"))]
        sub = G[:, k]; called = sub != 3
        return np.where(called, sub, 0).sum(1).astype(float), 2.0 * called.sum(1)

    def position(g, drop_id, snpset):
        cz, czn = cnt("Czech", drop_id); lt, ltn = cnt("Lithuanian", drop_id)
        ok = auto & snpset & (g != 3) & (czn >= 10) & (ltn >= 10)
        with np.errstate(invalid="ignore", divide="ignore"):
            pc, pl = cz / czn, lt / ltn
            # E[(pl - pc_hat)(g/2 - pc_hat)] carries + Var(pc_hat); its unbiased estimate is pc(1-pc)/(n-1)
            num = (pl - pc) * (g / 2.0 - pc) - pc * (1 - pc) / (czn - 1)
            den = (pl - pc) ** 2 - pc * (1 - pc) / (czn - 1) - pl * (1 - pl) / (ltn - 1)
        return jack_ratio(np.nan_to_num(num), np.nan_to_num(den), ok, block)

    sets = {"all HO SNPs": np.ones(nsnp, bool)}
    if "JW" in targets:
        sets["JW's SNPs only"] = targets["JW"] != 3
    rows = []
    for sname, sset in sets.items():
        print(f"\n## Position on the Czech (0) - Lithuanian (1) cline, {sname}\n")
        vals = {}
        for name, g in people.items():
            pop, iid = name.split(":")
            if pop in ("Estonian", "Hungarian", "Croatian") or (iid.endswith(".DG") and pop != "Polish"):
                continue
            est, se, n = position(g, iid, sset)
            vals.setdefault(pop, []).append(est)
            rows.append((sname, pop, iid, est, se, n))
        print("| population | n | mean | SD of individuals | min | max |")
        print("|---|---|---|---|---|---|")
        for pop, v in sorted(vals.items(), key=lambda kv: np.mean(kv[1])):
            v = np.array(v)
            print(f"| {pop} | {len(v)} | {v.mean():+.2f} | {v.std(ddof=1) if len(v) > 1 else float('nan'):.2f} | {v.min():+.2f} | {v.max():+.2f} |")
        print("\n| target | position | jackknife SE | SNPs | percentile among Ukrainians + Belarusians |")
        print("|---|---|---|---|---|")
        ref = np.array(vals["Ukrainian"] + vals["Ukrainian_North"] + vals["Belarusian"])
        for L, g in targets.items():
            est, se, n = position(g, None, sset)
            print(f"| {L} | {est:+.2f} | {se:.2f} | {n} | {100 * (ref < est).mean():.0f} |")
            rows.append((sname, "TARGET", L, est, se, n))

    print("\n## f4(Mbuti, Lithuanian; X, target)   positive = the target is closer to Lithuanians than X is\n")
    print("| target | X | n | f4 ×1000 | SE | Z | SNPs |")
    print("|---|---|---|---|---|---|---|")
    mb, mbn = cnt("Mbuti"); lt, ltn = cnt("Lithuanian")
    for L, g in targets.items():
        for X in ("Polish", "Czech", "Ukrainian", "Ukrainian_North", "Belarusian", "Russian"):
            x, xn = cnt(X)
            ok = auto & (g != 3) & (mbn > 0) & (ltn > 0) & (xn > 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                v = np.nan_to_num((mb / mbn - lt / ltn) * (x / xn - g / 2.0))
            est, se, n = jack_ratio(v, np.ones(nsnp), ok, block)
            print(f"| {L} | {X} | {int((np.array(plab) == X).sum())} | {1000 * est:+.3f} | {1000 * se:.3f} | {est / se:+.1f} | {n} |")
    with open(os.path.join(out, "persons_cline.tsv"), "w") as f:
        f.write("snp_set\tpopulation\tid\tposition\tse\tsnps\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
