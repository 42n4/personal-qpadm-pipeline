#!/usr/bin/env python3
"""
50_f4_sites_baltic.py -- f4(Outgroup, Balts; Saale Slavs, X) for medieval Polish sites.

The allele-level twin of the G25 table "site mean = Saale Slavs + Iron Age Balts": a positive f4
means X shares more drift with the Baltic group than the Saale Slavs do. One number with a block
jackknife Z per site; individuals with few SNPs still contribute, because frequencies are counted
per SNP over whoever is called.

Inputs (no merge is needed: the two packages have the same 1,233,013 SNPs in the same order):
  AADR v66 1240K, decoded GENO (packedancestrymap) + .ind + .snp + .anno
  Poseidon package 2025_Gretzinger_SlavicPeriod (PLINK bed/bim/fam + .janno)
Output: a table on stdout and OUT/f4_sites.tsv.

Usage: python3 -I 50_f4_sites_baltic.py OUTDIR [--min-snps 30000]
"""
import argparse
import csv
import os
import sys

import numpy as np

AADR = "/mnt/db1/aadr/v66.p1"
GENO = f"{AADR}/geno/v66.p1_1240K.aadr.patch.PUB.geno"
IND = f"{AADR}/v66.p1_1240K.aadr.patch.PUB.ind"
SNP = f"{AADR}/v66.p1_1240K.aadr.patch.PUB.snp"
ANNO = f"{AADR}/v66.p1_1240K.aadr.PUB.anno"
GRZ = "/mnt/db5/aadr_extra/gretzinger_slavic/2025_Gretzinger_SlavicPeriod"

SITES = [  # label, locality substring (AADR anno), all Poland, 900-1400 CE
    ("Santok", "Santok"), ("Niemcza", "Niemcza"), ("Ląd", "Ląd ("), ("Milicz", "Milicz"),
    ("Poznań-Śródka", "Poznań-Śródka"), ("Dziekanowice", "Dziekanowice"), ("Markowice", "Markowice"),
    ("Końskie", "Końskie"), ("Ostrów Lednicki", "Ostrów Lednicki"), ("Płońsk", "Płońsk"),
    ("Łęgowo", "Łęgowo"), ("Obłaczkowo", "Oblaczkowo"),
]


def read_aadr_meta(min_snps):
    ids = [l.split()[0] for l in open(IND)]
    lab = [l.split()[2] for l in open(IND)]
    pos = {g: i for i, g in enumerate(ids)}
    groups = {"Mbuti": [i for i, l in enumerate(lab) if l == "Mbuti"]}
    balt_lt, balt_lv, german = [], [], []
    for k, _ in SITES:
        groups[k] = []
    with open(ANNO, encoding="utf-8", errors="replace") as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        for p in r:
            if p[0] not in pos:
                continue
            try:
                bp, n = float(p[10]), int(p[26])
            except ValueError:
                continue
            if p[16] == "Poland" and 550 <= bp <= 1100 and n >= min_snps and "-o" not in p[14]:
                for k, sub in SITES:
                    if sub in p[15]:
                        groups[k].append(pos[p[0]])
            if p[16] == "Germany" and 850 <= bp <= 1550 and n >= 100000 and "-o" not in p[14]:
                german.append(pos[p[0]])
            if 1350 <= bp <= 2750 and n >= 100000:
                if p[16] == "Lithuania":
                    balt_lt.append(pos[p[0]])
                elif p[16] == "Latvia":
                    balt_lv.append(pos[p[0]])
    groups["Balts LT+LV 800 BCE-600 CE"] = balt_lt + balt_lv
    groups["Balts LT only (1st-6th c.)"] = balt_lt
    groups["Balts LV only (Kivutkalns, capture)"] = balt_lv
    groups["Germany 400-1100 CE (AADR)"] = german
    # deeper individuals only, to test coverage as a confounder on alleles
    deep = {}
    with open(ANNO, encoding="utf-8", errors="replace") as f:
        r = csv.reader(f, delimiter="\t"); next(r)
        for p in r:
            try:
                deep[p[0]] = int(p[26])
            except (ValueError, IndexError):
                pass
    for k, _ in SITES:
        groups[k + " >=200k"] = [i for i in groups[k] if deep.get(ids[i], 0) >= 200000]
    return len(ids), groups


def aadr_counts(nind, groups, nsnp):
    """ref-allele count and called-chromosome count per group and SNP."""
    rl = max(48, -(-nind // 4))
    mm = np.memmap(GENO, dtype=np.uint8, mode="r", offset=rl, shape=(nsnp, rl))
    cols = sorted({i for g in groups.values() for i in g})
    cidx = {c: k for k, c in enumerate(cols)}
    byte = np.array([c // 4 for c in cols])
    shift = np.array([6 - 2 * (c % 4) for c in cols], dtype=np.uint8)
    g = np.empty((nsnp, len(cols)), dtype=np.uint8)
    for a in range(0, nsnp, 100000):
        g[a:a + 100000] = (np.asarray(mm[a:a + 100000][:, byte]) >> shift) & 3
    out = {}
    for k, idx in groups.items():
        sub = g[:, [cidx[i] for i in idx]]
        called = sub != 3
        out[k] = (np.where(called, sub, 0).sum(1).astype(np.float64), 2.0 * called.sum(1))
    per_ind = {c: g[:, cidx[c]] for c in cols}
    return out, per_ind


def plink_counts(groups_by_iid, ref):
    fam = [l.split() for l in open(GRZ + ".fam")]
    n = len(fam)
    bim = [l.split() for l in open(GRZ + ".bim")]
    nsnp = len(bim)
    rl = -(-n // 4)
    mm = np.memmap(GRZ + ".bed", dtype=np.uint8, mode="r", offset=3, shape=(nsnp, rl))
    a1 = np.array([b[4] for b in bim]); a2 = np.array([b[5] for b in bim])
    same, flip = a1 == ref, a2 == ref
    iid = {f[1]: i for i, f in enumerate(fam)}
    out = {}
    for k, ids in groups_by_iid.items():
        cols = [iid[x] for x in ids if x in iid]
        byte = np.array([c // 4 for c in cols]); shift = np.array([2 * (c % 4) for c in cols], dtype=np.uint8)
        cnt = np.zeros(nsnp); tot = np.zeros(nsnp)
        for a in range(0, nsnp, 100000):
            code = (np.asarray(mm[a:a + 100000][:, byte]) >> shift) & 3     # 0 hom A1, 1 missing, 2 het, 3 hom A2
            called = code != 1
            a1c = np.where(code == 0, 2, np.where(code == 2, 1, 0))
            cnt[a:a + 100000] = np.where(called, a1c, 0).sum(1); tot[a:a + 100000] = 2.0 * called.sum(1)
        refc = np.where(same, cnt, np.where(flip, tot - cnt, 0.0))
        tot = np.where(same | flip, tot, 0.0)
        out[k] = (refc, tot, len(cols))
    return out, int((~(same | flip)).sum())


def f4_jack(pO, pB, pS, pX, ok, block):
    v = (pO - pB) * (pS - pX)
    v, b = v[ok], block[ok]
    ub, inv = np.unique(b, return_inverse=True)
    s = np.bincount(inv, weights=v); m = np.bincount(inv).astype(float)
    n, tot = m.sum(), s.sum()
    est = tot / n
    loo = (tot - s) / (n - m)                       # leave-one-block-out estimates
    g = len(ub)
    jk = g * est - ((n - m) / n * loo).sum()        # weighted jackknife (Busing et al.)
    h = n / m
    tau = h * est - (h - 1) * loo
    var = ((tau - jk) ** 2 / (h - 1)).sum() / g
    return est, np.sqrt(var), int(n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--min-snps", type=int, default=30000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    snp = [l.split() for l in open(SNP)]
    nsnp = len(snp)
    chrom = np.array([int(s[1]) if s[1].isdigit() else 99 for s in snp])
    gpos = np.array([float(s[2]) for s in snp])
    ref = np.array([s[4] for s in snp]); alt = np.array([s[5] for s in snp])
    auto = chrom <= 22
    tv = ~(((ref == "C") & (alt == "T")) | ((ref == "T") & (alt == "C")) | ((ref == "G") & (alt == "A")) | ((ref == "A") & (alt == "G")))
    block = chrom * 10000 + np.floor(gpos / 0.05).astype(int)      # 5 cM blocks

    nind, groups = read_aadr_meta(a.min_snps)
    A, per_ind = aadr_counts(nind, groups, nsnp)

    jan = list(csv.reader(open(GRZ + ".janno", encoding="utf-8"), delimiter="\t"))
    h = {c: i for i, c in enumerate(jan[0])}
    rng = np.random.default_rng(7)
    saale = [r[h["Poseidon_ID"]] for r in jan[1:] if r[h["Group_Name"]] == "EGermany_SP" and int(r[h["Nr_SNPs"]] or 0) >= 100000]
    perm = rng.permutation(len(saale))
    gz = {"Saale (EGermany_SP)": saale,
          "Saale half A": [saale[i] for i in perm[: len(saale) // 2]],
          "Saale half B": [saale[i] for i in perm[len(saale) // 2:]],
          "Gródek (Poland_SP)": [r[h["Poseidon_ID"]] for r in jan[1:] if r[h["Group_Name"]] == "Poland_SP" and int(r[h["Nr_SNPs"]] or 0) >= 100000],
          "Volhynia (NWUkraine_SP)": [r[h["Poseidon_ID"]] for r in jan[1:] if r[h["Group_Name"]] == "NWUkraine_SP" and int(r[h["Nr_SNPs"]] or 0) >= 100000],
          "Saale-region Migration period (EGermany_MP)": [r[h["Poseidon_ID"]] for r in jan[1:] if r[h["Group_Name"]] == "EGermany_MP" and int(r[h["Nr_SNPs"]] or 0) >= 100000]}
    P, bad = plink_counts(gz, ref)
    print(f"SNPs {nsnp}, autosomal {int(auto.sum())}, transversions {int((auto & tv).sum())}; allele mismatch between packages: {bad}")

    def freq(c, t):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(t > 0, c / t, np.nan)
    F = {k: freq(*v) for k, v in A.items()}
    F.update({k: freq(v[0], v[1]) for k, v in P.items()})
    N = {k: len(v) for k, v in groups.items()}
    N.update({k: v[2] for k, v in P.items()})

    rows = []
    for balt in ("Balts LT+LV 800 BCE-600 CE", "Balts LT only (1st-6th c.)"):
        print(f"\n## f4(Mbuti n={N['Mbuti']}, {balt} n={N[balt]}; Saale n={N['Saale (EGermany_SP)']}, X)   positive = X is closer to the Balts than the Saale Slavs are\n")
        print("| X | n | f4 ×1000 | SE | Z | SNPs | Z, transversions only | individuals with f4 > 0 |")
        print("|---|---|---|---|---|---|---|---|")
        tests = [(k, "Saale (EGermany_SP)") for k, _ in SITES] + [("Gródek (Poland_SP)", "Saale (EGermany_SP)"),
                 ("Volhynia (NWUkraine_SP)", "Saale (EGermany_SP)"), ("Saale-region Migration period (EGermany_MP)", "Saale (EGermany_SP)"),
                 ("Saale half B", "Saale half A")]
        for x, s in tests:
            if N[x] < 2:
                continue
            pO, pB, pS, pX = F["Mbuti"], F[balt], F[s], F[x]
            ok = auto & ~np.isnan(pO) & ~np.isnan(pB) & ~np.isnan(pS) & ~np.isnan(pX)
            est, se, n = f4_jack(pO, pB, pS, pX, ok, block)
            est2, se2, _n2 = f4_jack(pO, pB, pS, pX, ok & tv, block)
            pos = ""
            if x in groups:                                     # per individual, point estimate only
                k = 0
                for i in groups[x]:
                    gi = per_ind[i]
                    oki = ok & (gi != 3)
                    k += float(((pO - pB) * (pS - gi / 2.0))[oki].mean()) > 0
                pos = f"{k} of {len(groups[x])}"
            lab = x if s.startswith("Saale (") else f"{x} (against half A: the null)"
            print(f"| {lab} | {N[x]} | {1000 * est:+.3f} | {1000 * se:.3f} | {est / se:+.1f} | {n} | {est2 / se2:+.1f} | {pos} |")
            rows.append((balt, lab, N[x], est, se, est / se, n, est2 / se2, pos))
    # ---- same-technology variant: both populations of the right pair are Stolarek 2023 shotgun data.
    # The Saale Slavs are 1240K capture and the Polish sites shotgun; capture and shotgun data attract
    # their own kind in f-statistics, so the block above mixes ancestry with data type.
    core = ("Ostrów Lednicki", "Płońsk", "Łęgowo", "Obłaczkowo")
    rest = ("Santok", "Niemcza", "Ląd")

    def pooled(keys):
        c = sum(A[k][0] for k in keys); t = sum(A[k][1] for k in keys)
        return freq(c, t), sum(N[k] for k in keys)
    F["Piast core, pooled"], N["Piast core, pooled"] = pooled(core)
    for balt in ("Balts LT only (1st-6th c.)", "Balts LT+LV 800 BCE-600 CE"):
        print(f"\n## Same data type: f4(Mbuti, {balt} n={N[balt]}; Santok+Niemcza+Ląd, X), all Polish sites shotgun (Stolarek 2023)\n")
        print("| X | n | reference pool n | f4 ×1000 | SE | Z | SNPs | Z, transversions only |")
        print("|---|---|---|---|---|---|---|---|")
        for x in ("Piast core, pooled",) + core + rest + ("Milicz", "Poznań-Śródka", "Dziekanowice", "Markowice", "Końskie"):
            if N[x] < 2:
                continue
            pS, nS = pooled([k for k in rest if k != x])
            pO, pB, pX = F["Mbuti"], F[balt], F[x]
            ok = auto & ~np.isnan(pO) & ~np.isnan(pB) & ~np.isnan(pS) & ~np.isnan(pX)
            est, se, n = f4_jack(pO, pB, pS, pX, ok, block)
            est2, se2, _n2 = f4_jack(pO, pB, pS, pX, ok & tv, block)
            print(f"| {x} | {N[x]} | {nS} | {1000 * est:+.3f} | {1000 * se:.3f} | {est / se:+.1f} | {n} | {est2 / se2:+.1f} |")
            rows.append((balt + " / same data type", x, N[x], est, se, est / se, n, est2 / se2, ""))
    # ---- checks of 2026-11-22: (1) reference pool without Niemcza, (2) mirror test with a Central European
    # group in place of the Balts, (3) deeper individuals only on both sides
    def pooled2(keys, suf=""):
        c = sum(A[k + suf][0] for k in keys); t = sum(A[k + suf][1] for k in keys)
        return freq(c, t), sum(N[k + suf] for k in keys)
    seconds = [("Balts LT only (1st-6th c.)", "Balts, Lithuania (shotgun)"), ("Balts LT+LV 800 BCE-600 CE", "Balts, Lithuania + Latvia"),
               ("Saale-region Migration period (EGermany_MP)", "Central Europe: Saale region, Migration period (capture)"),
               ("Germany 400-1100 CE (AADR)", "Central Europe: Germany 400-1100 CE (AADR)")]
    print("\n## Checks: f4(Mbuti, P2; reference, Piast core pooled)   positive = Piast core is closer to P2 than the reference is\n")
    print("| P2 | n | reference | ref n | core n | f4 ×1000 | SE | Z | Z, transversions |")
    print("|---|---|---|---|---|---|---|---|---|")
    for p2, lab in seconds:
        for refk, suf, rl in ((rest, "", "Santok + Niemcza + Ląd"), (("Santok", "Ląd"), "", "Santok + Ląd"), (("Niemcza",), "", "Niemcza alone"),
                              (rest, " >=200k", "Santok + Niemcza + Ląd, >= 200k SNPs both sides"), (("Santok", "Ląd"), " >=200k", "Santok + Ląd, >= 200k SNPs both sides")):
            pS, nS = pooled2(refk, suf); pX, nX = pooled2(core, suf)
            pO, pB = F["Mbuti"], F[p2]
            ok = auto & ~np.isnan(pO) & ~np.isnan(pB) & ~np.isnan(pS) & ~np.isnan(pX)
            est, se, n = f4_jack(pO, pB, pS, pX, ok, block)
            est2, se2, _ = f4_jack(pO, pB, pS, pX, ok & tv, block)
            print(f"| {lab} | {N[p2]} | {rl} | {nS} | {nX} | {1000 * est:+.3f} | {1000 * se:.3f} | {est / se:+.1f} | {est2 / se2:+.1f} |")
            rows.append((lab + " / check", rl, nX, est, se, est / se, n, est2 / se2, ""))
    # ---- cline statistic without an outgroup: f4(Central Europe, Balts; reference, X).
    # With Mbuti in the first place a positive value can also mean "X is simply better data" (closer to every
    # West Eurasian group); the mirror test showed that. Here the first pair spans the cline itself.
    print("\n## Cline: f4(Central Europe, Balts; reference, X)   positive = X lies further toward the Balts on the Baltic - Central European cline\n")
    print("| Central Europe | Balts | reference | X | n | f4 ×1000 | SE | Z | Z, transversions |")
    print("|---|---|---|---|---|---|---|---|---|")
    for ce in ("Germany 400-1100 CE (AADR)", "Saale-region Migration period (EGermany_MP)"):
        for balt in ("Balts LT only (1st-6th c.)", "Balts LT+LV 800 BCE-600 CE"):
            for refk, rl in ((rest, "Santok + Niemcza + Ląd"), (("Santok", "Ląd"), "Santok + Ląd")):
                xs = [("Piast core, pooled", None)] + [(k, k) for k in core] + [(k, k) for k in ("Niemcza", "Dziekanowice", "Końskie")]
                for x, key in xs:
                    pS, nS = pooled2([k for k in refk if k != x])
                    pX = F[x]
                    ok = auto & ~np.isnan(F[ce]) & ~np.isnan(F[balt]) & ~np.isnan(pS) & ~np.isnan(pX)
                    est, se, n = f4_jack(F[ce], F[balt], pS, pX, ok, block)
                    est2, se2, _ = f4_jack(F[ce], F[balt], pS, pX, ok & tv, block)
                    print(f"| {ce.split(' (')[0]} | {balt.split(' ')[1]} {balt.split(' ')[2] if '+' in balt else ''} | {rl} | {x} | {N[x]} | {1000 * est:+.3f} | {1000 * se:.3f} | {est / se:+.1f} | {est2 / se2:+.1f} |")
                    rows.append((f"cline {ce} / {balt}", f"{rl} -> {x}", N[x], est, se, est / se, n, est2 / se2, ""))
    # ---- estimate and error side by side, all SNPs against transversions, three Baltic groups
    print("\n## Cline, Piast core pooled against Santok + Ląd: estimate and error separately\n")
    print("| Central Europe | Balts | n Balts | f4 ×1000, all SNPs | SE | Z | f4 ×1000, transversions | SE | Z |")
    print("|---|---|---|---|---|---|---|---|---|")
    pS, nS = pooled2(("Santok", "Ląd")); pX = F["Piast core, pooled"]
    for ce in ("Germany 400-1100 CE (AADR)", "Saale-region Migration period (EGermany_MP)"):
        for balt in ("Balts LT only (1st-6th c.)", "Balts LV only (Kivutkalns, capture)", "Balts LT+LV 800 BCE-600 CE"):
            ok = auto & ~np.isnan(F[ce]) & ~np.isnan(F[balt]) & ~np.isnan(pS) & ~np.isnan(pX)
            e1, s1, _ = f4_jack(F[ce], F[balt], pS, pX, ok, block)
            e2, s2, _ = f4_jack(F[ce], F[balt], pS, pX, ok & tv, block)
            print(f"| {ce.split(' (')[0]} | {balt} | {N[balt]} | {1000 * e1:+.3f} | {1000 * s1:.3f} | {e1 / s1:+.1f} | {1000 * e2:+.3f} | {1000 * s2:.3f} | {e2 / s2:+.1f} |")
    # ---- length of the cline for each Baltic group: f4 is not normalised, so a more distant second
    # population gives a larger f4 for the same position. f2 with the sampling terms removed; the data are
    # pseudo-haploid, so the number of independent chromosomes is half the coded count.
    def f2(ka, kb):
        ca, ta = (A[ka] if ka in A else P[ka][:2]); cb, tb = (A[kb] if kb in A else P[kb][:2])
        na, nb = ta / 2.0, tb / 2.0
        ok = auto & (na >= 4) & (nb >= 4)
        with np.errstate(invalid="ignore", divide="ignore"):
            pa, pb = ca / ta, cb / tb
            v = (pa - pb) ** 2 - pa * (1 - pa) / (na - 1) - pb * (1 - pb) / (nb - 1)
        return float(np.nanmean(v[ok]))
    print("\n## Length of the cline (f2 ×1000) and the position f4 / f2 of the Piast core relative to Santok + Ląd\n")
    print("| Central Europe | Balts | f2 ×1000 | f4 ×1000 | position = f4 / f2 |")
    print("|---|---|---|---|---|")
    for ce in ("Germany 400-1100 CE (AADR)",):
        for balt in ("Balts LT only (1st-6th c.)", "Balts LV only (Kivutkalns, capture)", "Balts LT+LV 800 BCE-600 CE"):
            ok = auto & ~np.isnan(F[ce]) & ~np.isnan(F[balt]) & ~np.isnan(pS) & ~np.isnan(pX)
            e1, s1, _ = f4_jack(F[ce], F[balt], pS, pX, ok, block)
            l = f2(ce, balt)
            print(f"| {ce.split(' (')[0]} | {balt} | {1000 * l:.3f} | {1000 * e1:+.3f} ± {1000 * s1:.3f} | {e1 / l:+.3f} ± {s1 / l:.3f} |")
    # ---- f4 ratio with the Baltic group split in two: the own drift of the Baltic group cancels.
    #   position = f4(CE, Balts_A; ref, X) / f4(CE, Balts_A; ref, Balts_B)     0 = the reference, 1 = the Balts
    # (f4 / f2 above is too small, because f2 holds the Baltic group's own drift and the numerator does not.)
    def sub_freq(idx):
        sub = np.stack([per_ind[i] for i in idx], axis=1); called = sub != 3
        t = 2.0 * called.sum(1)
        return freq(np.where(called, sub, 0).sum(1).astype(float), t)

    def ratio_jack(p1, p2, pS, pX, pB, ok):
        num = ((p1 - p2) * (pS - pX))[ok]; den = ((p1 - p2) * (pS - pB))[ok]
        ub, inv = np.unique(block[ok], return_inverse=True)
        sn, sd = np.bincount(inv, weights=num), np.bincount(inv, weights=den)
        est = sn.sum() / sd.sum(); loo = (sn.sum() - sn) / (sd.sum() - sd); g = len(ub)
        return est, float(np.sqrt((g - 1) / g * ((loo - loo.mean()) ** 2).sum())), sd.sum() / ok.sum(), int(ok.sum())
    print("\n## f4 ratio: position between the reference (0) and the Balts (1), Baltic group split in two halves\n")
    print("| Balts | reference | X | n | position (mean of 20 random splits, both directions) | jackknife SE of one split (median) | spread over splits (SD) |")
    print("|---|---|---|---|---|---|---|")
    rng2 = np.random.default_rng(11)
    ce = F["Germany 400-1100 CE (AADR)"]
    for bl in ("Balts LT only (1st-6th c.)", "Balts LT+LV 800 BCE-600 CE", "Balts LV only (Kivutkalns, capture)"):
        ids_b = groups[bl]
        for refk, rl in ((("Santok", "Ląd"), "Santok + Ląd"), (rest, "Santok + Niemcza + Ląd")):
            pSr, _ = pooled2(refk)
            for x in ("Piast core, pooled", "Ostrów Lednicki", "Płońsk", "Dziekanowice", "Końskie", "Niemcza"):
                if x == "Niemcza" and "Niemcza" in refk:
                    continue
                ests, ses = [], []
                for _ in range(20):
                    perm = rng2.permutation(len(ids_b)); h = len(ids_b) // 2
                    for ia, ib in ((perm[:h], perm[h:]), (perm[h:], perm[:h])):
                        pA, pB = sub_freq([ids_b[i] for i in ia]), sub_freq([ids_b[i] for i in ib])
                        ok = auto & ~np.isnan(ce) & ~np.isnan(pA) & ~np.isnan(pB) & ~np.isnan(pSr) & ~np.isnan(F[x])
                        e, se, _d, _n = ratio_jack(ce, pA, pSr, F[x], pB, ok)
                        ests.append(e); ses.append(se)
                print(f"| {bl} | {rl} | {x} | {N[x]} | {np.mean(ests):+.3f} | {np.median(ses):.3f} | {np.std(ests):.3f} |")
                rows.append((f"f4 ratio / {bl}", f"{rl} -> {x}", N[x], np.mean(ests), np.median(ses), np.mean(ests) / np.median(ses), 0, 0, ""))
    with open(os.path.join(a.out, "f4_sites.tsv"), "w") as f:
        f.write("balts\tX\tn\tf4\tse\tz\tsnps\tz_transversions\tindividuals_positive\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
