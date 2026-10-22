#!/usr/bin/env python3
"""30_build_pca.py -- transparent G25-style PCA from AADR + projection of targets.

Basis: present-day West-Eurasian Human Origins individuals of the AADR (".HO",
date 0, lat 10-75, lon -25-60, no outlier/QC labels), at most CAP per population.
PCA: SNPs standardised smartpca-style, x = (g - 2p) / sqrt(2p(1-p)), missing = 0;
eigen-decomposition of the n x n Gram matrix; 25 PCs.
Projection (smartpca lsqproject): every other individual -- ancient pseudo-
haploid, modern, or a personal genome -- gets s = argmin ||y_o - V_o s||^2 over
ITS OWN observed SNPs only, so low coverage widens the error but does not pull
the sample toward the origin (the shrinkage bug of mean-imputed projection).
Coordinates are reported like G25 "scaled" values: the score itself (variance =
eigenvalue) divided by sqrt(M) so magnitudes are comparable across runs.

Reads the SNP-major GENO directly (2 bits per genotype, high bits first,
0/1/2 = copies of the .snp column-5 allele, 3 = missing). Never point this at a
TGENO file -- it checks the header.

Usage:
  python3 -I 30_build_pca.py --prefix $AADR_GENO/v66.p1_HO.aadr.patch.PUB \
      --anno $AADR_GENO/v66.p1_HO.aadr.PUB.anno \
      --target PW=./PW_WGS_HO/pw:wgs \
      --target JW=./JW_MH/pw \
      --relabel relabel_regions.tsv --out ./PW_WGS_HO/pca
Targets are step-01/01b PED/MAP stems (<stem>.ped, <stem>.map). Append ":wgs" for a
forward-strand WGS call set (strand-ambiguous A/T, C/G SNPs kept); without it the
target is treated as a consumer array (ambiguous SNPs dropped, others strand-flipped
to the panel alleles when needed).
Outputs (in --out): individuals.csv, populations.csv (G25-like
"label,PC1..PC25"), targets.csv, eigen.tsv, basis.tsv, summary.txt.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

K = 25
CHUNK = 8192
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
BAD_LABEL = ("-o", "Ignore", "outlier", "Contam", "lowcov", "Dup", "Rel", "QCremove")


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def read_anno(path: str) -> dict:
    out = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 28:
                continue
            def num(x):
                try:
                    return float(x)
                except ValueError:
                    return math.nan
            out[c[0]] = dict(group=c[14], country=c[16], lat=num(c[17]), lon=num(c[18]),
                             bp=num(c[10]), ho_snps=num(c[27]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--anno", required=True)
    ap.add_argument("--target", action="append", default=[], help="LABEL=ped/map stem")
    ap.add_argument("--relabel", default="", help="TSV IID<TAB>label overriding AADR groups")
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--min-snps", type=int, default=15000, help="min HO SNPs for projected ancients")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    geno = a.prefix + ".geno"
    with open(geno, "rb") as f:
        head = f.read(48).split(b"\0")[0].split()
    if head[0] != b"GENO":
        die(f"{geno} header is {head[0]!r}, not GENO (convert TGENO with aadr/tgeno2geno.py)")
    nind, nsnp = int(head[1]), int(head[2])
    rlen = max(48, math.ceil(nind / 4))
    if os.path.getsize(geno) != rlen * (nsnp + 1):
        die("GENO size does not match its header")

    ind = [l.split() for l in open(a.prefix + ".ind")]
    snp = [l.split() for l in open(a.prefix + ".snp")]
    if len(ind) != nind or len(snp) != nsnp:
        die(".ind/.snp sizes do not match the GENO header")
    iid = [r[0] for r in ind]
    pop = [r[2] for r in ind]
    chrom = np.array([r[1] for r in snp])
    auto = np.isin(chrom, [str(i) for i in range(1, 23)])
    a1 = [r[4] for r in snp]
    a2 = [r[5] for r in snp]
    snp_index = {r[0]: i for i, r in enumerate(snp)}

    anno = read_anno(a.anno)
    relabel = {}
    if a.relabel:
        for line in open(a.relabel):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != "IID":
                relabel[p[0]] = p[1]

    # ---- basis selection --------------------------------------------------------
    by_pop = defaultdict(list)
    for i, x in enumerate(iid):
        m = anno.get(x)
        if not m or not x.endswith(".HO") or m["bp"] != 0:
            continue
        if not (10 <= m["lat"] <= 75 and -25 <= m["lon"] <= 60):
            continue
        if any(t in pop[i] for t in BAD_LABEL):
            continue
        by_pop[pop[i]].append(i)
    basis = sorted(i for p, v in by_pop.items() for i in v[: a.cap])
    # ---- projection selection: everything West-Eurasian worth projecting ------
    proj = []
    for i, x in enumerate(iid):
        m = anno.get(x)
        if not m or any(t in pop[i] for t in ("Ignore", "Contam")):
            continue
        if not (25 <= m["lat"] <= 75 and -25 <= m["lon"] <= 80):
            continue
        if m["bp"] > 0 and not (m["ho_snps"] >= a.min_snps):
            continue
        proj.append(i)
    proj = sorted(set(proj) | set(basis))
    sel = np.array(proj)
    pos_in_sel = {i: k for k, i in enumerate(proj)}
    bsel = np.array([pos_in_sel[i] for i in basis])
    print(f"basis: {len(basis)} individuals in {len(by_pop)} populations; projected: {len(proj)}")

    mm = np.memmap(geno, dtype=np.uint8, mode="r", offset=rlen, shape=(nsnp, rlen))
    byte = sel // 4
    shift = (6 - 2 * (sel % 4)).astype(np.uint8)

    def chunk_geno(s0: int, s1: int) -> np.ndarray:
        blk = np.asarray(mm[s0:s1])[:, byte]
        return (blk >> shift) & 3                           # snps x proj

    # ---- pass 1: allele frequencies + Gram matrix on the basis -------------------
    nb = len(basis)
    keep = np.zeros(nsnp, bool)
    p_all = np.zeros(nsnp)
    gram = np.zeros((nb, nb))
    for s0 in range(0, nsnp, CHUNK):
        s1 = min(s0 + CHUNK, nsnp)
        g = chunk_geno(s0, s1)[:, bsel].astype(np.float64)
        obs = g != 3
        n_obs = obs.sum(1)
        p = np.where(n_obs > 0, np.where(obs, g, 0).sum(1) / np.maximum(2 * n_obs, 1), 0)
        ok = auto[s0:s1] & (n_obs >= 0.95 * nb) & (p >= 0.01) & (p <= 0.99)
        keep[s0:s1] = ok
        p_all[s0:s1] = p
        if ok.any():
            sd = np.sqrt(2 * p[ok] * (1 - p[ok]))
            x = np.where(obs[ok], (g[ok] - 2 * p[ok][:, None]) / sd[:, None], 0.0)
            gram += x.T @ x
    M = int(keep.sum())
    print(f"SNPs kept for PCA: {M} of {nsnp}")
    w, U = np.linalg.eigh(gram)
    order = np.argsort(w)[::-1][:K]
    lam, U = w[order], U[:, order]

    # ---- pass 2: SNP loadings + projection of everyone ---------------------------
    nproj = len(proj)
    A = np.zeros((nproj, K * K))
    b = np.zeros((nproj, K))
    nobs = np.zeros(nproj)

    targets = []
    tvec = []
    for spec in a.target:
        lab, stem = spec.split("=", 1)
        wgs = stem.endswith(":wgs")
        stem = stem[:-4] if wgs else stem
        g = np.full(nsnp, 3, np.int8)
        ped = open(stem + ".ped").read().split()
        al = ped[6:]
        flips = amb = bad = 0
        for k, line in enumerate(open(stem + ".map")):
            sid = line.split()[1]
            j = snp_index.get(sid)
            if j is None:
                continue
            x, y = al[2 * k], al[2 * k + 1]
            if "0" in (x, y):
                continue
            pair = {a1[j], a2[j]}
            if pair in ({"A", "T"}, {"C", "G"}) and not wgs:
                amb += 1                       # strand-ambiguous on an array: skip
                continue
            if not {x, y} <= pair:
                x, y = COMP.get(x, "?"), COMP.get(y, "?")
                if not {x, y} <= pair:
                    bad += 1
                    continue
                flips += 1
            g[j] = (x == a1[j]) + (y == a1[j])
        targets.append(lab)
        tvec.append(g)
        print(f"target {lab}: {int((g != 3).sum())} panel SNPs (flipped {flips}, ambiguous skipped {amb}, unmatched {bad})")
    T = np.array(tvec) if tvec else np.zeros((0, nsnp), np.int8)
    At = np.zeros((len(targets), K * K))
    bt = np.zeros((len(targets), K))
    ntobs = np.zeros(len(targets))

    for s0 in range(0, nsnp, CHUNK):
        s1 = min(s0 + CHUNK, nsnp)
        ok = keep[s0:s1]
        if not ok.any():
            continue
        g = chunk_geno(s0, s1)[ok].astype(np.float32)
        p = p_all[s0:s1][ok].astype(np.float32)
        sd = np.sqrt(2 * p * (1 - p))
        obs = g != 3
        y = np.where(obs, (g - 2 * p[:, None]) / sd[:, None], 0.0).astype(np.float32)
        xb = y[:, bsel]
        V = (xb @ U.astype(np.float32)) / np.sqrt(lam).astype(np.float32)   # m x K, orthonormal-ish
        VV = (V[:, :, None] * V[:, None, :]).reshape(len(V), K * K)
        of = obs.astype(np.float32)
        A += of.T @ VV
        b += y.T @ V
        nobs += of.sum(0)
        if len(targets):
            tg = T[:, s0:s1][:, ok].astype(np.float32)
            tobs = tg != 3
            ty = np.where(tobs, (tg - 2 * p[None, :]) / sd[None, :], 0.0).astype(np.float32)
            tof = tobs.astype(np.float32)
            At += tof @ VV
            bt += ty @ V
            ntobs += tof.sum(1)

    def solve(Ai, bi):
        return np.linalg.solve(Ai.reshape(K, K), bi)

    scale = 1 / math.sqrt(M)
    S = np.array([solve(A[i], b[i]) for i in range(nproj)]) * scale
    St = np.array([solve(At[i], bt[i]) for i in range(len(targets))]) * scale if len(targets) else np.zeros((0, K))

    # ---- check: projected basis = PCA scores --------------------------------------
    direct = (U * np.sqrt(lam)) * scale
    dev = np.abs(np.abs(S[bsel]) - np.abs(direct)).max()

    def lab_of(i):
        return relabel.get(iid[i], pop[i])

    cols = ",".join(f"PC{k+1}" for k in range(K))
    with open(os.path.join(a.out, "individuals.csv"), "w") as f:
        f.write(f"label:iid,{cols},snps\n")
        for k, i in enumerate(proj):
            f.write(f"{lab_of(i)}:{iid[i]}," + ",".join(f"{v:.6f}" for v in S[k]) + f",{int(nobs[k])}\n")
    groups = defaultdict(list)
    for k, i in enumerate(proj):
        groups[lab_of(i)].append(k)
    with open(os.path.join(a.out, "populations.csv"), "w") as f:
        f.write(f"label,{cols},n\n")
        for g, ks in sorted(groups.items()):
            f.write(g + "," + ",".join(f"{v:.6f}" for v in S[ks].mean(0)) + f",{len(ks)}\n")
    with open(os.path.join(a.out, "targets.csv"), "w") as f:
        f.write(f"label,{cols},snps\n")
        for k, t in enumerate(targets):
            f.write(t + "," + ",".join(f"{v:.6f}" for v in St[k]) + f",{int(ntobs[k])}\n")
    with open(os.path.join(a.out, "eigen.tsv"), "w") as f:
        f.write("pc\teigenvalue\tvariance_share\n")
        for k in range(K):
            f.write(f"{k+1}\t{lam[k]:.4f}\t{lam[k]/np.trace(gram):.5f}\n")
    with open(os.path.join(a.out, "basis.tsv"), "w") as f:
        for pname, v in sorted(by_pop.items()):
            f.write(f"{pname}\t{min(len(v), a.cap)}\n")
    msg = (f"basis {nb} individuals / {len(by_pop)} populations; SNPs {M}; projected {nproj}; "
           f"targets {', '.join(targets)}; PC1-3 variance share "
           f"{', '.join(f'{lam[k]/np.trace(gram):.4f}' for k in range(3))}; "
           f"max |projected - direct| basis score deviation {dev:.2e}")
    open(os.path.join(a.out, "summary.txt"), "w").write(msg + "\n")
    print(msg)


if __name__ == "__main__":
    main()
