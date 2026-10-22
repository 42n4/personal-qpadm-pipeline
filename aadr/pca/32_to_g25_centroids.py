#!/usr/bin/env python3
"""32_to_g25_centroids.py -- feed the 30_build_pca.py coordinates to g25_nmonte_v10.py.

The nMonte runner has a fixed source catalogue (GroupDef entries) whose ids are
Davidski/EYD cluster names (Country_Period_Site_..._c1). This script rebuilds
each id it can from AADR individuals projected by 30_build_pca.py:
  country (first token) + period token (date window) + site token (matched in the
  AADR locality or group label, accent-insensitive). Sub-cluster ids of one site
  (_c1, _c2, ...) cannot be reproduced and get the whole-site mean -- the runner's
  alias dedup then treats them as one source. Ids with no AADR match are left out;
  the runner skips missing ids.

Usage:
  python3 -I 32_to_g25_centroids.py --pca DIR --anno ANNO --runner g25_nmonte_v10.py --out DIR
Outputs: centroids_ownpca.txt (runner format: id,PC1..PC25), mapping.tsv (audit:
id, n, sites used, rule), <TARGET>_ownpca_scaled.txt per target.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import unicodedata
from collections import defaultdict

import numpy as np

PERIODS = {   # token -> (min BP, max BP)
    "Neo": (5000, 9000), "CW": (4300, 4900), "BB": (4000, 4600), "BA": (2800, 4500),
    "IA": (1500, 2800), "Roman": (1500, 2100), "Mig": (1300, 1700), "Viking": (900, 1250),
    "Medieval": (450, 1450), "PostMed": (50, 500),
}
COUNTRY = {"Czech": "Czechia", "Czechia": "Czechia", "UK": "United Kingdom",
           "England": "United Kingdom", "Scotland": "United Kingdom"}
SKIP = {"England", "Scotland", "Flanders", "Czech", "Republic", "Bohemia", "Hun", "Medieval",
        "Sint", "Saint", "St", "Greater", "Lesser", "Lower", "Upper", "Saxony", "Zealand", "Jutland", "Viking"}
LABEL_PERIOD = {"Neo": "N", "Mig": "Migration"}
PERIOD_WORD = {"Neo": " n", "CW": "cordedware", "BB": "bellbeaker", "BA": " ba", "IA": " ia",
               "Roman": "roman", "Mig": "migration", "Viking": "viking", "Medieval": "medieval",
               "PostMed": "postmedieval"}
GENERIC = re.compile(r"^(c\d+)+$|^broad$|^b$|^\d+$")


def fold(s: str) -> str:
    s = s.replace("ł", "l").replace("Ł", "L")   # NFKD does not decompose ł
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()
    return s.lower().replace("-", " ").replace("_", " ")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pca", required=True)
    ap.add_argument("--anno", required=True)
    ap.add_argument("--runner", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    src = open(a.runner, encoding="utf-8").read()
    ids = sorted(set(i for grp in re.findall(r'GroupDef\(ids=\(([^)]*)\)', src)
                     for i in re.findall(r'"([^"]+)"', grp)))

    meta = {}
    with open(a.anno, encoding="utf-8", errors="replace") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            try:
                bp = float(c[10])
            except ValueError:
                continue
            meta[c[0]] = (c[16], bp, fold(c[15]), c[14])

    coords = {}
    with open(os.path.join(a.pca, "individuals.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            iid = row[0].split(":", 1)[1]
            coords[iid] = np.array([float(x) for x in row[1:26]])

    rows, audit = [], []
    for cid in ids:
        tok = cid.split("_")
        orig0 = tok[0]
        if tok[0] == "United" and len(tok) > 2:
            tok = ["UK"] + tok[2:]
        country = COUNTRY.get(tok[0], tok[0])
        if len(tok) < 2 or tok[1] not in PERIODS:
            audit.append((cid, 0, "", "period token not recognised"))
            continue
        lo, hi = PERIODS[tok[1]]
        cand = [t for t in tok[2:] if not GENERIC.match(t) and t not in SKIP]
        pool = [(iid, loc, fold(grp)) for iid, (ctry, bp, loc, grp) in meta.items()
                if iid in coords and ctry == country and lo <= bp <= hi and "-o" not in grp]
        site, members, rule = "", [], ""
        if cand:
            # Only the FIRST site token names the site; later tokens are county /
            # region words that match other sites (Gyor, Ternopil, Hrubieszow ...).
            t = cand[0]
            members = [iid for iid, loc, g in pool if fold(t) in loc or fold(t) in g]
            site, rule = t, f"site '{t}'"
        else:
            # Pure cluster/aggregate id (Poland_Viking_c1, Poland_BA_broad): AADR
            # labels starting with the same Country_Period only.
            pre = fold(f"{orig0}_{LABEL_PERIOD.get(tok[1], tok[1])}")
            members = [iid for iid, loc, g in pool if g == pre or g.startswith(pre + " ")]
            rule = f"label prefix '{orig0}_{LABEL_PERIOD.get(tok[1], tok[1])}'"
        sites = {meta[i][2].split(" (")[0][:30] for i in members}
        if not members:
            audit.append((cid, 0, "", f"no AADR match ({country}, {tok[1]}, {rule})"))
            continue
        m = np.mean([coords[i] for i in members], axis=0)
        rows.append((cid, m))
        audit.append((cid, len(members), "; ".join(sorted(sites))[:120], rule))

    with open(os.path.join(a.out, "centroids_ownpca.txt"), "w") as f:
        for cid, m in rows:
            f.write(cid + "," + ",".join(f"{v:.6f}" for v in m) + "\n")
    with open(os.path.join(a.out, "mapping.tsv"), "w") as f:
        f.write("id\tn\tsites\trule\n")
        for t in audit:
            f.write("\t".join(map(str, t)) + "\n")
    with open(os.path.join(a.pca, "targets.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            with open(os.path.join(a.out, f"{row[0]}_ownpca_scaled.txt"), "w") as g:
                g.write(f"{row[0]}_ownpca_scaled," + ",".join(row[1:26]) + "\n")
    print(f"catalogue ids: {len(ids)}; rebuilt from AADR: {len(rows)}; not available: {len(ids) - len(rows)}")


if __name__ == "__main__":
    main()
