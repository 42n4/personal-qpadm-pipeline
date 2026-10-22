#!/usr/bin/env python3
"""tgeno2geno.py -- convert an EIGENSTRAT packed TGENO file to packed GENO.

AADR v66.p1 ships its .geno files in the TRANSPOSED packed format (header
"TGENO": one record per INDIVIDUAL). admixtools 2.0.10 read_packedancestrymap()
only understands SNP-major "GENO" and silently decodes a TGENO file as if it
were GENO -- every genotype it returns is garbage (2026-10-21: Mbuti and Han
came out with identical allele frequencies, correlation 0.03, 38 % missing).

Format (both variants, EIGENSOFT packedancestrymap):
  header record: ASCII "<magic> nind nsnp ihash shash", zero-padded to rlen
  data: 2 bits per genotype, 4 per byte, FIRST item in the HIGH bits;
        0/1/2 = copies of the first allele, 3 = missing
  TGENO: a FIXED 48-byte header, then nind records of ceil(nsnp/4) bytes
         (checked on v66.p1 HO and 1240K: size = 48 + nind*ceil(nsnp/4))
  GENO : nsnp records of rlen_g = max(48, ceil(nind/4)) bytes
The ind/snp hashes are copied unchanged (same .ind/.snp files).

Usage:
  python3 -I tgeno2geno.py <in.geno (TGENO)> <out.geno>
Then point a prefix at it, e.g. a directory with out.geno plus symlinks to the
original .snp/.ind/.anno under the same stem.
"""
import math
import os
import sys

import numpy as np


def main() -> None:
    if len(sys.argv) != 3 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    if os.path.exists(dst):
        sys.exit(f"refusing to overwrite {dst}")

    with open(src, "rb") as f:
        head = f.read(48).split(b"\0")[0].split()
    if not head or head[0] != b"TGENO":
        sys.exit(f"{src}: header is {head[:1]}, not TGENO -- nothing to convert")
    nind, nsnp = int(head[1]), int(head[2])
    ihash, shash = head[3].decode(), head[4].decode()

    rlen_t = math.ceil(nsnp / 4)
    rlen_g = max(48, math.ceil(nind / 4))
    expect = 48 + rlen_t * nind
    if os.path.getsize(src) != expect:
        sys.exit(f"{src}: size {os.path.getsize(src)} != expected {expect} for TGENO {nind}x{nsnp}")

    mm = np.memmap(src, dtype=np.uint8, mode="r", offset=48, shape=(nind, rlen_t))
    tmp = dst + ".part"
    chunk_bytes = 2048                      # 8192 SNPs per pass
    pad_ind = rlen_g * 4 - nind
    with open(tmp, "wb") as out:
        hdr = f"GENO {nind:7d} {nsnp:7d} {ihash} {shash}".encode()
        out.write(hdr + b"\0" * (rlen_g - len(hdr)))
        done = 0
        for b0 in range(0, math.ceil(nsnp / 4), chunk_bytes):
            b1 = min(b0 + chunk_bytes, math.ceil(nsnp / 4))
            blk = np.asarray(mm[:, b0:b1])                      # nind x nb
            g = np.empty((nind, (b1 - b0) * 4), dtype=np.uint8)
            g[:, 0::4] = (blk >> 6) & 3
            g[:, 1::4] = (blk >> 4) & 3
            g[:, 2::4] = (blk >> 2) & 3
            g[:, 3::4] = blk & 3
            n_here = min(g.shape[1], nsnp - b0 * 4)
            t = np.ascontiguousarray(g[:, :n_here].T)           # snp x nind
            if pad_ind:
                t = np.hstack([t, np.full((n_here, pad_ind), 3, dtype=np.uint8)])
            packed = (t[:, 0::4] << 6) | (t[:, 1::4] << 4) | (t[:, 2::4] << 2) | t[:, 3::4]
            out.write(packed.astype(np.uint8).tobytes())
            done += n_here
        if done != nsnp:
            sys.exit(f"converted {done} SNPs, expected {nsnp}")
    want = rlen_g * (nsnp + 1)
    if os.path.getsize(tmp) != want:
        sys.exit(f"{tmp}: size {os.path.getsize(tmp)} != expected {want}")
    os.replace(tmp, dst)
    print(f"{dst}: GENO {nind} individuals x {nsnp} SNPs, {want} bytes")


if __name__ == "__main__":
    main()
