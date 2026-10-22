#!/usr/bin/env python3
"""Quick pre-flight check on the first N reads of every FASTQ (no full read-through).

Reports per file: header fields (instrument/run/flowcell/lane/index), read length,
mean base quality, %Q30, %N, and whether R1/R2 read names are in sync.
Output: 02_fastq_check/fastq_headcheck.tsv

NOTE: OneDrive "Files On-Demand" may download the whole file when it is first opened.
Mark the folder "Always keep on this device" before running the pipeline anyway.
"""
import argparse
import gzip
import sys
from itertools import islice
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--n-reads", type=int, default=100_000)
args = ap.parse_args()

raw = PROJECT_DIR / "01_raw_fastq"
out = PROJECT_DIR / "02_fastq_check" / "fastq_headcheck.tsv"


def stats(path: Path, n: int):
    names, lens, qsum, qn, q30, nn = [], [], 0, 0, 0, 0
    header = None
    with gzip.open(path, "rt") as fh:
        while len(names) < n:
            rec = list(islice(fh, 4))
            if len(rec) < 4:
                break
            h, seq, _, qual = (x.rstrip("\n") for x in rec)
            if header is None:
                header = h
            names.append(h.split()[0])
            lens.append(len(seq))
            nn += seq.upper().count("N")
            for c in qual:
                q = ord(c) - 33
                qsum += q
                q30 += q >= 30
            qn += len(qual)
    f = header[1:].split()[0].split(":") if header else []
    idx = header.split()[1].split(":")[-1] if header and len(header.split()) > 1 else ""
    return {
        "instrument": f[0] if len(f) > 0 else "", "run": f[1] if len(f) > 1 else "",
        "flowcell": f[2] if len(f) > 2 else "", "lane": f[3] if len(f) > 3 else "", "index": idx,
        "reads_checked": len(names), "read_len_min": min(lens) if lens else 0,
        "read_len_max": max(lens) if lens else 0,
        "mean_q": round(qsum / qn, 2) if qn else 0, "pct_q30": round(100 * q30 / qn, 2) if qn else 0,
        "pct_N": round(100 * nn / max(1, sum(lens)), 3),
    }, names


rows = []
for r1 in sorted(raw.glob("*_R1_001.fastq.gz")):
    r2 = Path(str(r1).replace("_R1_001", "_R2_001"))
    s1, n1 = stats(r1, args.n_reads)
    s2, n2 = stats(r2, args.n_reads)
    in_sync = all(a.rsplit("/", 1)[0] == b.rsplit("/", 1)[0] for a, b in zip(n1, n2))
    for name, s in ((r1.name, s1), (r2.name, s2)):
        rows.append({"file": name, **s, "pairs_in_sync": in_sync})
    print(f"[headcheck] {r1.name}: len {s1['read_len_max']}, Q30 {s1['pct_q30']}%/{s2['pct_q30']}%, sync={in_sync}")

if not rows:
    sys.exit("[headcheck] no R1/R2 pairs under 01_raw_fastq")
cols = list(rows[0].keys())
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
write_provenance(out, params={"n_reads": args.n_reads})
print(f"[headcheck] wrote {out}")
