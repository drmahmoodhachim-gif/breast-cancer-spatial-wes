#!/usr/bin/env python3
"""Pad (+/- BED_PADDING bp) and merge the Illumina Exome 2.0 Plus hg38 target BED.

Also reports raw and padded target size (Mb) -> used as TMB denominator.
Output: 03_resources/exome2plus_targets_hg38.pad100.bed and target_size.tsv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_env, write_provenance  # noqa: E402

env = load_env()
raw = Path(env["TARGET_BED_RAW"])
out = Path(env["TARGET_BED_PADDED"])
pad = int(env.get("BED_PADDING", 100))

if not raw.exists():
    sys.exit(f"[ERROR] Target BED not found: {raw}\n"
             "Download 'Twist for Illumina Exome 2.0 Plus Panel hg38 BED' from Illumina support "
             "and save it at that path.")

intervals = []
for line in raw.read_text().splitlines():
    if not line.strip() or line.startswith(("track", "browser", "#")):
        continue
    f = line.split("\t")
    chrom, start, end = f[0], int(f[1]), int(f[2])
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom  # GATK.GRCh38 uses chr-prefixed contigs
    intervals.append((chrom, start, end))

raw_bp = sum(e - s for _, s, e in intervals)


def chrom_key(c):
    c = c.replace("chr", "")
    return (0, int(c)) if c.isdigit() else (1, {"X": 23, "Y": 24, "M": 25}.get(c, 99), c)


intervals.sort(key=lambda x: (chrom_key(x[0]), x[1]))
merged = []
for chrom, s, e in intervals:
    s, e = max(0, s - pad), e + pad
    if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
        merged[-1][2] = max(merged[-1][2], e)
    else:
        merged.append([chrom, s, e])

with open(out, "w") as fh:
    for chrom, s, e in merged:
        fh.write(f"{chrom}\t{s}\t{e}\n")
pad_bp = sum(e - s for _, s, e in merged)

size_tsv = PROJECT_DIR / "03_resources" / "target_size.tsv"
size_tsv.write_text(
    "metric\tvalue\n"
    f"raw_intervals\t{len(intervals)}\nraw_target_mb\t{raw_bp/1e6:.3f}\n"
    f"padded_merged_intervals\t{len(merged)}\npadded_target_mb\t{pad_bp/1e6:.3f}\n"
    f"padding_bp\t{pad}\n")
write_provenance(out, [raw], {"padding": pad})
print(f"[bed] {len(intervals)} targets ({raw_bp/1e6:.2f} Mb) -> {len(merged)} padded ({pad_bp/1e6:.2f} Mb)")
