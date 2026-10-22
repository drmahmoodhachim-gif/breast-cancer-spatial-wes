#!/usr/bin/env python3
"""Collect per-sample sequencing QC from sarek outputs into one table.

Inputs (found by glob under SAREK_OUTDIR, robust to minor version path changes):
  fastp JSON, MarkDuplicates metrics, mosdepth summary + region distribution.
Output: 04_qc/qc_summary.tsv  (one row per sample; feeds Fig 1 and Table S1)
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_env, parse_sample_id, write_provenance  # noqa: E402

env = load_env()
res = Path(env["SAREK_OUTDIR"])
out = PROJECT_DIR / "04_qc" / "qc_summary.tsv"
rows: dict[str, dict] = {}


def row(path: Path) -> dict | None:
    p, s = parse_sample_id(str(path))
    if p is None:
        return None
    return rows.setdefault(f"{p}_{s}", {"sample_id": f"{p}_{s}", "patient": p, "specimen": s})


used = []
for f in res.rglob("*fastp*.json"):
    r = row(f)
    if r is None:
        continue
    j = json.loads(f.read_text())["summary"]
    r["read_pairs_raw"] = j["before_filtering"]["total_reads"] // 2
    r["read_pairs_trimmed"] = j["after_filtering"]["total_reads"] // 2
    r["pct_q30_trimmed"] = round(100 * j["after_filtering"]["q30_rate"], 2)
    used.append(f)

for f in res.rglob("*.metrics"):
    if "markduplicates" not in str(f).lower() and ".md." not in f.name:
        continue
    r = row(f)
    if r is None:
        continue
    lines = f.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY") and i + 1 < len(lines):
            hdr, val = line.split("\t"), lines[i + 1].split("\t")
            d = dict(zip(hdr, val))
            r["pct_duplication"] = round(100 * float(d.get("PERCENT_DUPLICATION", "nan") or "nan"), 2)
            break
    used.append(f)

# Prefer the final (recalibrated) BAM's mosdepth if several exist
for f in sorted(res.rglob("*.mosdepth.summary.txt"), key=lambda x: "recal" in x.name):
    r = row(f)
    if r is None:
        continue
    df = pd.read_csv(f, sep="\t")
    reg = df[df["chrom"] == "total_region"]
    r["mean_target_depth"] = round(float((reg if len(reg) else df[df["chrom"] == "total"])["mean"].iloc[0]), 1)
    used.append(f)

for f in sorted(res.rglob("*.mosdepth.region.dist.txt"), key=lambda x: "recal" in x.name):
    r = row(f)
    if r is None:
        continue
    df = pd.read_csv(f, sep="\t", header=None, names=["chrom", "depth", "frac"])
    tot = df[df["chrom"] == "total"].set_index("depth")["frac"]
    for d in (20, 50, 100):
        r[f"pct_target_ge{d}x"] = round(100 * float(tot.get(d, float("nan"))), 1)

if not rows:
    sys.exit(f"[ERROR] No QC files found under {res}. Has sarek finished?")

qc = pd.DataFrame(rows.values()).sort_values(["patient", "specimen"])
qc.to_csv(out, sep="\t", index=False)
write_provenance(out, used[:50])
print(qc.to_string(index=False))

# Copy NGSCheckMate (sample identity) output alongside
for f in res.rglob("*ngscheckmate*/*"):
    if f.is_file() and f.suffix in {".txt", ".pdf", ".png"}:
        (PROJECT_DIR / "04_qc" / f.name).write_bytes(f.read_bytes())
