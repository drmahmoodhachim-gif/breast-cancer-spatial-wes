#!/usr/bin/env python3
"""Flag segmental duplications and low-MAPQ reads. Draw KMT2C pileup snapshots.

Mappability tracks are not in the offline resources, so mappability stays NA.
Flags are added to the variant table. Nothing is removed.
"""
import sys
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_env, write_provenance  # noqa: E402

env = load_env()
RES = Path(env["SAREK_OUTDIR"])
SEG = PROJECT_DIR / "03_resources" / "genomicSuperDups.txt.gz"
URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"
VAR = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"


def load_segdups():
    if not SEG.exists() or SEG.stat().st_size < 1000:
        print(f"[artefact] downloading {URL}")
        urllib.request.urlretrieve(URL, SEG)
    raw = pd.read_csv(SEG, sep="\t", header=None, usecols=[1, 2, 3])
    raw.columns = ["chrom", "start", "end"]
    raw["chrom"] = raw["chrom"].astype(str)
    return {c: g.sort_values("start") for c, g in raw.groupby("chrom")}


def overlaps(chrom, pos, trees):
    g = trees.get(str(chrom))
    if g is None:
        g = trees.get("chr" + str(chrom).replace("chr", ""))
    if g is None:
        return False
    # UCSC start is 0-based; VCF pos is 1-based.
    hit = g[(g["start"] < pos) & (pos <= g["end"])]
    return not hit.empty


def find_bam(sid):
    path = RES / "preprocessing" / "recalibrated" / sid / f"{sid}.recal.bam"
    return path if path.exists() else None


def mapq_frac(bam, chrom, pos):
    n = low = 0
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter",
                          min_base_quality=0, ignore_overlaps=True, max_depth=100000):
        for pr in col.pileups:
            aln = pr.alignment
            if aln.is_duplicate or aln.is_secondary or aln.is_supplementary or aln.is_qcfail or pr.is_refskip or pr.is_del:
                continue
            n += 1
            if aln.mapping_quality < 50:
                low += 1
    return n, (low / n if n else float("nan"))


def snapshot(bams, chrom, pos, title, dest):
    fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharey=True)
    for ax, (sid, bam) in zip(axes, bams):
        counts = {b: 0 for b in "ACGT"}
        n = low = 0
        for col in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter",
                              min_base_quality=0, ignore_overlaps=True, max_depth=100000):
            for pr in col.pileups:
                aln = pr.alignment
                if aln.is_duplicate or aln.is_secondary or aln.is_supplementary or aln.is_qcfail or pr.is_del or pr.is_refskip:
                    continue
                qpos = pr.query_position
                if qpos is None:
                    continue
                base = aln.query_sequence[qpos].upper()
                if base in counts:
                    counts[base] += 1
                n += 1
                if aln.mapping_quality < 50:
                    low += 1
        ax.bar(list(counts), list(counts.values()), color=["#2a78d6", "#1baf7a", "#eb6834", "#52514e"])
        ax.set_title(f"{sid}\n{n} reads, {low} MAPQ<50", fontsize=8)
    fig.suptitle(title, fontsize=10)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140, bbox_inches="tight")
    plt.close(fig)


v = pd.read_csv(VAR, sep="\t")
trees = load_segdups()
v["segdup"] = [overlaps(c, int(p), trees) for c, p in zip(v["chrom"], v["pos"])]
v["mappability"] = "NA"
v["low_mappability"] = "NA"

fracs = {("T", i): float("nan") for i in v.index}
# filled below per row
t_frac, j_frac, n_frac = [], [], []
open_bams = {}
for patient in sorted(v["patient"].unique()):
    handles = {}
    for spec in ("N", "J", "T"):
        path = find_bam(f"{patient}_{spec}")
        if path is None:
            print(f"[artefact] missing BAM {patient}_{spec}")
            continue
        handles[spec] = pysam.AlignmentFile(str(path))
        open_bams[(patient, spec)] = handles[spec]
    sub = v[v["patient"] == patient]
    for idx, r in sub.iterrows():
        chrom = str(r["chrom"])
        if chrom not in (handles.get("T").references if "T" in handles else []):
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
        got = {}
        for spec, bam in handles.items():
            n, frac = mapq_frac(bam, chrom, int(r["pos"]))
            got[spec] = frac
        t_frac.append((idx, got.get("T", float("nan"))))
        j_frac.append((idx, got.get("J", float("nan"))))
        n_frac.append((idx, got.get("N", float("nan"))))
    print(f"[artefact] MAPQ {patient}: {len(sub)} sites")

v["frac_mapq_lt50_T"] = pd.Series({i: x for i, x in t_frac})
v["frac_mapq_lt50_J"] = pd.Series({i: x for i, x in j_frac})
v["frac_mapq_lt50_N"] = pd.Series({i: x for i, x in n_frac})
present_frac = v["frac_mapq_lt50_T"].where(v.get("T_present", False), v["frac_mapq_lt50_J"])
v["mapq_lt50_flag"] = (v[["frac_mapq_lt50_T", "frac_mapq_lt50_J"]].max(axis=1) > 0.10)
v["artefact_flag"] = v["segdup"] | v["mapq_lt50_flag"].fillna(False)

kmt = v[(v["patient"] == "DH00023") & (v["gene"] == "KMT2C")].copy()
snap_dir = PROJECT_DIR / "06_variants" / "kmt2c_snapshots"
fasta = env.get("REFERENCE_FASTA", "")
bam_paths = [str(find_bam(f"DH00023_{s}")) for s in ("N", "J", "T")]
batch = ["new", f"genome {fasta}", f"snapshotDirectory {snap_dir}"]
batch += [f"load {p}" for p in bam_paths if p != "None"]
for _, r in kmt.iterrows():
    chrom = str(r["chrom"])
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    label = f"{r['gene']}_{r.get('hgvsp_short', r['pos'])}".replace("*", "Ter").replace("/", "_")
    bams = [(f"DH00023_{s}", open_bams[("DH00023", s)]) for s in ("N", "J", "T") if ("DH00023", s) in open_bams]
    snapshot(bams, chrom, int(r["pos"]), f"DH00023 {label}", snap_dir / f"{label}.png")
    batch.append(f"goto {chrom}:{int(r['pos']) - 40}-{int(r['pos']) + 40}")
    batch.append(f"snapshot {label}.png")
(snap_dir / "igv_batch.txt").write_text("\n".join(batch) + "\n")
for bam in open_bams.values():
    bam.close()

v.to_csv(VAR, sep="\t", index=False)
write_provenance(VAR, [SEG], {"mappability": "NA, no offline track", "mapq_flag": ">10% reads MAPQ<50"})
print("[artefact] KMT2C")
cols = [c for c in ["gene", "hgvsp_short", "chrom", "pos", "segdup", "mappability", "frac_mapq_lt50_N",
                    "frac_mapq_lt50_J", "frac_mapq_lt50_T", "mapq_lt50_flag", "artefact_flag"] if c in kmt.columns or c in v.columns]
print(v.loc[kmt.index, cols].to_string(index=False))
print(f"[artefact] segdup {int(v['segdup'].sum())} mapq_flag {int(v['mapq_lt50_flag'].fillna(False).sum())}")
