#!/usr/bin/env python3
"""Force-count alleles at every candidate somatic locus in ALL 15 BAMs (look-back).

Candidates = variants passing call filters in T and/or J (06_variants/somatic_filtered.tsv).
Counting all 15 samples gives (i) the owner patient's N/J/T evidence and (ii) a background
error panel from the 12 samples of the 4 unrelated patients at the same position.

Read filters: no duplicates/secondary/supplementary/QC-fail, MAPQ >= lookback_min_mapq,
base quality >= lookback_min_baseq (SNVs). Per-strand counts and alt-read position are
kept for artefact flags.

Output: 08_lookback_N/allele_counts.tsv
"""
import argparse
import statistics
import sys
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PATIENTS, PROJECT_DIR, load_env, load_filters, write_provenance  # noqa: E402

env = load_env()
F = load_filters()
ap = argparse.ArgumentParser()
ap.add_argument("--variants", default=str(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"))
ap.add_argument("--bam-dir", default=env["SAREK_OUTDIR"], help="searched recursively for <sample>.recal.bam/cram")
ap.add_argument("--reference", default=None, help="FASTA, required only for CRAM input")
ap.add_argument("--out", default=str(PROJECT_DIR / "08_lookback_N" / "allele_counts.tsv"))
args = ap.parse_args()

MINMQ, MINBQ = F["lookback_min_mapq"], F["lookback_min_baseq"]


def find_bam(sid: str) -> Path | None:
    for pat in (f"{sid}.recal.bam", f"{sid}.recal.cram", f"{sid}.bam"):
        hits = list(Path(args.bam_dir).rglob(pat))
        if hits:
            return hits[0]
    return None


def count(bam, chrom, pos, ref, alt):
    """Return dict of counts at a VCF-style variant (1-based pos)."""
    is_snv = len(ref) == 1 and len(alt) == 1
    indel_len = len(alt) - len(ref)
    c = dict(depth=0, alt=0, alt_fwd=0, alt_rev=0, ref_fwd=0, ref_rev=0)
    dists, bqs = [], []
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter",
                          min_base_quality=0, ignore_overlaps=True, max_depth=100000):
        for pr in col.pileups:
            aln = pr.alignment
            if (aln.is_duplicate or aln.is_secondary or aln.is_supplementary or aln.is_qcfail
                    or aln.mapping_quality < MINMQ or pr.is_refskip):
                continue
            if pr.is_del:   # read spans a deletion at this base: uninformative here
                continue
            qpos = pr.query_position
            bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 40
            if is_snv:
                if bq < MINBQ:
                    continue
                base = aln.query_sequence[qpos].upper()
                is_alt, is_ref = base == alt.upper(), base == ref.upper()
            else:
                is_alt = pr.indel == indel_len and indel_len != 0
                is_ref = pr.indel == 0
            if not (is_alt or is_ref):
                c["depth"] += 1        # other allele: counts towards depth only
                continue
            c["depth"] += 1
            strand = "rev" if aln.is_reverse else "fwd"
            if is_alt:
                c["alt"] += 1
                c[f"alt_{strand}"] += 1
                dists.append(min(qpos, aln.query_length - qpos - 1))
                bqs.append(bq)
            else:
                c[f"ref_{strand}"] += 1
    c["alt_median_dist_end"] = statistics.median(dists) if dists else None
    c["alt_mean_bq"] = round(sum(bqs) / len(bqs), 1) if bqs else None
    return c


v = pd.read_csv(args.variants, sep="\t")
loci = v[["var_id", "patient", "chrom", "pos", "ref", "alt"]].drop_duplicates("var_id")
rows, used = [], []
for p in PATIENTS:
    for s in ("N", "J", "T"):
        sid = f"{p}_{s}"
        path = find_bam(sid)
        if path is None:
            print(f"[WARN] BAM not found for {sid}")
            continue
        used.append(path)
        bam = pysam.AlignmentFile(str(path), reference_filename=args.reference)
        for r in loci.itertuples():
            c = count(bam, r.chrom, int(r.pos), r.ref, r.alt)
            rows.append({"var_id": r.var_id, "variant_patient": r.patient, "sample_id": sid,
                         "sample_patient": p, "specimen": s, **c})
        bam.close()
        print(f"[count] {sid}: {len(loci)} loci")

out = Path(args.out)
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
write_provenance(out, [Path(args.variants)], {"min_mapq": MINMQ, "min_baseq": MINBQ})
print(f"[count] wrote {out}")
