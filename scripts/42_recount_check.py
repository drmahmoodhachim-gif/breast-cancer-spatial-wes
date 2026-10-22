#!/usr/bin/env python3
"""Recount the one failing SNV without quality filters, and the indel misses with pysam.

Does not edit somatic_filtered.tsv. Prints whether each recount matches the table.
"""
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_env  # noqa: E402

ENV = load_env()
bam_root = Path(ENV.get("SAREK_OUTDIR", ""))
V = __import__("pandas").read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")


def bam_for(patient, spec):
    hits = list(bam_root.rglob(f"{patient}_{spec}.recal.bam"))
    return hits[0] if hits else None


def mpileup_snv(bam, chrom, pos, alt):
    import shutil
    import subprocess
    sam = shutil.which("samtools") or "samtools"
    cmd = [sam, "mpileup", "-r", f"{chrom}:{pos}-{pos}", str(bam)]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0 or not out.stdout.strip():
        return None, None
    parts = out.stdout.split("\t")
    depth = int(parts[3])
    bases = parts[4].upper()
    return bases.count(str(alt).upper()), depth


def indel_change(ref, alt):
    ref, alt = str(ref).upper(), str(alt).upper()
    if len(ref) > len(alt) and ref.startswith(alt):
        return -(len(ref) - len(alt))
    if len(alt) > len(ref) and alt.startswith(ref):
        return len(alt) - len(ref)
    return None


def pysam_indel(bam, chrom, pos, delta, min_mq, min_bq):
    pos0 = int(pos) - 1
    alt_n = depth = 0
    with pysam.AlignmentFile(str(bam), "rb") as af:
        for col in af.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="samtools",
                             min_mapping_quality=min_mq, min_base_quality=min_bq,
                             ignore_overlaps=False, ignore_orphans=True):
            if col.reference_pos != pos0:
                continue
            for pu in col.pileups:
                if pu.is_refskip or pu.is_del:
                    continue
                depth += 1
                if pu.indel == delta:
                    alt_n += 1
    return alt_n, depth


snv = V[V["var_id"] == "DH00018:chr11:85664329:C>A"].iloc[0]
jb = bam_for("DH00018", "J")
alt_n, depth = mpileup_snv(jb, snv.chrom, int(snv.pos), snv.alt)
tv = int(snv.J_alt) / int(snv.J_depth)
pv = alt_n / depth if depth else 0
print(f"SNV unfiltered J table {int(snv.J_alt)}/{int(snv.J_depth)} ({tv:.3f}) pileup {alt_n}/{depth} ({pv:.3f}) gap {abs(pv-tv):.3f}")
print("SNV_MATCH" if abs(pv - tv) <= 0.02 else "SNV_NO_MATCH")

indels = V[V["ref"].str.len() != V["alt"].str.len()]
# The audit flags a specimen when the quality-filtered mpileup alt is 0 and the table alt is > 0.
# Recount every indel specimen in the 30-site draw is unnecessary: recount every primary indel
# specimen whose table alt is > 0, and report those the previous audit listed.
wanted = {
    "DH00018:chr14:49655049:AAATTACCACACAG>A",
    "DH00023:chr16:23988552:TTC>T",
    "DH00018:chr1:153799594:GGCCGCTTGTCTCCTGGA>G",
    "DH00018:chr9:133568451:TG>T",
}
rows = indels[indels["var_id"].isin(wanted)]
ok = True
for r in rows.itertuples(index=False):
    delta = indel_change(r.ref, r.alt)
    for spec in ("N", "J", "T"):
        t_alt = int(getattr(r, f"{spec}_alt"))
        t_dep = int(getattr(r, f"{spec}_depth"))
        if t_alt <= 0:
            continue
        path = bam_for(r.patient, spec)
        a20, d20 = pysam_indel(path, r.chrom, r.pos, delta, 20, 20)
        a0, d0 = pysam_indel(path, r.chrom, r.pos, delta, 0, 0)
        match = a20 == t_alt or a0 == t_alt
        if not match:
            ok = False
        print(f"{r.var_id} {spec} table {t_alt}/{t_dep} pysamQ20 {a20}/{d20} pysam0 {a0}/{d0} {'MATCH' if match else 'NO'}")
print("INDEL_MATCH" if ok else "INDEL_NO_MATCH")
