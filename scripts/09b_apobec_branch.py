#!/usr/bin/env python3
"""APOBEC-context SNVs by spatial branch.

An SNV is APOBEC-context when it is C>T or C>G in a TCW trinucleotide (W = A or T),
strand-collapsed, so G>A or G>C in WGA counts as well. Context bases come from the
reference FASTA. Each private branch is compared with the shared branch by Fisher's exact test.
Output: 10_tmb_signatures/apobec_by_branch.tsv
"""
import os
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

BRANCHES = ("shared", "T_private", "J_private")


def load_fasta():
    import pysam
    fasta = os.environ.get("REFERENCE_FASTA", "")
    if not fasta:
        env = PROJECT_DIR / "config" / "project.env"
        for line in env.read_text().splitlines():
            if line.startswith("REFERENCE_FASTA="):
                fasta = line.split("=", 1)[1].strip().strip('"')
    if not fasta or not Path(fasta).exists():
        raise SystemExit(f"REFERENCE_FASTA not found: {fasta}")
    return pysam.FastaFile(fasta)


def is_apobec(ref, alt, tri):
    ref, alt, tri = ref.upper(), alt.upper(), tri.upper()
    if len(tri) != 3 or tri[1] != ref or ref not in "ACGT" or alt not in "ACGT":
        return False
    if ref == "C" and alt in ("T", "G"):
        return tri[0] == "T" and tri[2] in ("A", "T")
    if ref == "G" and alt in ("A", "C"):
        return tri[0] in ("A", "T") and tri[2] == "A"
    return False


def context(fa, chrom, pos):
    name = str(chrom)
    if name not in fa.references and not name.startswith("chr") and f"chr{name}" in fa.references:
        name = f"chr{name}"
    if name not in fa.references and name.startswith("chr") and name[3:] in fa.references:
        name = name[3:]
    if name not in fa.references or pos < 2:
        return ""
    try:
        tri = fa.fetch(name, pos - 2, pos + 1)
    except (ValueError, KeyError):
        return ""
    return tri.upper()


V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
snv = V[(V["filter"].astype(str) == "PASS") & (V["variant_type"] == "SNV") & V["spatial_class"].isin(BRANCHES)].copy()
fa = load_fasta()
flags = []
for r in snv.itertuples(index=False):
    flags.append(is_apobec(r.ref, r.alt, context(fa, r.chrom, int(r.pos))))
snv["apobec"] = flags

rows = []
for patient, g in snv.groupby("patient"):
    counts = {}
    for branch in BRANCHES:
        sub = g[g["spatial_class"] == branch]
        n = int(len(sub))
        k = int(sub["apobec"].sum()) if n else 0
        counts[branch] = (k, n)
        rows.append({"patient": patient, "branch": branch, "n_snv": n, "n_apobec": k,
                     "frac": (k / n) if n else float("nan"), "fisher_p": float("nan")})
    shared_k, shared_n = counts["shared"]
    for row in rows:
        if row["patient"] != patient or row["branch"] == "shared":
            continue
        k, n = counts[row["branch"]]
        if n == 0 or shared_n == 0:
            continue
        table = [[k, n - k], [shared_k, shared_n - shared_k]]
        _, pval = fisher_exact(table, alternative="two-sided")
        row["fisher_p"] = pval

out = PROJECT_DIR / "10_tmb_signatures" / "apobec_by_branch.tsv"
tab = pd.DataFrame(rows)
tab.to_csv(out, sep="\t", index=False)
write_provenance(out, [PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"],
                 {"motif": "TCW C>T/C>G, strand-collapsed", "test": "Fisher private vs shared"})
print(tab.to_string(index=False))
print(f"[apobec] wrote {out}")
