#!/usr/bin/env python3
"""Add gnomAD exome allele frequency from the existing VEP 116 output.

Uses the same transcript as the MANE annotation (MANE Select, else canonical).
Does not drop any of the 380 rows.
"""
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib.machinery import SourceFileLoader

ann = SourceFileLoader("mane_ann", str(Path(__file__).resolve().parent / "37_mane_annotation.py")).load_module()
from common import PROJECT_DIR, write_provenance  # noqa: E402

src = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
v = pd.read_csv(src, sep="\t")
if len(v) != 380:
    raise SystemExit(f"[gnomad] expected 380 rows, found {len(v)}")
caches = {p: ann.load_vep(p) for p in v["patient"].unique()}
afs = []
unmatched = []
for r in v.itertuples(index=False):
    key = (str(r.chrom), int(r.pos), str(r.ref), str(r.alt))
    blocks = caches[r.patient].get(key)
    if blocks is None:
        blocks = caches[r.patient].get((key[0], key[1] - 1, key[2], key[3]))
    if not blocks:
        unmatched.append(r.var_id)
        afs.append("")
        continue
    chosen, _rule = ann.choose(blocks)
    afs.append(chosen.get("gnomADe_AF") or "")
if unmatched:
    raise SystemExit(f"[gnomad] unmatched {unmatched[:5]}")
v["gnomAD_AF"] = afs
v.to_csv(src, sep="\t", index=False)
write_provenance(src, [PROJECT_DIR / "06_variants" / "vep_hgvs" / "DH00018.pass.vep.vcf.gz"],
                 {"field": "gnomADe_AF", "transcript": "MANE Select else canonical"})
num = pd.to_numeric(v["gnomAD_AF"], errors="coerce")
hi = v.loc[num > 1e-3, ["var_id", "patient", "gene", "gnomAD_AF"]]
print(f"[gnomad] column added. above 0.001: {len(hi)} of 380")
if len(hi):
    print(hi.to_string(index=False))
