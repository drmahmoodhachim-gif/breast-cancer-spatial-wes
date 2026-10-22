#!/usr/bin/env python3
"""Add the autosomal fraction with major copy number (tcn - lcn) >= 2."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

AUTO = {str(i) for i in range(1, 23)}
pp_path = PROJECT_DIR / "07_cnv_purity" / "purity_ploidy.tsv"
pp = pd.read_csv(pp_path, sep="\t")
segs = pd.concat(
    pd.read_csv(f, sep="\t") for f in (PROJECT_DIR / "07_cnv_purity" / "facets").glob("*/*.segments.tsv")
)
segs["chrom"] = segs["chrom"].astype(str).str.replace("chr", "", regex=False)
segs = segs[segs["chrom"].isin(AUTO)].copy()
segs["length"] = (segs["end"] - segs["start"]).clip(lower=0)
segs["major"] = segs["tcn.em"] - segs["lcn.em"]


def frac(g):
    total = g["length"].sum()
    if total <= 0:
        return float("nan")
    hit = g.loc[g["major"] >= 2, "length"].sum()
    return float(hit / total)


values = {sid: frac(g) for sid, g in segs.groupby("sample_id")}
pp["frac_autosome_major_cn_ge2"] = pp["sample_id"].map(values)
pp.to_csv(pp_path, sep="\t", index=False)
write_provenance(pp_path, list((PROJECT_DIR / "07_cnv_purity" / "facets").glob("*/*.segments.tsv")),
                 {"column": "frac_autosome_major_cn_ge2", "rule": "autosomal bases with tcn-lcn >= 2"})
print(pp[["sample_id", "purity", "ploidy", "WGD_heuristic", "frac_autosome_major_cn_ge2"]].to_string(index=False))
