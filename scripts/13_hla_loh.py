#!/usr/bin/env python3
"""HLA class I loss of heterozygosity from FACETS allele-specific copy number (exome resolution).

For each T and J sample, the FACETS segment covering HLA-A, -B and -C (GRCh38) is read:
  LOH         minor copy number = 0 with total copy number > 0  -> one HLA haplotype lost
  imbalance   minor CN > 0 but major CN > minor CN (allelic imbalance, partial loss / gain)
If the patient is HLA-homozygous at a gene (from 14_immune/hla/hla_types.tsv), LOH is not
informative for that gene and is reported as 'homozygous'.
Allele-specific HLA LOH (LOHHLA/DASH) is more precise and is listed as an optional follow-up.
Output: 14_immune/hla/hla_loh.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

HLA = {"HLA-A": 29_942_000, "HLA-C": 31_270_000, "HLA-B": 31_355_000}   # GRCh38, chr6
D = PROJECT_DIR / "07_cnv_purity" / "facets"
types_f = PROJECT_DIR / "14_immune" / "hla" / "hla_types.tsv"
types = pd.read_csv(types_f, sep="\t") if types_f.exists() else pd.DataFrame()

rows = []
for f in sorted(D.glob("*/*.segments.tsv")):
    s = pd.read_csv(f, sep="\t")
    sid = s["sample_id"].iloc[0]
    s["chrom"] = s["chrom"].astype(str).str.replace("chr", "", regex=False)
    for gene, pos in HLA.items():
        hit = s[(s["chrom"] == "6") & (s["start"] <= pos) & (s["end"] >= pos)]
        tcn = float(hit["tcn.em"].iloc[0]) if len(hit) else np.nan
        lcn = float(hit["lcn.em"].iloc[0]) if len(hit) and pd.notna(hit["lcn.em"].iloc[0]) else np.nan
        homo = False
        if len(types):
            n = types[types["sample_id"] == sid[:7] + "_N"]
            if len(n):
                k = gene[-1]
                homo = n.iloc[0][f"{k}1"] == n.iloc[0][f"{k}2"]
        if homo:
            state = "homozygous (not informative)"
        elif np.isnan(tcn):
            state = "no segment"
        elif not np.isnan(lcn) and lcn == 0 and tcn > 0:
            state = "LOH"
        elif not np.isnan(lcn) and (tcn - lcn) > lcn:
            state = "allelic imbalance"
        else:
            state = "balanced"
        rows.append({"sample_id": sid, "patient": sid[:7], "specimen": sid[-1], "gene": gene,
                     "total_cn": tcn, "minor_cn": lcn, "state": state})

out = PROJECT_DIR / "14_immune" / "hla" / "hla_loh.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
t = pd.DataFrame(rows)
t.to_csv(out, sep="\t", index=False)
write_provenance(out, list(D.glob("*/*.segments.tsv"))[:20] + ([types_f] if types_f.exists() else []), {"positions": HLA})
print(t.pivot_table(index="sample_id", columns="gene", values="state", aggfunc="first").to_string() if len(t) else "no FACETS segments")
