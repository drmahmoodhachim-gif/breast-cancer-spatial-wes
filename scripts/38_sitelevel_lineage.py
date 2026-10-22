#!/usr/bin/env python3
"""Post hoc site-level test of DH00023 lineages in adjacent normal.

Counts sites with at least one alternate read, and the pooled reads, for
J-lineage versus T-lineage. This is not the beta-binomial branch test.
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

NOTE = "post hoc. Not a replacement for the beta-binomial branch test (D042)."
lin = pd.read_csv(PROJECT_DIR / "09_ccf_clonality" / "DH00023_lineages.tsv", sep="\t")
lin = lin[lin["record"] == "variant"][["var_id", "lineage", "J_vaf", "T_vaf"]].drop_duplicates("var_id")
var = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
df = lin.merge(var[["var_id", "N_alt", "N_depth"]], on="var_id", how="inner")
df = df[df["lineage"].isin(["J-lineage", "T-lineage"])].copy()
df["N_alt"] = pd.to_numeric(df["N_alt"], errors="coerce").fillna(0)
df["N_depth"] = pd.to_numeric(df["N_depth"], errors="coerce").fillna(0)
df["has_alt"] = df["N_alt"] >= 1

rows = []
site = {}
read = {}
for lineage in ("J-lineage", "T-lineage"):
    g = df[df["lineage"] == lineage]
    source = g["J_vaf"] if lineage == "J-lineage" else g["T_vaf"]
    site[lineage] = (int(g["has_alt"].sum()), int(len(g)))
    read[lineage] = (int(g["N_alt"].sum()), int(g["N_depth"].sum()))
    rows.append({
        "test": "source_vaf",
        "lineage": lineage,
        "n_sites": int(len(g)),
        "n_sites_with_N_alt": int(g["has_alt"].sum()),
        "N_alt_reads": int(g["N_alt"].sum()),
        "N_depth_reads": int(g["N_depth"].sum()),
        "source_vaf_median": float(pd.to_numeric(source, errors="coerce").median()),
        "odds_ratio": "",
        "p_fisher": "",
        "note": NOTE,
    })

for name, table in (("sites", site), ("pooled_reads", read)):
    # rows: J then T. columns: event, non-event.
    j_hit, j_n = table["J-lineage"]
    t_hit, t_n = table["T-lineage"]
    oddsratio, pvalue = fisher_exact([[j_hit, j_n - j_hit], [t_hit, t_n - t_hit]])
    rows.append({
        "test": name,
        "lineage": "J-lineage vs T-lineage",
        "n_sites": "",
        "n_sites_with_N_alt": f"{j_hit}/{j_n} vs {t_hit}/{t_n}",
        "N_alt_reads": j_hit if name == "pooled_reads" else "",
        "N_depth_reads": j_n if name == "pooled_reads" else "",
        "source_vaf_median": "",
        "odds_ratio": oddsratio,
        "p_fisher": pvalue,
        "note": NOTE,
    })

out = pd.DataFrame(rows)
path = PROJECT_DIR / "08_lookback_N" / "branch_test_sitelevel.tsv"
out.to_csv(path, sep="\t", index=False)
write_provenance(path, [
    PROJECT_DIR / "09_ccf_clonality" / "DH00023_lineages.tsv",
    PROJECT_DIR / "06_variants" / "somatic_filtered.tsv",
], {"note": NOTE})
print(out.to_string(index=False))
