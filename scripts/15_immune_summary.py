#!/usr/bin/env python3
"""Combine the immune analyses into one per-sample table.

Inputs (each optional; missing parts are left empty):
  14_immune/tcell/tcell_fraction.tsv             T cell ExTRECT
  14_immune/hla/hla_types.tsv, hla_loh.tsv       OptiType, FACETS-based HLA LOH
  14_immune/neoantigens/<patient>/pvacseq/MHC_Class_I/*.all_epitopes.tsv
  06_variants/somatic_filtered.tsv                spatial classes
Outputs:
  14_immune/neoantigen_variants.tsv   one row per variant: best mutant IC50, binder (<cutoff), spatial class
  14_immune/immune_summary.tsv        one row per sample: T-cell fraction, neoantigen burden, HLA LOH
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PATIENTS, PROJECT_DIR, load_env, write_provenance  # noqa: E402

env = load_env()
CUT = float(env.get("NEOAG_IC50_CUTOFF", 500))
I = PROJECT_DIR / "14_immune"
V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")

# ---- neoantigens ----------------------------------------------------------------
ep_rows = []
for p in PATIENTS:
    for f in (I / "neoantigens" / p).glob("pvacseq/MHC_Class_I/*.all_epitopes.tsv"):
        e = pd.read_csv(f, sep="\t", low_memory=False)
        ic = next((c for c in ["Median MT IC50 Score", "Best MT IC50 Score", "MHCflurry MT IC50 Score"] if c in e.columns),
                  next((c for c in e.columns if "MT IC50" in c), None))
        if ic is None or e.empty:
            continue
        e = e.rename(columns={ic: "mt_ic50"})
        e["patient"] = p
        ep_rows.append(e[["patient", "Chromosome", "Start", "Stop", "Reference", "Variant", "Gene Name",
                          "HLA Allele", "MT Epitope Seq", "mt_ic50"]])
NE = pd.DataFrame()
if ep_rows:
    E = pd.concat(ep_rows, ignore_index=True)
    E["mt_ic50"] = pd.to_numeric(E["mt_ic50"], errors="coerce")
    best = E.sort_values("mt_ic50").groupby(["patient", "Chromosome", "Start", "Stop"], as_index=False).first()
    v = V[["patient", "chrom", "pos", "var_id", "gene", "spatial_class", "T_present", "J_present"]].copy()
    m1 = best.merge(v, left_on=["patient", "Chromosome", "Stop"], right_on=["patient", "chrom", "pos"], how="inner")
    m2 = best.merge(v, left_on=["patient", "Chromosome", "Start"], right_on=["patient", "chrom", "pos"], how="inner")
    NE = pd.concat([m1, m2]).drop_duplicates("var_id")
    NE["binder"] = NE["mt_ic50"] < CUT
    NE.to_csv(I / "neoantigen_variants.tsv", sep="\t", index=False)

# ---- per-sample summary ---------------------------------------------------------
rows = []
tc_f = I / "tcell" / "tcell_fraction.tsv"
TC = pd.read_csv(tc_f, sep="\t").set_index("sample") if tc_f.exists() else pd.DataFrame()
loh_f = I / "hla" / "hla_loh.tsv"
LOH = pd.read_csv(loh_f, sep="\t") if loh_f.exists() else pd.DataFrame()
for p in PATIENTS:
    for s in ("N", "J", "T"):
        sid = f"{p}_{s}"
        r = {"sample_id": sid, "patient": p, "specimen": s}
        if sid in TC.index:
            for c in ["TCRA.tcell.fraction", "TCRA.tcell.fraction.lwr", "TCRA.tcell.fraction.upr", "qcFit"]:
                r[c] = TC.loc[sid, c] if c in TC.columns else np.nan
            adj = [c for c in TC.columns if "adj" in c.lower() and "fraction" in c.lower()]
            r["tcell_fraction_purity_adjusted"] = TC.loc[sid, adj[0]] if adj else np.nan
        if s != "N":
            pv = V[(V["patient"] == p) & V[f"{s}_present"]]
            r["somatic_variants_present"] = len(pv)
            if len(NE):
                nb = NE[(NE["patient"] == p) & NE[f"{s}_present"] & NE["binder"]]
                r["neoantigenic_variants"] = len(nb)
                r["neoantigenic_shared"] = int((nb["spatial_class"] == "shared").sum())
                r["neoantigenic_private"] = int((nb["spatial_class"] == f"{s}_private").sum())
            if len(LOH):
                l = LOH[LOH["sample_id"] == sid]
                r["HLA_LOH_genes"] = ",".join(l.loc[l["state"] == "LOH", "gene"]) or "none"
        rows.append(r)
S = pd.DataFrame(rows)
out = I / "immune_summary.tsv"
S.to_csv(out, sep="\t", index=False)
write_provenance(out, [p for p in [tc_f, loh_f] if p.exists()], {"ic50_cutoff": CUT})
print(S.to_string(index=False))
