#!/usr/bin/env python3
"""Reports for the rescued look-back and AKT1, plus DH00023 lineage fractions and CN events."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

ROOT = PROJECT_DIR
V = pd.read_csv(ROOT / "06_variants" / "somatic_filtered.tsv", sep="\t")
R = pd.read_csv(ROOT / "06_variants" / "rescued_normal_artifact.tsv", sep="\t")
L = pd.read_csv(ROOT / "08_lookback_N" / "lookback_results.tsv", sep="\t")
N = L[(L["specimen"] == "N") & L["var_id"].isin(R["var_id"])].copy()
hi = N[(N["evidence_class"] == "confident") & (N["vaf"] >= 0.05)].copy()
hi = hi.merge(R[["var_id", "gnomAD_AF", "n_other_patients_N_with_alt>=2", "credible_rescue", "hgvsp_short"]],
              on="var_id", how="left")
cols = ["variant_patient", "gene", "hgvsp_short", "vaf", "T_vaf", "J_vaf", "gnomAD_AF",
        "n_other_patients_N_with_alt>=2", "credible_rescue", "var_id"]
out45 = ROOT / "08_lookback_N" / "rescued_n_vaf_ge5.tsv"
hi[cols].sort_values(["variant_patient", "vaf"], ascending=[True, False]).to_csv(out45, sep="\t", index=False)
write_provenance(out45, [ROOT / "06_variants" / "rescued_normal_artifact.tsv",
                         ROOT / "08_lookback_N" / "lookback_results.tsv"],
                 {"rule": "rescued sites with look-back class confident and N VAF >= 0.05"})
print(f"[report] confident rescued N VAF>=5%: {len(hi)}")

A = pd.read_csv(ROOT / "06_variants" / "somatic_all.tsv", sep="\t")
akt = A[A["gene"] == "AKT1"].copy()
primary_ids = set(V["var_id"])
rescued_ids = set(R["var_id"])
akt["set"] = np.where(akt["var_id"].isin(primary_ids), "primary",
                      np.where(akt["var_id"].isin(rescued_ids), "rescued", "excluded"))
akt_cols = ["patient", "set", "filter", "hgvsp", "hgvsp_short", "N_alt", "N_depth", "J_alt", "J_depth",
            "T_alt", "T_depth", "N_vaf", "var_id"]
akt_out = ROOT / "08_lookback_N" / "akt1_sites.tsv"
akt[akt_cols].to_csv(akt_out, sep="\t", index=False)
print(akt[akt_cols].to_string(index=False))

# ---- DH00023 lineages (D043) -------------------------------------------------
g = V[V["patient"] == "DH00023"].copy()
both = (g["T_vaf"] > 0) & (g["J_vaf"] > 0)
ratio = np.where(both, g["T_vaf"] / g["J_vaf"], np.nan)
g["vaf_ratio"] = ratio
g["lineage"] = "other"
g.loc[(g["spatial_class"] == "shared") & g["vaf_ratio"].between(0.6, 1.6), "lineage"] = "trunk"
g.loc[(g["spatial_class"] == "shared") & (g["vaf_ratio"] < 0.6), "lineage"] = "J-lineage"
g.loc[g["spatial_class"] == "J_private", "lineage"] = "J-lineage"
g.loc[g["spatial_class"] == "T_private", "lineage"] = "T-lineage"

pur = pd.read_csv(ROOT / "07_cnv_purity" / "purity_reconciled.tsv", sep="\t")
tin = pd.read_csv(ROOT / "08_lookback_N" / "tumour_in_normal_estimates.tsv", sep="\t")
t_pur = float(pur.loc[pur["sample"] == "DH00023_T", "used_purity"].iloc[0])
j_pur = float(pur.loc[pur["sample"] == "DH00023_J", "used_purity"].iloc[0])
slope = float(tin.loc[tin["patient"] == "DH00023", "TiN_shared"].iloc[0])
purity = {"T": t_pur, "J": j_pur, "N": slope * t_pur}
vaf_col = {"T": "T_vaf", "J": "J_vaf", "N": "N_vaf"}

trunk = g[g["lineage"] == "trunk"]
rows = []
for region, col in vaf_col.items():
    med_trunk = float(trunk[col].median()) if len(trunk) else np.nan
    rel = {}
    for lineage in ("trunk", "J-lineage", "T-lineage"):
        sub = g[g["lineage"] == lineage]
        med = float(sub[col].median()) if len(sub) else np.nan
        ccf = med / med_trunk if np.isfinite(med) and np.isfinite(med_trunk) and med_trunk > 0 else np.nan
        if lineage == "trunk":
            ccf = 1.0 if len(sub) else np.nan
        ccf = float(min(1.0, ccf)) if np.isfinite(ccf) else np.nan
        rel[lineage] = 0.0 if not np.isfinite(ccf) else ccf
        rows.append({"record": "lineage", "region": region, "lineage": lineage, "n_variants": len(sub),
                     "median_vaf": med, "ccf_vs_trunk": ccf, "purity": purity[region]})
    j_rel, t_rel = rel["J-lineage"], rel["T-lineage"]
    if j_rel + t_rel > 1:
        scale = j_rel + t_rel
        j_rel, t_rel = j_rel / scale, t_rel / scale
        trunk_only = 0.0
    else:
        trunk_only = max(0.0, 1.0 - j_rel - t_rel)
    cancer = purity[region]
    fractions = {
        "normal": 1.0 - cancer,
        "trunk-only": trunk_only * cancer,
        "J-lineage": j_rel * cancer,
        "T-lineage": t_rel * cancer,
    }
    for name, frac in fractions.items():
        rows.append({"record": "fraction", "region": region, "lineage": name, "n_variants": np.nan,
                     "median_vaf": np.nan, "ccf_vs_trunk": np.nan, "purity": cancer,
                     "sample_fraction": frac})

var_rows = g[["var_id", "gene", "hgvsp_short", "spatial_class", "lineage", "T_vaf", "J_vaf", "N_vaf", "vaf_ratio"]].copy()
var_rows.insert(0, "record", "variant")
var_rows.insert(1, "region", "")
fra = pd.DataFrame(rows)
lin = pd.concat([var_rows, fra], ignore_index=True)
lin_path = ROOT / "09_ccf_clonality" / "DH00023_lineages.tsv"
lin.to_csv(lin_path, sep="\t", index=False)
write_provenance(lin_path, [ROOT / "06_variants" / "somatic_filtered.tsv"],
                 {"decision": "D043", "trunk_ratio": "0.6-1.6", "N_purity": "TiN slope * T purity"})
print("[lineage] counts")
print(g["lineage"].value_counts().to_string())
kmt = g[g["gene"] == "KMT2C"][["gene", "hgvsp_short", "lineage", "T_vaf", "J_vaf", "vaf_ratio"]]
print(kmt.to_string(index=False))

# ---- copy-number events (D044) ----------------------------------------------
ARMS = {
    "1q gain": ("1", 125100001, 248956422, "gain"),
    "8q gain": ("8", 45900001, 145138636, "gain"),
    "16q loss": ("16", 36300001, 90338345, "loss"),
    "17p LOH": ("17", 1, 22800000, "loh"),
}
GENES = {
    "CCND1 amp": ("11", 69641156, 69654474),
    "FGFR1 amp": ("8", 38411139, 38468834),
}
seg_files = list((ROOT / "07_cnv_purity" / "facets").rglob("*segments.tsv"))
seg = pd.concat([pd.read_csv(p, sep="\t") for p in seg_files], ignore_index=True)
seg["chrom"] = seg["chrom"].astype(str).str.replace("^chr", "", regex=True)
pl = pd.read_csv(ROOT / "07_cnv_purity" / "purity_ploidy.tsv", sep="\t")
samples = sorted(seg["sample_id"].unique())


def overlap(g, chrom, start, end):
    m = g[(g["chrom"] == chrom) & (g["end"] >= start) & (g["start"] <= end)].copy()
    if m.empty:
        return m
    m["w"] = (np.minimum(m["end"], end) - np.maximum(m["start"], start)).clip(lower=1)
    return m


def wmedian(values, weights):
    if len(values) == 0 or np.nansum(weights) <= 0:
        return np.nan
    order = np.argsort(values)
    v, w = np.asarray(values, float)[order], np.asarray(weights, float)[order]
    cut = np.cumsum(w) >= (w.sum() / 2)
    return float(v[np.argmax(cut)])


cn_rows = []
for sid in samples:
    g0 = seg[seg["sample_id"] == sid]
    ploidy = float(pl.loc[pl["sample_id"] == sid, "ploidy"].iloc[0])
    base = int(round(ploidy))
    for event, (chrom, start, end) in GENES.items():
        m = overlap(g0, chrom, start, end)
        if m.empty:
            cn_rows.append({"sample": sid, "event": event, "locus": f"chr{chrom}:{start}-{end}",
                            "called": False, "tcn": np.nan, "lcn": np.nan, "ploidy": ploidy,
                            "segment_size": np.nan})
            continue
        # The segment that covers the most of the gene.
        top = m.sort_values("w", ascending=False).iloc[0]
        size = int(top["end"] - top["start"])
        called = bool(pd.notna(top["tcn.em"]) and top["tcn.em"] >= 2 * ploidy and size <= 10_000_000)
        cn_rows.append({"sample": sid, "event": event, "locus": f"chr{chrom}:{start}-{end}",
                        "called": called, "tcn": top["tcn.em"], "lcn": top["lcn.em"],
                        "ploidy": ploidy, "segment_size": size})
    for event, (chrom, start, end, kind) in ARMS.items():
        m = overlap(g0, chrom, start, end)
        arm = end - start
        if m.empty:
            cn_rows.append({"sample": sid, "event": event, "locus": f"chr{chrom}:{start}-{end}",
                            "called": False, "tcn": np.nan, "lcn": np.nan, "ploidy": ploidy,
                            "segment_size": arm})
            continue
        tcn = wmedian(m["tcn.em"], m["w"])
        lcn = wmedian(m["lcn.em"].dropna(), m.loc[m["lcn.em"].notna(), "w"]) if m["lcn.em"].notna().any() else np.nan
        if kind == "gain":
            hit = m["tcn.em"].notna() & (m["tcn.em"] > base)
        elif kind == "loss":
            hit = m["tcn.em"].notna() & (m["tcn.em"] < base)
        else:
            hit = m["lcn.em"].notna() & (m["lcn.em"] <= 0)
        frac = float(m.loc[hit, "w"].sum() / arm) if arm else 0
        cn_rows.append({"sample": sid, "event": event, "locus": f"chr{chrom}:{start}-{end}",
                        "called": frac >= 0.5, "tcn": tcn, "lcn": lcn, "ploidy": ploidy,
                        "segment_size": arm, "fraction": round(frac, 3)})

cn = pd.DataFrame(cn_rows)
cn_path = ROOT / "07_cnv_purity" / "cn_events.tsv"
cn.to_csv(cn_path, sep="\t", index=False)
write_provenance(cn_path, [ROOT / "07_cnv_purity" / "purity_ploidy.tsv"],
                 {"decision": "D044", "amp": "tcn >= 2*ploidy and segment <= 10 Mb",
                  "arm": ">=50% of bases above/below round(ploidy)", "loh": "lcn 0 on >=50% of the arm"})
focus = cn[cn["event"].isin(["FGFR1 amp", "16q loss"])][["sample", "event", "called", "tcn", "lcn", "ploidy", "segment_size"]]
print(focus.to_string(index=False))
