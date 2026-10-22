#!/usr/bin/env python3
"""Copy the audited numbers and the figures that match them into 13_paper/.

Reads existing tables. Does not rerun alignment, variant calling, FACETS, HLA, or dilution.
"""
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_paper"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

V = pd.read_csv(ROOT / "06_variants" / "somatic_filtered.tsv", sep="\t")
counts = (
    V.groupby(["patient", "spatial_class"]).size().unstack(fill_value=0).reset_index()
)
for col in ("shared", "T_private", "J_private", "T_unresolved", "J_unresolved"):
    if col not in counts.columns:
        counts[col] = 0
counts["n"] = counts[["shared", "T_private", "J_private", "T_unresolved", "J_unresolved"]].sum(axis=1)
counts = counts[["patient", "shared", "T_private", "J_private", "T_unresolved", "J_unresolved", "n"]]
counts.to_csv(TAB / "table_variant_counts.tsv", sep="\t", index=False)

drv = pd.read_csv(ROOT / "12_tables" / "main" / "table2_driver_alterations.tsv", sep="\t")
keep = ["patient", "gene", "hgvsp", "consequence", "chrom", "pos", "ref", "alt",
        "spatial_class", "J_vaf", "T_vaf", "J_ccf", "T_ccf"]
drv = drv[keep].copy()
drv["oncoprint_row"] = drv["gene"].where(drv["gene"] != "BRCA2", "BRCA2 VUS")
for col in ("J_vaf", "T_vaf"):
    drv[col] = drv[col].round(3)
for col in ("J_ccf", "T_ccf"):
    drv[col] = drv[col].round(2)
drv.to_csv(TAB / "table_oncoprint_variants.tsv", sep="\t", index=False)

pur = pd.read_csv(ROOT / "07_cnv_purity" / "purity_reconciled.tsv", sep="\t")
pur["used_purity"] = pur["used_purity"].round(3)
pur.to_csv(TAB / "table_purity.tsv", sep="\t", index=False)

tin = pd.read_csv(ROOT / "08_lookback_N" / "tumour_in_normal_estimates.tsv", sep="\t")
tin_out = pd.DataFrame({
    "patient": tin["patient"],
    "n_shared": tin["n_shared"],
    "TiN": tin["TiN_shared"].round(4),
    "TiN_ci_lo": tin["ci_lo"].round(4),
    "TiN_ci_hi": tin["ci_hi"].round(4),
    "T_private_n": tin["T_private_n"],
    "T_private_N_alt_expected": tin["T_private_N_alt_expected"].round(1),
    "T_private_N_alt_observed": tin["T_private_N_alt_observed"],
    "T_private_test": tin["T_private_p_label"],
    "J_private_n": tin["J_private_n"],
    "J_private_N_alt_expected": tin["J_private_N_alt_expected"].round(1),
    "J_private_N_alt_observed": tin["J_private_N_alt_observed"],
    "J_private_test": tin["J_private_p_label"],
})
tin_out.to_csv(TAB / "table_tumour_in_normal.tsv", sep="\t", index=False)

dil = pd.read_csv(ROOT / "08_lookback_N" / "dilution_test.tsv", sep="\t")
slope_cols = ["baseline_slope", "expected_slope", "observed_slope", "expected_increment",
              "recovered", "ci_lo", "ci_hi", "fit_slope", "fit_intercept", "fit_r2"]
for col in slope_cols:
    dil[col] = dil[col].round(4)
dil[["patient", "spiked_fraction", "n_shared"] + slope_cols].to_csv(
    TAB / "table_dilution.tsv", sep="\t", index=False)

hla = pd.read_csv(ROOT / "14_immune" / "hla" / "hla_types.tsv", sep="\t")
hla["Reads"] = hla["Reads"].round(0).astype(int)
hla.drop(columns=["Objective"], errors="ignore").to_csv(TAB / "table_hla.tsv", sep="\t", index=False)

cn = pd.read_csv(ROOT / "07_cnv_purity" / "cn_events.tsv", sep="\t")
cn = cn[cn["event"].isin(["FGFR1 amp", "16q loss"])].copy()
cn["tcn"] = cn["tcn"].round(1)
cn["lcn"] = cn["lcn"].round(1)
cn["ploidy"] = cn["ploidy"].round(2)
cn["segment_Mb"] = (cn["segment_size"] / 1e6).round(2)
cn["fraction"] = cn["fraction"].round(2)
cn[["sample", "event", "called", "tcn", "lcn", "ploidy", "segment_Mb", "fraction"]].to_csv(
    TAB / "table_fgfr1_16q.tsv", sep="\t", index=False)

lin = pd.read_csv(ROOT / "09_ccf_clonality" / "DH00023_lineages.tsv", sep="\t")
frac_path = ROOT / "09_ccf_clonality" / "dh00023_lineage_fractions_d065.tsv"
frac = pd.read_csv(frac_path, sep="\t")
for col in ("median_vaf", "vaf_ratio_to_trunk", "median_uncapped_ccf"):
    frac[col] = frac[col].round(3)
frac.to_csv(TAB / "table_dh00023_lineage_fractions.tsv", sep="\t", index=False)
n_by = lin[lin["record"] == "variant"].groupby("lineage").size()
n_by.to_csv(TAB / "table_dh00023_lineage_counts.tsv", sep="\t", header=["n"])

tmb = pd.read_csv(ROOT / "12_tables" / "main" / "tmb.tsv", sep="\t")
tmb.to_csv(TAB / "table_tmb.tsv", sep="\t", index=False)

ped = pd.read_csv(ROOT / "16_ancestry" / "peddy.peddy.ped", sep="\t")
ped = ped[["sample_id", "ancestry-prediction", "PC1", "PC2"]].copy()
ped["PC1"] = ped["PC1"].round(3)
ped["PC2"] = ped["PC2"].round(3)
ped.to_csv(TAB / "table_ancestry.tsv", sep="\t", index=False)

ident = pd.read_csv(ROOT / "04_qc" / "sample_identity.tsv", sep="\t")
same = ident["same_patient"].astype(str).str.lower().isin(["true", "1", "yes"])
rows = []
for label, mask in (("within_patient", same), ("between_patient", ~same)):
    sub = ident.loc[mask, "correlation"]
    rows.append({
        "comparison": label,
        "n_pairs": int(mask.sum()),
        "min_correlation": round(float(sub.min()), 4),
        "max_correlation": round(float(sub.max()), 4),
        "scale": "NGSCheckMate SNP-VAF correlation",
    })
pd.DataFrame(rows).to_csv(TAB / "table_identity.tsv", sep="\t", index=False)

locked_figs = [
    "main/Fig1_design_QC",
    "main/Fig2_oncoprint",
    "main/Fig3_TJ_sharing",
    "main/Fig3_concordance",
    "main/Fig4_CCF_T_vs_J",
    "main/Fig5_CNV_purity",
    "main/Fig5_DH00023",
    "main/Fig6_lookback_N",
    "main/Fig7_tumour_in_normal",
    "main/Fig9_clonal_sharing_per_patient_PRELIMINARY",
    "supplementary/FigS2_TMB",
    "supplementary/FigS_ancestry",
    "supplementary/FigS_dilution",
    "supplementary/FigS_immune",
    "supplementary/FigS_rescued",
]
copied = []
for rel in locked_figs:
    for ext in (".pdf", ".png"):
        src = ROOT / "11_figures" / f"{rel}{ext}"
        if src.exists():
            dst = FIG / src.name
            shutil.copy2(src, dst)
            copied.append(src.name)

print(f"[lock] {len(V)} primary variants")
print(counts.to_string(index=False))
print("[lock] identity")
print(pd.DataFrame(rows).to_string(index=False))
print(f"[lock] figures copied: {len(copied)}")
for name in copied:
    print(" ", name)
