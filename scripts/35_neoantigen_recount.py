#!/usr/bin/env python3
"""Recount candidate neoantigens from existing pVACseq tables. Does not re-run pVACseq.

Keeps peptides from the MANE Select transcript, otherwise the canonical transcript,
otherwise Transcript Support Level 1. Variant type must be missense, frameshift, or
an in-frame indel. One row per patient and variant, at the lowest mutant IC50.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PATIENTS, PROJECT_DIR, write_provenance  # noqa: E402

I = PROJECT_DIR / "14_immune"
ALLOWED = {"missense", "FS", "inframe_ins", "inframe_del"}
PROT_RE = "missense_variant|frameshift_variant|inframe_insertion|inframe_deletion"


def truthy(series):
    return series.astype(str).str.lower().isin(["true", "yes"])


parts = []
for patient in PATIENTS:
    hits = list((I / "neoantigens" / patient).glob("pvacseq/MHC_Class_I/*.all_epitopes.tsv"))
    if not hits:
        continue
    e = pd.read_csv(hits[0], sep="\t", low_memory=False)
    e["patient"] = patient
    parts.append(e)
E = pd.concat(parts, ignore_index=True)
E["mane"] = truthy(E["MANE Select"])
E["canon"] = truthy(E["Canonical"])
E["tsl1"] = pd.to_numeric(E["Transcript Support Level"], errors="coerce").eq(1)
E["vtype_ok"] = E["Variant Type"].isin(ALLOWED)
E["mt"] = pd.to_numeric(E["MHCflurry MT IC50 Score"], errors="coerce")
E["wt"] = pd.to_numeric(E["MHCflurry WT IC50 Score"], errors="coerce")
E["pct"] = pd.to_numeric(E["MHCflurry MT Percentile"], errors="coerce")
key = ["patient", "Chromosome", "Start", "Reference", "Variant"]

kept = []
for _, g in E.groupby(key, sort=False):
    ok = g[g["vtype_ok"]]
    if ok["mane"].any():
        use = ok[ok["mane"]]
        rule = "MANE"
    elif ok["canon"].any():
        use = ok[ok["canon"]]
        rule = "canonical"
    elif ok["tsl1"].any():
        use = ok[ok["tsl1"]]
        rule = "TSL1"
    else:
        continue
    row = use.sort_values("mt", na_position="last").iloc[0].copy()
    row["transcript_rule"] = rule
    kept.append(row)
K = pd.DataFrame(kept)

V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
V["var_key_start"] = V["patient"] + "|" + V["chrom"].astype(str) + "|" + V["pos"].astype(str) + "|" + V["ref"].astype(str) + "|" + V["alt"].astype(str)
K["var_key_start"] = K["patient"] + "|" + K["Chromosome"].astype(str) + "|" + K["Start"].astype(int).astype(str) + "|" + K["Reference"].astype(str) + "|" + K["Variant"].astype(str)
K["var_key_stop"] = K["patient"] + "|" + K["Chromosome"].astype(str) + "|" + K["Stop"].astype(int).astype(str) + "|" + K["Reference"].astype(str) + "|" + K["Variant"].astype(str)
ann_cols = ["var_id", "patient", "chrom", "pos", "ref", "alt", "gene", "consequence", "hgvsp", "is_coding", "spatial_class"]
A = V[ann_cols].copy()
A["pos"] = A["pos"].astype(int)
K["Start"] = K["Start"].astype(int)
K["Stop"] = K["Stop"].astype(int)
hits = []
for rec in K.itertuples(index=False):
    cand = A[(A["patient"] == rec.patient) & (A["chrom"] == rec.Chromosome) & (A["ref"] == rec.Reference) & (A["alt"] == rec.Variant)]
    cand = cand[cand["pos"].isin([rec.Start, rec.Stop, rec.Start + 1])]
    if len(cand) != 1:
        raise SystemExit(f"[neoag] {rec.patient} {rec.Chromosome}:{rec.Start} {rec.Reference}>{rec.Variant} matched {len(cand)}")
    hits.append(cand.iloc[0]["var_id"])
K["var_id"] = hits
M = K.merge(A.drop(columns=["patient"]), on="var_id", how="left")

lin = pd.read_csv(PROJECT_DIR / "09_ccf_clonality" / "DH00023_lineages.tsv", sep="\t")
lin = lin[lin["record"] == "variant"][["var_id", "lineage"]].drop_duplicates("var_id")
M = M.merge(lin, on="var_id", how="left")


def branch(r):
    if r["patient"] == "DH00023" and pd.notna(r["lineage"]):
        return {"trunk": "trunk", "J-lineage": "J", "T-lineage": "T"}.get(r["lineage"], "unassigned")
    return {"shared": "trunk", "J_private": "J", "T_private": "T"}.get(r["spatial_class"], "unassigned")


M["branch"] = M.apply(branch, axis=1)
M["mt_lt_500"] = M["mt"] < 500
M["mt_lt_50"] = M["mt"] < 50
M["pct_lt_2"] = M["pct"] < 2
M["wt_also_lt_500"] = M["wt"] < 500
M["differential"] = M["mt_lt_500"] & ~M["wt_also_lt_500"].fillna(False)

# The 37 binders whose stored consequence is not a protein-changing SO term.
old = pd.read_csv(I / "neoantigen_variants.tsv", sep="\t")
old = old.merge(V[["var_id", "consequence", "gene", "hgvsp"]], on="var_id", suffixes=("", "_s"))
prot_old = old["consequence"].fillna("").str.contains(PROT_RE)
bad_ids = set(old.loc[old["binder"] & ~prot_old, "var_id"])
lookup = {
    (r.patient, r.chrom, int(r.pos), r.ref, r.alt): r.var_id
    for r in A.itertuples(index=False)
}

def ep_var_id(r):
    for pos in (int(r.Start), int(r.Stop), int(r.Start) + 1):
        hit = lookup.get((r.patient, r.Chromosome, pos, r.Reference, r.Variant))
        if hit:
            return hit
    return np.nan

E["var_id"] = E.apply(ep_var_id, axis=1)
bad_rows = E[E["var_id"].isin(bad_ids) & (E["mt"] < 500)].copy()
report = []
for vid, g in bad_rows.groupby("var_id"):
    report.append({
        "var_id": vid,
        "n_peptides_lt_500": len(g),
        "n_noncanonical": int((~g["mane"]).sum()),
        "any_mane": bool(g["mane"].any()),
        "genes": ",".join(sorted(set(g["Gene Name"].dropna().astype(str)))),
        "variant_types": ",".join(sorted(set(g["Variant Type"].dropna().astype(str)))),
        "transcripts": ",".join(sorted(set(g["Transcript"].dropna().astype(str)))),
    })
R = pd.DataFrame(report)
# variants among the 37 with no IC50<500 peptide row still need a line
seen = set(R["var_id"]) if len(R) else set()
for vid in sorted(bad_ids - seen):
    g = E[E["var_id"] == vid]
    R = pd.concat([R, pd.DataFrame([{
        "var_id": vid,
        "n_peptides_lt_500": 0,
        "n_noncanonical": int((~g["mane"]).sum()) if len(g) else 0,
        "any_mane": bool(g["mane"].any()) if len(g) else False,
        "genes": ",".join(sorted(set(g["Gene Name"].dropna().astype(str)))) if len(g) else "",
        "variant_types": ",".join(sorted(set(g["Variant Type"].dropna().astype(str)))) if len(g) else "",
        "transcripts": ",".join(sorted(set(g["Transcript"].dropna().astype(str)))) if len(g) else "",
    }])], ignore_index=True)
som = V.set_index("var_id")
if not len(R):
    R = pd.DataFrame(columns=["var_id", "n_peptides_lt_500", "n_noncanonical", "any_mane", "genes", "variant_types", "transcripts"])
R["somatic_gene"] = R["var_id"].map(som["gene"])
R["somatic_consequence"] = R["var_id"].map(som["consequence"])
R["somatic_hgvsp"] = R["var_id"].map(som["hgvsp"])
R["from_noncanonical_transcript"] = ~R["any_mane"]
n_noncanon = int(R["from_noncanonical_transcript"].sum())
print(f"[37] stored non-coding binders {len(bad_ids)}; from a non-canonical transcript: {n_noncanon}")

# Denominator: protein-altering on the stored SO term, plus MANE peptide variants the SO term missed.
V["protein_so"] = V["consequence"].fillna("").str.contains(PROT_RE)
mane_ids = set(M["var_id"])
V["protein_altering"] = V["protein_so"] | V["var_id"].isin(mane_ids)
V = V.merge(M[["var_id", "branch"]], on="var_id", how="left")
need = V["protein_altering"] & V["branch"].isna()
def branch_v(r):
    if r["patient"] == "DH00023":
        hit = lin.loc[lin["var_id"] == r["var_id"], "lineage"]
        if len(hit):
            return {"trunk": "trunk", "J-lineage": "J", "T-lineage": "T"}.get(hit.iloc[0], "unassigned")
    return {"shared": "trunk", "J_private": "J", "T_private": "T"}.get(r["spatial_class"], "unassigned")
if need.any():
    V.loc[need, "branch"] = V.loc[need].apply(branch_v, axis=1)

rows = []
for patient in PATIENTS:
    pv = V[V["patient"] == patient]
    pm = M[M["patient"] == patient]
    b500 = pm[pm["mt_lt_500"]]
    diff = pm[pm["differential"]]
    rec = {
        "patient": patient,
        "protein_altering": int(pv["protein_altering"].sum()),
        "is_coding": int(pv["is_coding"].isin([True, "True"]).sum()),
        "mt_ic50_lt_500": int(pm["mt_lt_500"].sum()),
        "mt_ic50_lt_50": int(pm["mt_lt_50"].sum()),
        "percentile_lt_2": int(pm["pct_lt_2"].sum()),
        "wt_also_lt_500": int((pm["mt_lt_500"] & pm["wt_also_lt_500"]).sum()),
        "differential_lt_500": int(pm["differential"].sum()),
        "trunk": int((diff["branch"] == "trunk").sum()),
        "J_branch": int((diff["branch"] == "J").sum()),
        "T_branch": int((diff["branch"] == "T").sum()),
        "unassigned": int((diff["branch"] == "unassigned").sum()),
        "lt500_trunk": int((b500["branch"] == "trunk").sum()),
        "lt500_J_branch": int((b500["branch"] == "J").sum()),
        "lt500_T_branch": int((b500["branch"] == "T").sum()),
        "lt500_unassigned": int((b500["branch"] == "unassigned").sum()),
    }
    for br in ("trunk", "J", "T", "unassigned"):
        n_b = int((b500["branch"] == br).sum())
        n_d = int(((pv["protein_altering"]) & (pv["branch"] == br)).sum())
        if n_b > n_d:
            raise SystemExit(f"[neoag] {patient} {br}: binders {n_b} > protein-altering {n_d}")
    rows.append(rec)
C = pd.DataFrame(rows)
C["label"] = "candidate neoantigens (no RNA expression data)"
print(C.to_string(index=False))
print("[pct] min/max", float(M["pct"].min()), float(M["pct"].max()))

keep_cols = [
    "patient", "var_id", "Gene Name", "Transcript", "transcript_rule", "Variant Type",
    "HGVSp", "mt", "wt", "pct", "mt_lt_500", "mt_lt_50", "pct_lt_2", "wt_also_lt_500",
    "differential", "spatial_class", "lineage", "branch", "consequence",
]
out_var = M[keep_cols].rename(columns={
    "Gene Name": "gene_mane", "mt": "best_mt_ic50", "wt": "wt_ic50", "pct": "mt_percentile",
})
out_var = out_var.sort_values(["patient", "best_mt_ic50"])
var_path = I / "candidate_neoantigens.tsv"
count_path = I / "candidate_neoantigen_counts.tsv"
bad_path = I / "candidate_neoantigen_noncoding_source.tsv"
out_var.to_csv(var_path, sep="\t", index=False)
C.to_csv(count_path, sep="\t", index=False)
R.to_csv(bad_path, sep="\t", index=False)

S = pd.read_csv(I / "immune_summary.tsv", sep="\t")
for col in [c for c in S.columns if c.startswith("neoantigenic") or c.startswith("patient_") or c.startswith("candidate_") or c in {
    "is_coding", "mt_ic50_lt_500", "mt_ic50_lt_50", "percentile_lt_2", "wt_also_lt_500",
    "differential_lt_500", "lt500_trunk", "lt500_J_branch", "lt500_T_branch", "lt500_unassigned",
}]:
    S.drop(columns=col, inplace=True)
add = C.drop(columns=["label"]).rename(columns={
    "protein_altering": "patient_protein_altering",
    "differential_lt_500": "candidate_neoantigens",
    "trunk": "candidate_trunk",
    "J_branch": "candidate_J_branch",
    "T_branch": "candidate_T_branch",
    "unassigned": "candidate_unassigned",
})
# Store the patient count on the tumour row only, so a sum across specimens does not triple it.
S = S.merge(add, on="patient", how="left")
patient_cols = [c for c in add.columns if c != "patient"]
S.loc[S["specimen"] != "T", patient_cols] = np.nan
S["candidate_neoantigen_note"] = np.where(
    S["specimen"] == "T",
    "candidate neoantigens (no RNA expression data); mutant IC50 < 500 nM and wild-type IC50 not < 500 nM; per patient, this row only",
    "",
)
sum_path = I / "immune_summary.tsv"
S.to_csv(sum_path, sep="\t", index=False)
inputs = list((I / "neoantigens").glob("DH000*/pvacseq/MHC_Class_I/*.all_epitopes.tsv"))
inputs.append(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv")
write_provenance(count_path, inputs, {
    "rule": "MANE Select, else canonical, else TSL 1; missense/FS/inframe; best MT IC50 per variant",
    "noncanonical_among_37": n_noncanon,
})
print("[neoag] wrote", count_path.name)
