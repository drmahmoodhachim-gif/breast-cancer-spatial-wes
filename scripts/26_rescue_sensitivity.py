#!/usr/bin/env python3
"""Keep normal_artifact rescues out of the primary set and score a sensitivity TiN.

Primary calls (PASS) stay in 06_variants/somatic_filtered.tsv.
Rescued calls are written to 06_variants/rescued_normal_artifact.tsv.
Expected normal VAF is the primary shared TiN slope times tumour VAF.
A credible rescue has N_vaf <= 3 x expected, alt>=2 in no other patient's N,
and gnomAD AF < 1e-4 (missing AF counts as not common).
The branch test on primary + credible rescues is written beside the main TiN table.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

OUT = PROJECT_DIR / "06_variants"
VAR = OUT / "somatic_filtered.tsv"
V = pd.read_csv(VAR, sep="\t")
if "call_set" not in V.columns:
    sys.exit("[rescue] call_set column missing; refusing to split the table")

primary = V[V["call_set"].eq("PASS")].copy()
rescued = V[V["call_set"].eq("rescued_normal_artifact")].copy()
if primary.empty or rescued.empty:
    sys.exit(f"[rescue] unexpected split: PASS {len(primary)} rescued {len(rescued)}")


def tin_on(df, vaf_col):
    if df.empty:
        return np.nan
    den = (df["N_depth"] * df[vaf_col]).sum()
    return float(df["N_alt"].sum() / den) if den > 0 else np.nan


slopes = {}
for patient, g in primary.groupby("patient"):
    shared = g[(g["spatial_class"] == "shared") & (g["T_vaf"] > 0)]
    slopes[patient] = tin_on(shared, "T_vaf")

counts = pd.read_csv(PROJECT_DIR / "08_lookback_N" / "allele_counts.tsv", sep="\t")
other_n = counts[(counts["specimen"] == "N") & (counts["alt"] >= 2) &
                 (counts["sample_patient"] != counts["variant_patient"])]
n_other = other_n.groupby("var_id")["sample_id"].nunique()


def gnomad_from_vep(path: Path) -> dict:
    if not path.exists():
        return {}
    fields, out = None, {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = line.split("Format: ")[1].split('"')[0].split("|")
                continue
            if line.startswith("#") or fields is None:
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, info = f[0], int(f[1]), f[3], f[4].split(","), f[7]
            csq = ""
            for kv in info.split(";"):
                if kv.startswith("CSQ="):
                    csq = kv[4:]
            idx = fields.index("gnomADe_AF") if "gnomADe_AF" in fields else None
            if idx is None or not csq:
                continue
            for alt in alts:
                af = ""
                for block in csq.split(","):
                    parts = block.split("|")
                    if len(parts) <= idx:
                        continue
                    picked = parts[fields.index("CANONICAL")] if "CANONICAL" in fields else ""
                    if parts[idx] and (picked == "YES" or not af):
                        af = parts[idx]
                out[(chrom, pos, ref, alt)] = af
    return out


af_map = {}
for patient in rescued["patient"].unique():
    af_map.update(gnomad_from_vep(OUT / "vep_hgvs" / f"{patient}.rescued.vep.vcf.gz"))

rows = []
for r in rescued.itertuples(index=False):
    slope = slopes.get(r.patient, np.nan)
    t_vaf = r.T_vaf if pd.notna(r.T_vaf) else np.nan
    expected = slope * t_vaf if pd.notna(slope) and pd.notna(t_vaf) else np.nan
    n_vaf = r.N_vaf if pd.notna(r.N_vaf) else np.nan
    ratio = n_vaf / expected if pd.notna(expected) and expected > 0 else np.nan
    raw_af = af_map.get((str(r.chrom), int(r.pos), r.ref, r.alt), "")
    af = float(raw_af) if raw_af not in ("", ".") else np.nan
    others = int(n_other.get(r.var_id, 0))
    af_ok = pd.isna(af) or af < 1e-4
    credible = bool(pd.notna(n_vaf) and pd.notna(expected) and n_vaf <= 3 * expected and others == 0 and af_ok)
    rows.append({
        "patient": r.patient, "var_id": r.var_id, "chrom": r.chrom, "pos": r.pos,
        "ref": r.ref, "alt": r.alt, "gene": r.gene, "hgvsp_short": getattr(r, "hgvsp_short", ""),
        "T_vaf": t_vaf, "N_vaf": n_vaf, "expected_N_vaf": expected, "ratio": ratio,
        "n_other_patients_N_with_alt>=2": others, "gnomAD_AF": af,
        "spatial_class": r.spatial_class, "TiN_slope": slope, "credible_rescue": credible,
    })

tab = pd.DataFrame(rows)
dest = OUT / "rescued_normal_artifact.tsv"
tab.to_csv(dest, sep="\t", index=False)
write_provenance(dest, [VAR, PROJECT_DIR / "08_lookback_N" / "allele_counts.tsv"],
                 {"credible_rescue": "N_vaf <= 3 * TiN_slope * T_vaf; no other N with alt>=2; gnomAD AF < 1e-4 or missing"})

primary.to_csv(VAR, sep="\t", index=False)
print(f"[rescue] primary {len(primary)} written back to {VAR.name}")
summary = []
for patient, g in tab.groupby("patient"):
    summary.append({
        "patient": patient,
        "total": len(g),
        "credible_rescue": int(g["credible_rescue"].sum()),
        "T_enriched": int((g["spatial_class"] == "T_private").sum()),
        "J_enriched": int((g["spatial_class"] == "J_private").sum()),
        "credible_T_enriched": int(((g["credible_rescue"]) & (g["spatial_class"] == "T_private")).sum()),
        "credible_J_enriched": int(((g["credible_rescue"]) & (g["spatial_class"] == "J_private")).sum()),
    })
print(pd.DataFrame(summary).to_string(index=False))

# Sensitivity branch test: primary PASS plus credible rescues. Not used by the main figure.
rng = np.random.default_rng(1234)
N_BOOT = 1000
cred = rescued.merge(tab.loc[tab["credible_rescue"], ["var_id"]], on="var_id")
sens = pd.concat([primary, cred], ignore_index=True)


def pct_ci(xs):
    xs = [x for x in xs if np.isfinite(x)]
    if not xs:
        return np.nan, np.nan
    lo, hi = np.percentile(xs, [2.5, 97.5])
    return float(lo), float(hi)


def depletion(shared, private, vaf_col):
    if len(shared) < 5 or private.empty:
        return np.nan
    n_ge = 0
    ns, nq = len(shared), len(private)
    for _ in range(N_BOOT):
        sb = shared.iloc[rng.integers(0, ns, size=ns)]
        qb = private.iloc[rng.integers(0, nq, size=nq)]
        ts, tq = tin_on(sb, vaf_col), tin_on(qb, vaf_col)
        if np.isfinite(ts) and np.isfinite(tq) and ts > 0 and tq >= ts:
            n_ge += 1
    return (1 + n_ge) / (N_BOOT + 1)


sens_rows = []
for patient, g in sens.groupby("patient"):
    shared = g[(g["spatial_class"] == "shared") & (g["T_vaf"] > 0)].reset_index(drop=True)
    sj = shared[shared["J_vaf"] > 0].reset_index(drop=True)
    rec = {"patient": patient, "n_shared": len(shared), "TiN_shared": tin_on(shared, "T_vaf"),
           "n_variants": len(g), "n_credible_added": int((g["var_id"].isin(set(cred["var_id"]))).sum())}
    boots = [tin_on(shared.iloc[rng.integers(0, len(shared), size=len(shared))], "T_vaf")
             for _ in range(N_BOOT)] if len(shared) >= 5 else []
    rec["ci_lo"], rec["ci_hi"] = pct_ci(boots)
    for cls, col, trunk in (("T_private", "T_vaf", shared), ("J_private", "J_vaf", sj)):
        q = g[g["spatial_class"] == cls].reset_index(drop=True)
        rec[f"{cls}_n"] = len(q)
        rec[f"{cls}_p_depleted"] = depletion(trunk, q, col)
    sens_rows.append(rec)

sens_path = PROJECT_DIR / "08_lookback_N" / "tumour_in_normal_sensitivity.tsv"
pd.DataFrame(sens_rows).to_csv(sens_path, sep="\t", index=False)
write_provenance(sens_path, [dest, VAR], {"set": "primary PASS plus credible_rescue", "bootstrap": N_BOOT})
print(pd.DataFrame(sens_rows).round(4).to_string(index=False))
