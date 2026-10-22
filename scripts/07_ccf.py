#!/usr/bin/env python3
"""Purity- and copy-number-adjusted cancer cell fraction (CCF) for T and J.

For variant with VAF v in a sample of purity p, local total CN (tumour) C and normal CN 2:
    expected cellular copies   u = v * (p*C + 2*(1-p)) / p
    multiplicity               m = max(1, round(u)) capped at the major allele CN
    CCF                        = u / m           (capped at 1 for display; raw kept)
95% CI from the Wilson interval on VAF, propagated through the same transform.
Clonal call: CCF upper bound >= 0.9 (reported as 'clonal'), else 'subclonal'.

Inputs : 06_variants/somatic_filtered.tsv, 07_cnv_purity/purity_ploidy.tsv,
         07_cnv_purity/facets/*/*.segments.tsv
Output : 09_ccf_clonality/ccf.tsv, 07_cnv_purity/purity_reconciled.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
PP = pd.read_csv(PROJECT_DIR / "07_cnv_purity" / "purity_ploidy.tsv", sep="\t").set_index("sample_id")
segs = pd.concat([pd.read_csv(f, sep="\t") for f in (PROJECT_DIR / "07_cnv_purity" / "facets").glob("*/*.segments.tsv")])
segs["chrom"] = segs["chrom"].astype(str).str.replace("chr", "", regex=False)
SEG = {sid: g for sid, g in segs.groupby("sample_id")}


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def local_cn(sid, chrom, pos):
    g = SEG.get(sid)
    if g is None:
        return np.nan, np.nan
    hit = g[(g["chrom"] == str(chrom).replace("chr", "")) & (g["start"] <= pos) & (g["end"] >= pos)]
    if hit.empty:
        return np.nan, np.nan
    return float(hit["tcn.em"].iloc[0]), float(hit["lcn.em"].iloc[0]) if pd.notna(hit["lcn.em"].iloc[0]) else np.nan


def ccf_from_vaf(v, p, C, major):
    if any(pd.isna(x) for x in (v, p, C)) or p <= 0:
        return np.nan, np.nan
    u = v * (p * C + 2 * (1 - p)) / p
    m = max(1, int(round(u)))
    if not pd.isna(major) and major >= 1:
        m = min(m, int(major))
    return u / m, m


def n_diploid_shared(patient, specimen, sh):
    """Shared mutations on a balanced diploid segment (tcn=2, lcn=1)."""
    n = 0
    sid = f"{patient}_{specimen}"
    for r in sh.itertuples(index=False):
        if pd.isna(getattr(r, f"{specimen}_vaf")):
            continue
        C, minor = local_cn(sid, r.chrom, r.pos)
        if pd.notna(C) and pd.notna(minor) and int(round(C)) == 2 and int(round(minor)) == 1:
            n += 1
    return n


def vaf_purity(patient, specimen, wgd):
    """QC purity from shared VAFs. Does not replace FACETS except for DH00024.

    With at least 10 shared mutations on tcn=2, lcn=1, clonal VAF is the median
    of shared VAFs at or above the sample median. Below that count, including
    DH00024, it is the median of every shared VAF (not the upper half).
    Diploid samples use p = 2v. Whole-genome-doubled samples use v = p/(2+p),
    so p = 2v/(1-v).
    """
    sh = V[(V["patient"] == patient) & (V["spatial_class"] == "shared")]
    if "call_set" in sh.columns:
        sh = sh[sh["call_set"].eq("PASS")]
    vafs = pd.to_numeric(sh[f"{specimen}_vaf"], errors="coerce").dropna()
    if vafs.empty:
        return np.nan
    if n_diploid_shared(patient, specimen, sh) < 10:
        v = float(vafs.median())
    else:
        upper = vafs[vafs >= vafs.median()]
        v = float(upper.median()) if len(upper) else float(vafs.median())
    if wgd:
        if v >= 1:
            return 1.0
        return float(min(1.0, (2 * v) / (1 - v)))
    return float(min(1.0, 2 * v))


# FACETS purity is the default. VAF overrides only the flat-genome samples DH00024_T and DH00024_J.
FLAT_VAF = {"DH00024_T", "DH00024_J"}
PURITY_SOURCE = {}
recon = []
for sid in [f"{p}_{s}" for p in sorted(V["patient"].unique()) for s in ("T", "J")]:
    patient, specimen = sid.split("_")
    facets = PP["purity"].get(sid, np.nan) if sid in PP.index else np.nan
    facets = float(facets) if pd.notna(facets) else np.nan
    raw_wgd = PP["WGD_heuristic"].get(sid, False) if sid in PP.index and "WGD_heuristic" in PP.columns else False
    wgd = str(raw_wgd).strip().lower() in {"true", "1", "yes"}
    vp = vaf_purity(patient, specimen, wgd)
    if sid in FLAT_VAF:
        used, source = vp, "vaf_override"
    else:
        used, source = facets, "facets"
    PURITY_SOURCE[sid] = source
    if sid in PP.index:
        PP.loc[sid, "purity"] = used
    recon.append({"sample": sid, "facets_purity": facets, "vaf_purity": vp,
                  "used_purity": used, "source": source})

recon_df = pd.DataFrame(recon)
recon_out = PROJECT_DIR / "07_cnv_purity" / "purity_reconciled.tsv"
recon_df.to_csv(recon_out, sep="\t", index=False)
write_provenance(recon_out, [PROJECT_DIR / "07_cnv_purity" / "purity_ploidy.tsv",
                             PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"],
                 {"rule": "FACETS purity except flat-genome DH00024_T/J. vaf_purity uses the upper half of shared VAFs when >=10 shared sites are tcn=2/lcn=1; otherwise the median of all shared VAFs"})
print(recon_df.round(3).to_string(index=False))

rows = []
for r in V.itertuples(index=False):
    rec = {"var_id": r.var_id}
    for s in ("T", "J"):
        sid = f"{r.patient}_{s}"
        p = PP["purity"].get(sid, np.nan)
        C, minor = local_cn(sid, r.chrom, r.pos)
        if pd.isna(C) and PURITY_SOURCE.get(sid) == "vaf_override":
            C, minor = 2.0, 1.0   # no FACETS fit: assume diploid
        major = C - minor if not pd.isna(minor) else C
        k, n = getattr(r, f"{s}_alt"), getattr(r, f"{s}_depth")
        v = k / n if n else np.nan
        ccf, m = ccf_from_vaf(v, p, C, major)
        lo_v, hi_v = wilson(k, n)
        lo = (lo_v * (p * C + 2 * (1 - p)) / p / m) if not pd.isna(ccf) else np.nan
        hi = (hi_v * (p * C + 2 * (1 - p)) / p / m) if not pd.isna(ccf) else np.nan
        rec.update({f"{s}_purity": p, f"{s}_purity_source": PURITY_SOURCE.get(sid, "missing"), f"{s}_tcn": C, f"{s}_lcn": minor, f"{s}_multiplicity": m,
                    f"{s}_ccf_raw": ccf, f"{s}_ccf": min(ccf, 1.0) if not pd.isna(ccf) else np.nan,
                    f"{s}_ccf_lo": min(lo, 1.0) if not pd.isna(lo) else np.nan,
                    f"{s}_ccf_hi": min(hi, 1.0) if not pd.isna(hi) else np.nan})
        if getattr(r, f"{s}_present"):
            rec[f"{s}_clonality"] = ("clonal" if rec[f"{s}_ccf_hi"] >= 0.9 else "subclonal") if not pd.isna(ccf) else "no_CN"
        else:
            rec[f"{s}_clonality"] = "not_present"
    rows.append(rec)

out = PROJECT_DIR / "09_ccf_clonality" / "ccf.tsv"
ccf = V.merge(pd.DataFrame(rows), on="var_id")
ccf.to_csv(out, sep="\t", index=False)
write_provenance(out, [PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", PROJECT_DIR / "07_cnv_purity" / "purity_ploidy.tsv"])
print(ccf.groupby(["patient", "T_clonality"]).size().unstack(fill_value=0))
