#!/usr/bin/env python3
"""Ascertainment check for the rescued normal-artifact sites.

The stored median ratios are the D054 simulation. Do not rerun this script to
refresh them: the look-back dispersion prior has since been refit. A rerun would
replace those medians. The observed column is per patient.

Simulates normal alternate reads at primary shared sites from a beta-binomial
whose mean is the tumour-in-normal slope times the tumour allele fraction,
using the look-back dispersion. Mutect2's normal_artifact flag is not a fixed
alt-count in this project, so the filter is applied as N alt >= k for k = 3..6.
The median normal/tumour allele-fraction ratio among sites that would be
filtered is compared with the observed ratio of the 27 credible rescued sites.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import betabinom

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

N_REP = 400
KS = (3, 4, 5, 6)
rng = np.random.default_rng(54)

prior = pd.read_csv(PROJECT_DIR / "08_lookback_N" / "background_prior.tsv", sep="\t").iloc[0]
conc = float(prior["alpha0"] + prior["beta0"])
tin = pd.read_csv(PROJECT_DIR / "08_lookback_N" / "tumour_in_normal_estimates.tsv", sep="\t")
slope = dict(zip(tin["patient"], tin["TiN_shared"]))
var = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
shared = var[(var["spatial_class"] == "shared") & (var["T_vaf"] > 0) & (var["N_depth"] > 0)].copy()
qc = pd.read_csv(PROJECT_DIR / "08_lookback_N" / "credible_rescue_qc.tsv", sep="\t")
per_patient = qc.groupby("patient")["n_over_tumour"].median()
pooled = float(qc["n_over_tumour"].median())

rows = []
for patient, g in list(shared.groupby("patient")) + [("all", shared)]:
    observed = pooled if patient == "all" else (float(per_patient[patient]) if patient in per_patient.index else np.nan)
    if patient == "all":
        sl = g["patient"].map(slope).to_numpy(float)
    else:
        sl = np.full(len(g), float(slope[patient]))
    depth = g["N_depth"].to_numpy(int)
    tv = np.clip(g["T_vaf"].to_numpy(float), 1e-6, 1)
    mu = np.clip(sl * tv, 1e-8, 1 - 1e-8)
    a = mu * conc
    b = (1 - mu) * conc
    medians = {k: [] for k in KS}
    fractions = {k: [] for k in KS}
    for _ in range(N_REP):
        alt = betabinom.rvs(depth, a, b, random_state=rng)
        ratio = (alt / depth) / tv
        for k in KS:
            keep = alt >= k
            fractions[k].append(float(keep.mean()))
            if keep.any():
                medians[k].append(float(np.median(ratio[keep])))
    for k in KS:
        vals = np.asarray(medians[k], float)
        rows.append({
            "patient": patient,
            "n_shared_sites": int(len(g)),
            "k_alt": k,
            "fraction_filtered_mean": float(np.mean(fractions[k])),
            "median_NT_ratio": float(np.median(vals)) if len(vals) else np.nan,
            "median_NT_ratio_lo": float(np.percentile(vals, 2.5)) if len(vals) else np.nan,
            "median_NT_ratio_hi": float(np.percentile(vals, 97.5)) if len(vals) else np.nan,
            "observed_rescued_median_NT": observed,
            "observed_rescued_median_NT_pooled": pooled,
            "note": "Illustration. observed_rescued_median_NT is this patient's rescued median; the all row is the pooled median. The simulation does not apply N VAF <= 3x expected, absence in other normals, or gnomAD.",
        })

out = pd.DataFrame(rows)
path = PROJECT_DIR / "08_lookback_N" / "rescued_ascertainment.tsv"
out.to_csv(path, sep="\t", index=False)
write_provenance(path, [
    PROJECT_DIR / "06_variants" / "somatic_filtered.tsv",
    PROJECT_DIR / "08_lookback_N" / "background_prior.tsv",
    PROJECT_DIR / "08_lookback_N" / "credible_rescue_qc.tsv",
], {"replicates": N_REP, "concentration": conc, "observed_median": observed})
print(out.round(4).to_string(index=False))
print(f"[ascertain] observed rescued median N/T = {observed:.4f}")
