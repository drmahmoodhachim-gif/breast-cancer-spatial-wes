#!/usr/bin/env python3
"""Site-level error model for the look-back of T/J variants into adjacent normal (N).

For each candidate locus discovered in T and/or J of patient P:
  1. Background = the same locus in samples of the OTHER four patients, excluding any
     sample whose patient also carries a somatic call there, and any background sample
     with VAF >= 0.2 (likely a germline allele in that person).
  2. A global beta-binomial prior (alpha0, beta0) is fitted by maximum likelihood to all
     background observations; each site's error distribution is that prior updated with the
     site's own background alt/ref counts (empirical Bayes).
  3. P-value for P's N (and J, T) = upper tail of the site's beta-binomial predictive.
  4. Benjamini–Hochberg FDR across loci within each specimen.
  5. Artefact flags: strand bias, alt reads clustered at read ends.

Classes for N (written to the paper exactly like this):
  N_confident          q < alpha, VAF >= lookback_confident_vaf, alt >= min_alt, no flags
  N_supporting_reads   q < alpha, but below the confident VAF/alt thresholds, no flags
  N_flagged_artefact   q < alpha but strand/position flags
  N_no_evidence        q >= alpha
  N_low_depth          depth < min_depth_normal

Outputs: 08_lookback_N/lookback_results.tsv, lookback_summary.tsv, background_prior.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import betabinom

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_filters, write_provenance  # noqa: E402

F = load_filters()
D = PROJECT_DIR / "08_lookback_N"
counts = pd.read_csv(D / "allele_counts.tsv", sep="\t")
var = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
counts = counts[counts["var_id"].isin(set(var["var_id"]))].copy()

var["locus"] = var["chrom"] + ":" + var["pos"].astype(str) + ":" + var["ref"] + ">" + var["alt"]
counts = counts.merge(var[["var_id", "locus"]], on="var_id", how="left")
called = set(zip(var["patient"], var["locus"]))
counts["vaf"] = counts["alt"] / counts["depth"].replace(0, np.nan)

own = counts[counts["sample_patient"] == counts["variant_patient"]].copy()
bg = counts[(counts["sample_patient"] != counts["variant_patient"]) & (counts["depth"] > 0)].copy()
bg = bg[[(p, l) not in called for p, l in zip(bg["sample_patient"], bg["locus"])]]
bg = bg[bg["vaf"].fillna(0) < 0.2]

# ---- global prior by ML -------------------------------------------------------
k, n = bg["alt"].to_numpy(int), bg["depth"].to_numpy(int)


def nll(theta):
    a, b = np.exp(theta)
    return -betabinom.logpmf(k, n, a, b).sum()


fit = minimize(nll, x0=np.log([1.0, 500.0]), method="Nelder-Mead")
a0, b0 = np.exp(fit.x)
pd.DataFrame({"alpha0": [a0], "beta0": [b0], "mean_error_rate": [a0 / (a0 + b0)],
              "n_background_obs": [len(bg)]}).to_csv(D / "background_prior.tsv", sep="\t", index=False)
print(f"[lookback] background prior alpha={a0:.3f} beta={b0:.1f} (mean error {a0/(a0+b0):.2e}, {len(bg)} obs)")

site = bg.groupby("var_id").agg(bg_alt=("alt", "sum"), bg_depth=("depth", "sum"), bg_n=("sample_id", "nunique"))
own = own.merge(site, on="var_id", how="left").fillna({"bg_alt": 0, "bg_depth": 0, "bg_n": 0})
own["site_alpha"] = a0 + own["bg_alt"]
own["site_beta"] = b0 + (own["bg_depth"] - own["bg_alt"])
own["site_error_rate"] = own["site_alpha"] / (own["site_alpha"] + own["site_beta"])
own["p_value"] = [float(betabinom.sf(x - 1, m, a, b)) if m > 0 and x > 0 else 1.0
                  for x, m, a, b in zip(own["alt"], own["depth"], own["site_alpha"], own["site_beta"])]


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


own["q_value"] = own.groupby("specimen")["p_value"].transform(bh)
own["flag_strand"] = (own["alt"] >= 4) & ((own["alt_fwd"] < F["lookback_min_alt_per_strand"]) |
                                           (own["alt_rev"] < F["lookback_min_alt_per_strand"]))
own["flag_read_end"] = own["alt_median_dist_end"].fillna(99) < F["lookback_min_median_dist_from_read_end"]


def classify(r):
    if r["depth"] < F["min_depth_normal"]:
        return "low_depth"
    if r["q_value"] >= F["lookback_alpha_fdr"]:
        return "no_evidence"
    if r["flag_strand"] or r["flag_read_end"]:
        return "flagged_artefact"
    if r["vaf"] >= F["lookback_confident_vaf"] and r["alt"] >= F["lookback_min_alt_confident"]:
        return "confident"
    return "supporting_reads"


own["evidence_class"] = own.apply(classify, axis=1)
res = own.merge(var[["var_id", "gene", "hgvsp", "consequence", "is_driver_gene", "spatial_class",
                     "T_vaf", "J_vaf"]], on="var_id", how="left")
res.to_csv(D / "lookback_results.tsv", sep="\t", index=False)

summary = (res[res["specimen"] == "N"].groupby(["variant_patient", "evidence_class"]).size()
           .unstack(fill_value=0).rename_axis("patient"))
summary.to_csv(D / "lookback_summary.tsv", sep="\t")
write_provenance(D / "lookback_results.tsv", [D / "allele_counts.tsv"], {**{k: F[k] for k in F if k.startswith("lookback")},
                                                                         "alpha0": a0, "beta0": b0})
print(summary)
