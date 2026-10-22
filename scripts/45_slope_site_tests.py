#!/usr/bin/env python3
"""D067. Slope and positive-site tests for the joint-rule DH00023 sets.

The beta-binomial slope-bootstrap p-values are written by Figure 7.
This script adds the Poisson read test and the fixed-slope beta-binomial
at the recorded slopes, plus the positive-site expectations already stored
in dh00023_detection_joint_rule.tsv. It rebuilds branch_test_refit.tsv
without dropping the reads-only rows.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import betabinom

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

ROOT = PROJECT_DIR
REFIT_A, REFIT_B = 0.015705149447423628, 22.300145625836034
PREV_A, PREV_B = 0.049785335388002075, 8.588851328794249
N_SIM = 10000
SLOPES = [
    ("shared_through_origin", 0.026174665261835),
    ("error_subtracted", 0.0178),
    ("trunk_only", 0.0235),
    ("trunk_only_bootstrap_low", 0.0101),
    ("trunk_only_bootstrap_high", 0.0384),
]


def poisson_lower(expected: float, observed: int) -> float:
    if not np.isfinite(expected) or expected < 0:
        return np.nan
    # P(X <= observed). For zero reads this is exp(-expected).
    from scipy.stats import poisson
    return float(poisson.cdf(observed, expected))


def fixed_beta_p(depth, vaf, obs, slope, alpha, beta, seed):
    depth = np.asarray(depth, int)
    mu = np.clip(slope * np.asarray(vaf, float), 1e-8, 1 - 1e-8)
    conc = alpha + beta
    rng = np.random.default_rng(seed)
    sim = np.empty(N_SIM)
    for i in range(N_SIM):
        sim[i] = betabinom.rvs(depth, mu * conc, (1 - mu) * conc, random_state=rng).sum()
    return float(np.mean(sim <= obs))


def main():
    V = pd.read_csv(ROOT / "06_variants" / "somatic_filtered.tsv", sep="\t")
    det = pd.read_csv(ROOT / "08_lookback_N" / "dh00023_detection_joint_rule.tsv", sep="\t")
    old = pd.read_csv(ROOT / "08_lookback_N" / "branch_test_refit.tsv", sep="\t")
    if "rule" not in old.columns:
        old.insert(0, "rule", "reads_only")
        old.insert(1, "test", "beta_binomial_slope_bootstrap")
    old = old[old["rule"] == "reads_only"]
    d = V[V["patient"] == "DH00023"]
    sets = {
        "T_private": d[d["spatial_class"] == "T_private"],
        "J_private": d[d["spatial_class"] == "J_private"],
    }
    base_expected = {
        "T_private": float(det.loc[det["group"] == "joint_T_enriched", "N_alt_expected_at_slope"].iloc[0]),
        "J_private": float(det.loc[det["group"] == "joint_J_enriched", "N_alt_expected_at_slope"].iloc[0]),
    }
    base_slope = 0.026174665261835
    rows = []
    for branch, g in sets.items():
        vaf_col = "T_vaf" if branch == "T_private" else "J_vaf"
        depth = g["N_depth"].to_numpy(int)
        vaf = g[vaf_col].fillna(0).to_numpy(float)
        obs = int(g["N_alt"].fillna(0).sum())
        for i, (name, slope) in enumerate(SLOPES):
            expected = base_expected[branch] * slope / base_slope
            rows.append({
                "rule": "joint",
                "test": "poisson_read_count",
                "patient": "DH00023",
                "branch": branch,
                "n": len(g),
                "N_alt_observed": obs,
                "N_alt_expected": expected,
                "p_previous": np.nan,
                "label_previous": "",
                "p_refit": poisson_lower(expected, obs),
                "label_refit": name,
                "alpha_previous": PREV_A,
                "beta_previous": PREV_B,
                "alpha_refit": REFIT_A,
                "beta_refit": REFIT_B,
                "primary": "no",
                "slope": slope,
                "note": "Poisson lower tail P(X<=observed). Ignores read-to-read variability. Not the reported test.",
            })
            if branch == "T_private":
                p_refit = fixed_beta_p(depth, vaf, obs, slope, REFIT_A, REFIT_B, 3000 + i)
                p_prev = fixed_beta_p(depth, vaf, obs, slope, PREV_A, PREV_B, 4000 + i)
                rows.append({
                    "rule": "joint",
                    "test": "beta_binomial_slope_fixed",
                    "patient": "DH00023",
                    "branch": branch,
                    "n": len(g),
                    "N_alt_observed": obs,
                    "N_alt_expected": expected,
                    "p_previous": p_prev,
                    "label_previous": f"p = {p_prev:.3g}",
                    "p_refit": p_refit,
                    "label_refit": f"p = {p_refit:.3g}",
                    "alpha_previous": PREV_A,
                    "beta_previous": PREV_B,
                    "alpha_refit": REFIT_A,
                    "beta_refit": REFIT_B,
                    "primary": "sensitivity",
                    "slope": slope,
                    "note": "Beta-binomial with the slope held fixed. The reported test is the slope-bootstrap row from Figure 7.",
                })
    for _, r in det.iterrows():
        if not str(r["group"]).startswith("joint_"):
            continue
        rows.append({
            "rule": "joint",
            "test": "positive_sites",
            "patient": "DH00023",
            "branch": r["group"],
            "n": r["n_sites"],
            "N_alt_observed": r["n_positive_sites"],
            "N_alt_expected": r["expected_positive_binomial"],
            "p_previous": r["expected_positive_betabinom"],
            "label_previous": "expected_positive_betabinom",
            "p_refit": r["poisson_binomial_p_lower"],
            "label_refit": "poisson_binomial_p_lower",
            "alpha_previous": PREV_A,
            "beta_previous": PREV_B,
            "alpha_refit": REFIT_A,
            "beta_refit": REFIT_B,
            "primary": "site_model",
            "slope": r["slope"],
            "note": (
                f"Observed positive sites in N_alt_observed. "
                f"N_alt_expected is the binomial expectation. "
                f"p_previous column here is the beta-binomial expectation of positive sites, not a p-value. "
                f"p_refit is the exact Poisson-binomial lower tail. "
                f"Beta-binomial positive-site p = {r['beta_binomial_positive_site_p']}. "
                f"Plus-error Poisson-binomial p = {r['poisson_binomial_p_lower_plus_error']}."
            ),
        })
    est_path = ROOT / "08_lookback_N" / "tumour_in_normal_estimates.tsv"
    if est_path.exists():
        est = pd.read_csv(est_path, sep="\t")
        # Only treat this file as the joint-rule draw when Figure 7 has been rerun
        # on the smaller enriched sets. The reads-only file has DH00023 T n = 46.
        if "T_private_n" in est.columns and int(est.loc[est["patient"] == "DH00023", "T_private_n"].iloc[0]) < 46:
            for _, r in est.iterrows():
                for branch, prefix in (("T_private", "T_private"), ("J_private", "J_private")):
                    rows.append({
                        "rule": "joint",
                        "test": "beta_binomial_slope_bootstrap",
                        "patient": r["patient"],
                        "branch": branch,
                        "n": r.get(f"{prefix}_n", np.nan),
                        "N_alt_observed": r.get(f"{prefix}_N_alt_observed", np.nan),
                        "N_alt_expected": r.get(f"{prefix}_N_alt_expected", np.nan),
                        "p_previous": r.get(f"{prefix}_p_previous", np.nan),
                        "label_previous": r.get(f"{prefix}_p_previous_label", ""),
                        "p_refit": r.get(f"{prefix}_p_depleted", np.nan),
                        "label_refit": r.get(f"{prefix}_p_label", ""),
                        "alpha_previous": PREV_A,
                        "beta_previous": PREV_B,
                        "alpha_refit": REFIT_A,
                        "beta_refit": REFIT_B,
                        "primary": "yes",
                        "slope": r.get("TiN_shared", np.nan),
                        "note": "Figure 7 method: 10,000 beta-binomial draws, slope resampled from the shared-site bootstrap. Reported test.",
                    })
    out = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    path = ROOT / "08_lookback_N" / "branch_test_refit.tsv"
    out.to_csv(path, sep="\t", index=False)
    write_provenance(path, inputs=[ROOT / "06_variants" / "somatic_filtered.tsv", est_path], params={"decision": "D067"})
    print(out[out["rule"] == "joint"][["test", "branch", "n", "N_alt_observed", "N_alt_expected", "p_refit", "label_refit"]].to_string(index=False))


if __name__ == "__main__":
    main()
