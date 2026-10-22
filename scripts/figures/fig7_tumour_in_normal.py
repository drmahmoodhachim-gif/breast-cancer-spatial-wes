#!/usr/bin/env python3
"""Figure — relative allele-fraction slope in adjacent normal (tumour-in-normal, TiN).

If tumour cells (or tumour DNA) are present in N, every somatic variant appears in N at
a VAF proportional to its VAF in the tumour:  VAF_N ≈ TiN × VAF_T.
a  Per patient: VAF in N vs VAF in T for every variant present in T (shared and T-private),
   with the weighted least-squares slope through the origin (= TiN estimate).
b  TiN estimate per patient with 95 % bootstrap CI (resampling variants), shared variants only.
Background: the mean site error rate from unrelated patients' samples (dashed).
Input: 06_variants/somatic_filtered.tsv, 08_lookback_N/background_prior.tsv
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import (CLASS_COLOR, GRID, INK2, PROJECT_DIR, WIDTH, panel_label, save, setup)  # noqa: E402
from scipy.stats import betabinom

setup()
vf = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
bf = PROJECT_DIR / "08_lookback_N" / "background_prior.tsv"
legacy = PROJECT_DIR / "08_lookback_N" / "tumour_in_normal_estimates_reads_only.tsv"
current_est = PROJECT_DIR / "08_lookback_N" / "tumour_in_normal_estimates.tsv"
if current_est.exists() and not legacy.exists():
    legacy.write_bytes(current_est.read_bytes())
V = pd.read_csv(vf, sep="\t")
prior = pd.read_csv(bf, sep="\t") if bf.exists() else pd.DataFrame()
bg = float(prior["mean_error_rate"].iloc[0]) if len(prior) else np.nan
a0 = float(prior["alpha0"].iloc[0]) if len(prior) and "alpha0" in prior.columns else 1.0
b0 = float(prior["beta0"].iloc[0]) if len(prior) and "beta0" in prior.columns else 500.0
rng = np.random.default_rng(1234)
N_BOOT = 1000
N_SIM = 10000


def tin_on(df, vaf_col):
    """Read-weighted slope through the origin: sum(N_alt) / sum(N_depth * VAF)."""
    if df.empty:
        return np.nan
    den = (df["N_depth"] * df[vaf_col]).sum()
    return df["N_alt"].sum() / den if den > 0 else np.nan


def tin(df):
    return tin_on(df, "T_vaf")


def wls_intercept(df):
    """Weighted N_vaf = a + b * T_vaf. Weights are N depth."""
    if len(df) < 3:
        return np.nan, np.nan
    x = df["T_vaf"].to_numpy(float)
    y = df["N_vaf"].to_numpy(float)
    w = df["N_depth"].to_numpy(float)
    if not np.isfinite(w).any() or np.nansum(w) <= 0:
        return np.nan, np.nan
    X = np.column_stack([np.ones(len(x)), x])
    xtw = X.T * w
    beta = np.linalg.lstsq(xtw @ X, xtw @ y, rcond=None)[0]
    return float(beta[0]), float(beta[1])


def pct_ci(xs):
    xs = [x for x in xs if np.isfinite(x)]
    if not xs:
        return np.nan, np.nan
    lo, hi = np.percentile(xs, [2.5, 97.5])
    return float(lo), float(hi)


PREV_A, PREV_B = 0.049785335388002075, 8.588851328794249


def parametric_branch(private, vaf_col, point_slope, slope_draws, alpha, beta):
    """Lower-tail test that N alt reads are below the TiN expectation. D042, D067."""
    empty = {"expected": np.nan, "p": np.nan, "label": "too few reads to test"}
    if private.empty or not np.isfinite(point_slope):
        return empty
    depth = private["N_depth"].to_numpy(int)
    vaf = np.clip(private[vaf_col].to_numpy(float), 0, 1)
    obs = int(private["N_alt"].sum())
    expected = float(np.sum(depth * point_slope * vaf))
    if not np.isfinite(expected) or expected < 5 or len(slope_draws) == 0:
        return {"expected": expected, "p": np.nan, "label": "too few reads to test"}
    conc = alpha + beta
    draws = rng.choice(np.asarray(slope_draws, float), size=N_SIM, replace=True)
    sim = np.empty(N_SIM)
    for i, sl in enumerate(draws):
        mu = np.clip(sl * vaf, 1e-8, 1 - 1e-8)
        sim[i] = betabinom.rvs(depth, mu * conc, (1 - mu) * conc, random_state=rng).sum()
    p = float(np.mean(sim <= obs))
    label = "p < 1e-4" if p < 1e-4 else f"p = {p:.3g}"
    return {"expected": expected, "p": p, "label": label}


rows = []
for p, g in V.groupby("patient"):
    s = g[(g["spatial_class"] == "shared") & (g["T_vaf"] > 0)].reset_index(drop=True)
    est = tin(s)
    intercept, slope_free = wls_intercept(s)
    boots_tin, boots_a, boots_b = [], [], []
    if len(s) >= 5:
        for _ in range(N_BOOT):
            sb = s.iloc[rng.integers(0, len(s), size=len(s))]
            boots_tin.append(tin(sb))
            a, b = wls_intercept(sb)
            boots_a.append(a)
            boots_b.append(b)
    lo, hi = pct_ci(boots_tin)
    alo, ahi = pct_ci(boots_a)
    blo, bhi = pct_ci(boots_b)
    priv = g[(g["spatial_class"] == "T_private") & (g["T_vaf"] > 0)]
    den = float((s["N_depth"] * s["T_vaf"]).sum())
    error_null = float((s["N_depth"] * bg).sum() / den) if den > 0 and np.isfinite(bg) else np.nan
    rec = {"patient": p, "n_shared": len(s), "TiN_shared": est, "ci_lo": lo, "ci_hi": hi,
           "error_null_slope": error_null,
           "intercept": intercept, "intercept_ci_lo": alo, "intercept_ci_hi": ahi,
           "TiN_intercept_model": slope_free, "TiN_intercept_ci_lo": blo, "TiN_intercept_ci_hi": bhi,
           "n_T_private": len(priv), "TiN_T_private": tin(priv) if len(priv) else np.nan}
    sj = s[s["J_vaf"] > 0].reset_index(drop=True)
    tin_j = tin_on(sj, "J_vaf")
    boots_j = [tin_on(sj.iloc[rng.integers(0, len(sj), size=len(sj))], "J_vaf")
               for _ in range(N_BOOT)] if len(sj) >= 5 else []
    slope_bank = {"T_private": [x for x in boots_tin if np.isfinite(x)],
                  "J_private": [x for x in boots_j if np.isfinite(x)]}
    point = {"T_private": est, "J_private": tin_j}
    for cls, ref in (("T_private", "T_vaf"), ("J_private", "J_vaf")):
        q = g[g["spatial_class"] == cls].reset_index(drop=True)
        tested = parametric_branch(q, ref, point[cls], slope_bank[cls], a0, b0)
        previous = parametric_branch(q, ref, point[cls], slope_bank[cls], PREV_A, PREV_B)
        rec[f"{cls}_n"] = len(q)
        rec[f"{cls}_N_alt_expected"] = tested["expected"]
        rec[f"{cls}_N_alt_observed"] = int(q["N_alt"].sum()) if len(q) else 0
        rec[f"{cls}_p_depleted"] = tested["p"]
        rec[f"{cls}_p_label"] = tested["label"]
        rec[f"{cls}_p_previous"] = previous["p"]
        rec[f"{cls}_p_previous_label"] = previous["label"]
    rows.append(rec)
T = pd.DataFrame(rows)
out_tab = PROJECT_DIR / "08_lookback_N" / "tumour_in_normal_estimates.tsv"
T.to_csv(out_tab, sep="\t", index=False)
print(T.round(4).to_string(index=False))

patients = list(T["patient"])
fig = plt.figure(figsize=(WIDTH["double"], 92 / 25.4))
gs = fig.add_gridspec(2, len(patients), height_ratios=[1, 0.85], hspace=0.75, wspace=0.12)
for k, p in enumerate(patients):
    ax = fig.add_subplot(gs[0, k])
    g = V[(V["patient"] == p) & (V["T_vaf"] > 0)]
    for c in ("T_private", "shared"):
        s = g[g["spatial_class"] == c]
        ax.scatter(s["T_vaf"], s["N_vaf"] * 100, s=9, color=CLASS_COLOR[c], edgecolor="white", lw=0.4,
                   zorder=3 if c == "shared" else 2)
    est = float(T.loc[T["patient"] == p, "TiN_shared"].iloc[0])
    xs = np.linspace(0, 0.7, 10)
    if np.isfinite(est):
        ax.plot(xs, xs * est * 100, color=INK2, lw=0.8, zorder=4)
        ax.text(0.03, 0.93, f"relative slope {est:.3f}", transform=ax.transAxes, fontsize=6, va="top", color="#1a1a1a")
    ax.set_xlim(0, 0.7)
    ax.set_ylim(-0.2, 5)
    ax.set_title(p, fontsize=6.5)
    ax.set_xlabel("Tumour VAF")
    ax.yaxis.grid(True, color=GRID, lw=0.4)
    if k == 0:
        ax.set_ylabel("Adjacent-normal VAF (%)")
        panel_label(ax, "a", x=-0.45)
    else:
        ax.set_yticklabels([])
ax.scatter([], [], s=12, color=CLASS_COLOR["shared"], label="Shared T+J")
ax.scatter([], [], s=12, color=CLASS_COLOR["T_private"], label="T-enriched")
ax.plot([], [], color=INK2, lw=0.8, label="Fit through origin")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=5.5)

b = fig.add_subplot(gs[1, :])
y = np.arange(len(T))
b.errorbar(T["TiN_shared"], y, xerr=[T["TiN_shared"] - T["ci_lo"], T["ci_hi"] - T["TiN_shared"]],
           fmt="o", color=CLASS_COLOR["shared"], ms=4, capsize=2, lw=0.8)
for i, r in T.iterrows():
    b.text((r["ci_hi"] if np.isfinite(r["ci_hi"]) else r["TiN_shared"]) + 0.0015, i,
           f"{r['TiN_shared']:.3f}  (n = {r['n_shared']})", va="center", fontsize=6, color=INK2)
b.scatter(T["error_null_slope"], y, marker="x", s=16, color=INK2, linewidths=0.6, zorder=4)
b.text(0.0, -0.72, "cross: error-only null on this slope scale", fontsize=5.5, color=INK2, va="center", clip_on=False)
b.set_yticks(y)
b.set_yticklabels(T["patient"])
b.invert_yaxis()
b.set_ylim(len(T) - 0.45, -0.85)
b.set_xlim(0, max(0.07, np.nanmax(T["ci_hi"]) + 0.025))
b.set_xlabel("Relative allele fraction in adjacent normal (normal / tumour)")
b.xaxis.grid(True, color=GRID, lw=0.4)
b.set_axisbelow(True)
panel_label(b, "b", x=-0.08)
dh = T[T["patient"] == "DH00023"].iloc[0]
fig.text(0.01, 0.008,
         "DH00023 joint presence rule, T-enriched: "
         f"{int(dh['T_private_N_alt_observed'])} reads vs {dh['T_private_N_alt_expected']:.1f} expected, "
         f"refit {dh['T_private_p_label']}; previous dispersion {dh['T_private_p_previous_label']}. "
         f"J-enriched: {int(dh['J_private_N_alt_observed'])} vs {dh['J_private_N_alt_expected']:.1f}, "
         f"refit {dh['J_private_p_label']}. "
         "The reads-only 46-site refit was p = 0.011 and is a different set. "
         "These are not cell fractions and do not show invasion.",
         fontsize=5.0, color=INK2)
save(fig, "Fig7_tumour_in_normal", inputs=[vf, bf],
     note="D067. Slopes use the joint-rule classes. Branch tests are 10,000 beta-binomial draws with slope bootstrap, refit prior and previous prior. The 46-site p = 0.011 is retained only as the reads-only comparison.")
