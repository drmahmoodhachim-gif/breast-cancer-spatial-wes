#!/usr/bin/env python3
"""Recovered tumour-in-normal slope against the expected increment from in-silico dilution.

recovered = observed slope − original count-based baseline.
expected increment = expected slope − that same baseline.
The identity line is recovered = expected increment.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import INK2, PROJECT_DIR, WIDTH, save, setup  # noqa: E402

setup()
tsv = PROJECT_DIR / "08_lookback_N" / "dilution_test.tsv"
df = pd.read_csv(tsv, sep="\t")
if "recovered" not in df.columns:
    raise SystemExit("dilution_test.tsv has no recovered column")
if (df["observed_slope"] > 0.2).any():
    raise SystemExit("a row looks like the rejected full-tumour merge")
colors = {"DH00018": "#0072B2", "DH00023": "#882255"}
fig, ax = plt.subplots(figsize=(WIDTH["onehalf"], WIDTH["onehalf"] * 0.9))
x = df["expected_increment"]
y_lo = df["ci_lo"] - df["baseline_slope"]
y_hi = df["ci_hi"] - df["baseline_slope"]
lo = min(x.min(), y_lo.min(), 0)
hi = max(x.max(), y_hi.max())
pad = 0.004
lims = [lo - pad, hi + pad]
ax.plot(lims, lims, color="#3d3d3a", lw=0.6, zorder=1)
for patient, g in df.groupby("patient"):
    yerr = [g["recovered"] - (g["ci_lo"] - g["baseline_slope"]),
            (g["ci_hi"] - g["baseline_slope"]) - g["recovered"]]
    ax.errorbar(
        g["expected_increment"], g["recovered"], yerr=yerr, fmt="o",
        color=colors[patient], ms=3.5, lw=0.6, capsize=2, capthick=0.5,
        label=patient, zorder=2,
    )
    fit = g.iloc[0]
    xs = pd.Series(lims)
    ax.plot(xs, fit["fit_intercept"] + fit["fit_slope"] * xs,
            color=colors[patient], lw=0.5, ls="--", zorder=1, label=f"{patient} fit")
ax.plot([], [], color="#3d3d3a", lw=0.6, label="Identity")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Expected increment")
ax.set_ylabel("Recovered")
leg = ax.legend(loc="upper left", fontsize=6)
leg.set_title("recovered = observed − original estimate;\n0% point shows the merge-only offset", prop={"size": 5.5})
ax.tick_params(colors=INK2)
save(fig, "FigS_dilution", kind="supplementary", inputs=[tsv],
     note="D047. Identity line is solid. x is expected slope minus the original baseline; y is observed slope minus that baseline. The 0% point is a header-only tumour BAM. Error bars are the observed-slope bootstrap shifted by the baseline.")
