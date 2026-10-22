#!/usr/bin/env python3
"""Credible rescued sites, coloured by patient, with no gene labels."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import GRID, PROJECT_DIR, SPECIMEN_COLOR, WIDTH, panel_label, save, setup  # noqa: E402

setup()
sites_path = PROJECT_DIR / "08_lookback_N" / "credible_rescue_sites.tsv"
R = pd.read_csv(sites_path, sep="\t")
N = R.copy()
N["vaf"] = N["N_alt"] / N["N_depth"].replace(0, np.nan)
N["T_vaf"] = N["T_alt"] / N["T_depth"].replace(0, np.nan)
N["J_vaf"] = N["J_alt"] / N["J_depth"].replace(0, np.nan)
N["tj_vaf"] = N[["T_vaf", "J_vaf"]].max(axis=1)
colors = {
    "DH00018": "#3D6B66", "DH00019": "#6B4C3B", "DH00022": "#6B6B3A",
    "DH00023": "#6B5478", "DH00024": "#4D4D4D",
}

fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH["double"], 78 / 25.4), gridspec_kw={"wspace": 0.38, "width_ratios": [1.2, 1]})
for patient, g in N.groupby("patient"):
    a.scatter(g["tj_vaf"], g["vaf"].fillna(0) + 1e-3, s=18, color=colors.get(patient, "#4D4D4D"),
              edgecolor="white", linewidth=0.3, label=f"{patient} (n={len(g)})")
a.set_yscale("log")
a.set_ylim(8e-4, 1)
a.set_xlim(0, 1)
a.axhline(0.05, ls="--", color="#52514e", lw=0.5)
a.text(0.98, 0.062, "dashed line: VAF = 0.05", ha="right", va="bottom", fontsize=5.5, color="#52514e")
a.set_xlabel("Max VAF in tumour / junction")
a.set_ylabel("Normal VAF (zeros drawn at 0.001)")
a.set_title(f"Credible rescued sites (n={len(N)})")
a.yaxis.grid(True, color=GRID, lw=0.4)
a.legend(fontsize=5.5, frameon=False)
panel_label(a, "a", x=-0.18)

patients = list(dict.fromkeys(R["patient"]))
x = np.arange(len(patients))
width = 0.24
for k, (col, lab) in enumerate((("in_N", "N"), ("in_J", "J"), ("in_T", "T"))):
    counts = [int(R.loc[R["patient"] == p, col].sum()) for p in patients]
    b.bar(x + (k - 1) * width, counts, width=width, color=SPECIMEN_COLOR[lab], label=lab)
b.set_xticks(x)
b.set_xticklabels([p.replace("DH000", "") for p in patients])
b.set_ylabel("Sites with alternate reads")
b.set_xlabel("Patient")
b.legend(ncol=3, fontsize=6, frameon=False, loc="upper right")
b.yaxis.grid(True, color=GRID, lw=0.4)
b.set_axisbelow(True)
panel_label(b, "b", x=-0.2)
save(fig, "FigS_rescued", kind="supplementary", inputs=[sites_path],
     note="Dashed line is VAF 0.05. Zeros are drawn at 0.001. Read presence is not an independent validation of a field mutation. Not part of the 380.")
