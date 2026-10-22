#!/usr/bin/env python3
"""Supplementary figures.

S2  Coding TMB per sample (T, J). Denominator is the unpadded 37.453 Mb target.
S3  Variant allele fraction distributions per sample (T, J), shared vs enriched.
The purity-versus-cellularity figure is not drawn: sequenced-section cellularity is not in the pathology table.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import (CLASS_COLOR, GRID, INK2, PROJECT_DIR, SPECIMEN_COLOR, WIDTH, placeholder,  # noqa: E402
                       save, setup)

setup()
P = PROJECT_DIR

# ---- S2 ---------------------------------------------------------------------
tf = P / "12_tables" / "main" / "tmb.tsv"
fig, ax = plt.subplots(figsize=(WIDTH["onehalf"], 60 / 25.4))
if tf.exists():
    t = pd.read_csv(tf, sep="\t")
    pts = sorted(t["patient"].unique())
    for k, s in enumerate(("J", "T")):
        sub = t[t["specimen"] == s].set_index("patient").reindex(pts)
        x = np.arange(len(pts)) + (k - 0.5) * 0.36
        ax.bar(x, sub["TMB_coding_per_mb"], width=0.34, color=SPECIMEN_COLOR[s], label=s)
    ax.set_xticks(range(len(pts)))
    ax.set_xticklabels(pts)
    ax.set_ylabel("Coding mutations / Mb")
    ax.yaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, frameon=False, borderaxespad=0)
    fig.text(0.01, 0.01, "coding mutations per unpadded target Mb (37.453)", fontsize=5.5, color=INK2)
else:
    placeholder(ax, "Awaiting TMB table")
save(fig, "FigS2_TMB", kind="supplementary", inputs=[tf])

# ---- S3 ---------------------------------------------------------------------
vf = P / "06_variants" / "somatic_filtered.tsv"
V = pd.read_csv(vf, sep="\t")
pts = sorted(V["patient"].unique())
fig, axes = plt.subplots(2, len(pts), figsize=(WIDTH["double"], 70 / 25.4), sharex=True, sharey="row")
axes = np.atleast_2d(axes)
lin_path = P / "09_ccf_clonality" / "DH00023_lineages.tsv"
lineage = {}
if lin_path.exists():
    lin = pd.read_csv(lin_path, sep="\t")
    lineage = dict(zip(lin.loc[lin["record"] == "variant", "var_id"], lin.loc[lin["record"] == "variant", "lineage"]))
bins = np.linspace(0, 0.7, 15)
legend_handles = {}
for j, p in enumerate(pts):
    for i, s in enumerate(("T", "J")):
        ax = axes[i, j]
        ax.set_xlim(0, 0.7)
        if p == "DH00024":
            ax.text(0.5, 0.5, "not evaluable", transform=ax.transAxes, ha="center", va="center", fontsize=6, color=INK2)
        else:
            sub = V[(V["patient"] == p) & V[f"{s}_present"]].copy()
            if p == "DH00023":
                sub["bucket"] = sub["var_id"].map(lineage).fillna(sub["spatial_class"])
                groups = [("trunk", CLASS_COLOR["shared"], "Trunk"),
                          ("J-lineage", CLASS_COLOR["J_private"], "J-lineage"),
                          ("T-lineage", CLASS_COLOR["T_private"], "T-lineage")]
            else:
                sub["bucket"] = sub["spatial_class"]
                groups = [("shared", CLASS_COLOR["shared"], "Shared"),
                          ("T_private", CLASS_COLOR["T_private"], "T-enriched"),
                          ("J_private", CLASS_COLOR["J_private"], "J-enriched")]
            arrays, colors, labels = [], [], []
            for key, color, label in groups:
                vals = sub.loc[sub["bucket"] == key, f"{s}_vaf"].dropna()
                if key == "T-lineage" and s == "J":
                    continue
                if key == "J-lineage" and s == "T" and p == "DH00023":
                    pass
                if len(vals) == 0 and key.endswith("private") and not key.startswith(s):
                    continue
                arrays.append(vals.to_numpy())
                colors.append(color)
                labels.append(label)
                legend_handles[label] = color
            if arrays:
                ax.hist(arrays, bins=bins, stacked=True, color=colors, label=labels)
        if i == 0:
            ax.set_title("DH00024" if p == "DH00024" else p, fontsize=6.5)
        if j == 0:
            ax.set_ylabel(f"{s}: variants")
        if i == 1:
            ax.set_xlabel("VAF")
axes[0, -1].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in legend_handles.values()],
                   labels=list(legend_handles), fontsize=5, loc="upper left", bbox_to_anchor=(1.02, 1))
save(fig, "FigS3_VAF_distributions", kind="supplementary", inputs=[vf])
