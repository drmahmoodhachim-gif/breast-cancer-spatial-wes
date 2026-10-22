#!/usr/bin/env python3
"""Figure 1 — Study design and sequencing quality.

a  Sampling schematic (adjacent normal / junction / tumour from each resection).
b  Read pairs (millions) per sample.     c  Duplicate rate (%).
d  Mean on-target depth (x).             e  % targets >= 50x.
Input: 04_qc/qc_summary.tsv
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import (INK2, PROJECT_DIR, SPECIMEN_COLOR, SPECIMEN_NAME, WIDTH, panel_label,  # noqa: E402
                       placeholder, save, setup)

setup()
qc_f = PROJECT_DIR / "04_qc" / "qc_summary.tsv"
qc = pd.read_csv(qc_f, sep="\t") if qc_f.exists() else pd.DataFrame()

fig = plt.figure(figsize=(WIDTH["double"], 110 / 25.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.1, 1, 1], hspace=0.55, wspace=0.35)

# ---- a: schematic ---------------------------------------------------------------
ax = fig.add_subplot(gs[:, 0])
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")
for i, s in enumerate(["N", "J", "T"]):
    x = 0.4 + i * 3.1
    ax.add_patch(Rectangle((x, 4.4), 2.8, 3.0, facecolor=SPECIMEN_COLOR[s], alpha=0.18,
                           edgecolor=SPECIMEN_COLOR[s], lw=1))
    ax.text(x + 1.4, 6.3, s, ha="center", va="center", fontsize=11, fontweight="bold", color=SPECIMEN_COLOR[s])
    ax.text(x + 1.4, 5.1, SPECIMEN_NAME[s].replace(" ", "\n"), ha="center", va="center", fontsize=6)
ax.add_patch(FancyArrowPatch((0.6, 3.6), (9.4, 3.6), arrowstyle="-|>", mutation_scale=8, color=INK2, lw=0.7))
ax.text(5, 3.0, "towards tumour", ha="center", fontsize=6, color=INK2)
ax.text(5, 8.6, "Breast cancer resections (n = 5)", ha="center", fontsize=7, fontweight="bold")
ax.text(5, 1.8, "15 whole exomes (Illumina Exome 2.0 Plus)\nPE sequencing · GRCh38 · nf-core/sarek",
        ha="center", fontsize=6, color=INK2)
panel_label(ax, "a", x=0.0, y=0.95)

metrics = [("read_pairs_raw", "Read pairs (M)", 1e-6, "b"), ("pct_duplication", "Duplicates (%)", 1, "c"),
           ("mean_target_depth", "Mean on-target depth (×)", 1, "d"), ("pct_target_ge50x", "Targets ≥50× (%)", 1, "e")]
pos = [(0, 1), (0, 2), (1, 1), (1, 2)]
for (col, label, scale, letter), (r, c) in zip(metrics, pos):
    ax = fig.add_subplot(gs[r, c])
    panel_label(ax, letter)
    if qc.empty or col not in qc:
        placeholder(ax, f"{label}\n(awaiting sarek QC)")
        continue
    patients = sorted(qc["patient"].unique())
    w = 0.26
    for k, s in enumerate(["N", "J", "T"]):
        sub = qc[qc["specimen"] == s].set_index("patient").reindex(patients)
        x = np.arange(len(patients)) + (k - 1) * (w + 0.02)
        ax.bar(x, sub[col] * scale, width=w, color=SPECIMEN_COLOR[s], label=s, linewidth=0)
    ax.set_xticks(np.arange(len(patients)))
    ax.set_xticklabels([p.replace("DH000", "DH") for p in patients], rotation=0)
    ax.set_ylabel(label)
    ax.yaxis.grid(True, color="#e4e3df", lw=0.4)
    ax.set_axisbelow(True)
    if letter == "b":
        ax.legend(title=None, ncol=3, loc="upper left", bbox_to_anchor=(0, 1.3), handlelength=1, columnspacing=0.8,
                  labels=[f"{s}  {SPECIMEN_NAME[s]}" for s in ["N", "J", "T"]])

save(fig, "Fig1_design_QC", inputs=[qc_f] if qc_f.exists() else [])
