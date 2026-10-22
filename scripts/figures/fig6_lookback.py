#!/usr/bin/env python3
"""Figure 6 — Targeted look-back of T/J variants into adjacent normal tissue.

a  For every candidate locus: VAF in N (y) vs max VAF in T/J (x), shape+colour = evidence class.
   Dashed line = 5 % VAF (confident-call floor at exome depth).
b  Per patient, number of loci in each N evidence class.
Supporting-read observations (< 5 % VAF) are hypotheses for orthogonal validation, not calls.
Input: 08_lookback_N/lookback_results.tsv restricted to the 380 PASS variants in somatic_filtered.tsv.
Rescued normal_artifact sites are not drawn here.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import GRID, PROJECT_DIR, SPECIMEN_COLOR, WIDTH, panel_label, save, setup  # noqa: E402

setup()
lf = PROJECT_DIR / "08_lookback_N" / "lookback_results.tsv"
primary = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
L = pd.read_csv(lf, sep="\t")
keep = set(pd.read_csv(primary, sep="\t")["var_id"])
L = L[L["var_id"].isin(keep)]
N = L[L["specimen"] == "N"].copy()
N["tj_vaf"] = N[["T_vaf", "J_vaf"]].max(axis=1)

style = {  # class: (label, colour, marker, size)
    "confident": ("Confident (q<0.05, VAF≥5 %)", SPECIMEN_COLOR["N"], "s", 16),
    "supporting_reads": ("Supporting reads (q<0.05, <5 %)", "#7fb6dc", "o", 12),
    "flagged_artefact": ("Flagged artefact", "#8a8984", "x", 12),
    "no_evidence": ("No evidence", "#c9c8c2", ".", 6),
    "low_depth": ("Low depth in N", "#e4e3df", ".", 6),
}
fig = plt.figure(figsize=(WIDTH["double"], 78 / 25.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.05, 0.72], wspace=0.55)
a = fig.add_subplot(gs[0])
b = fig.add_subplot(gs[1])
c = fig.add_subplot(gs[2])
for cls in ["no_evidence", "low_depth", "flagged_artefact", "supporting_reads", "confident"]:
    s = N[N["evidence_class"] == cls]
    lab, col, mk, sz = style[cls]
    a.scatter(s["tj_vaf"], s["vaf"].fillna(0) + 1e-3, marker=mk, s=sz, color=col,
              edgecolor="white" if mk in "so" else None, linewidth=0.3, label=f"{lab} (n={len(s)})")
for r in N[(N["evidence_class"].isin(["confident", "supporting_reads"])) & N["is_driver_gene"] & (N["vaf"].fillna(0) >= 0.02)].itertuples():
    a.annotate(f"{r.gene}", (r.tj_vaf, r.vaf + 1e-3), xytext=(3, 2), textcoords="offset points",
               fontsize=5.5, fontstyle="italic")
a.set_yscale("log")
a.set_ylim(8e-4, 1)
a.set_xlim(0, 1)
a.axhline(0.05, ls="--", color="#52514e", lw=0.5)
a.text(0.99, 0.058, "line at 0.05", ha="right", va="bottom", fontsize=5.5, color="#52514e")
a.set_xlabel("Max VAF in tumour / junction")
a.set_ylabel("VAF in adjacent normal (+0.001)")
a.yaxis.grid(True, color=GRID, lw=0.4)
panel_label(a, "a", x=-0.2)

order = ["confident", "supporting_reads", "flagged_artefact", "no_evidence", "low_depth"]
tab = N.groupby(["variant_patient", "evidence_class"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
left = np.zeros(len(tab))
for cls in order:
    b.barh(range(len(tab)), tab[cls], left=left, color=style[cls][1], height=0.6, edgecolor="white", linewidth=1)
    left += tab[cls].to_numpy()
b.set_yticks(range(len(tab)))
b.set_yticklabels(tab.index)
b.invert_yaxis()
b.set_xlabel("Candidate loci (n)")
b.xaxis.grid(True, color=GRID, lw=0.4)
b.set_axisbelow(True)
panel_label(b, "b", x=-0.25)
totals = [int((N["evidence_class"] == k).sum()) for k in order]
c.barh(range(len(order)), totals, color=[style[k][1] for k in order], height=0.62, edgecolor="white")
for i, v in enumerate(totals):
    c.text(v + max(totals) * 0.02, i, str(v), va="center", fontsize=6, color="#1a1a1a")
c.set_yticks(range(len(order)))
c.set_yticklabels([style[k][0].split(" (")[0] for k in order], fontsize=5.5)
c.invert_yaxis()
c.set_xlabel("Sites (n)")
c.set_xlim(0, max(totals) * 1.18)
c.xaxis.grid(True, color=GRID, lw=0.4)
c.set_axisbelow(True)
panel_label(c, "c", x=-0.45)
fig.legend(handles=[Patch(color=style[k][1], label=style[k][0]) for k in order],
           loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.01), fontsize=5.5, frameon=False)
fig.subplots_adjust(bottom=0.24)
save(fig, "Fig6_lookback_N", inputs=[lf, primary],
     note="Three panels. Dashed line is VAF 0.05. Zeros are drawn at 0.001 so they appear on the log axis. Confident means q<0.05 and VAF>=0.05 and is an operational class. PASS primary variants only.")
