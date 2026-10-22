#!/usr/bin/env python3
"""Figure — immune context across adjacent normal, junction and tumour.

a  T-cell fraction per sample (T cell ExTRECT, 95 % CI) for N, J, T of each patient.
b  Candidate neoantigens, one bar per patient, stacked by trunk / J branch / T branch.
   DH00023 uses the lineage table. Wild-type binders below 500 nM are excluded.
   Labelled candidate because there is no RNA expression data.
c  HLA class I loss of heterozygosity in J and T (FACETS), per gene.
Input: 14_immune/immune_summary.tsv, 14_immune/hla/hla_loh.tsv
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import (CLASS_COLOR, GRID, INK2, PROJECT_DIR, SPECIMEN_COLOR, WIDTH, panel_label,  # noqa: E402
                       placeholder, save, setup)

setup()
I = PROJECT_DIR / "14_immune"
S = pd.read_csv(I / "immune_summary.tsv", sep="\t")
loh_f = I / "hla" / "hla_loh.tsv"
LOH = pd.read_csv(loh_f, sep="\t") if loh_f.exists() else pd.DataFrame()
pts = list(dict.fromkeys(S["patient"]))

fig = plt.figure(figsize=(WIDTH["double"], 70 / 25.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.2, 0.8], wspace=0.45)

a = fig.add_subplot(gs[0])
if "TCRA.tcell.fraction" in S and S["TCRA.tcell.fraction"].notna().any():
    for k, s in enumerate(["N", "J", "T"]):
        d = S[S["specimen"] == s].set_index("patient").reindex(pts)
        x = np.arange(len(pts)) + (k - 1) * 0.25
        y = d["TCRA.tcell.fraction"] * 100
        lo = (d["TCRA.tcell.fraction"] - d["TCRA.tcell.fraction.lwr"]).clip(lower=0) * 100
        hi = (d["TCRA.tcell.fraction.upr"] - d["TCRA.tcell.fraction"]).clip(lower=0) * 100
        detected = y.fillna(0) > 0
        a.errorbar(x[detected.to_numpy()], y[detected.to_numpy()],
                   yerr=[lo[detected.to_numpy()], hi[detected.to_numpy()]],
                   fmt="o", ms=4, color=SPECIMEN_COLOR[s], capsize=1.5, lw=0.7, label=s)
        if (~detected).any():
            a.scatter(x[~detected.to_numpy()], np.zeros((~detected).sum()),
                      marker="o", s=28, facecolors="none", edgecolors=SPECIMEN_COLOR[s], lw=0.8, zorder=3)
        bad = d["qcFit"] > 4 if "qcFit" in d else pd.Series(False, index=d.index)
        if bad.any():
            a.scatter(x[bad.to_numpy()], y[bad.to_numpy()].fillna(0),
                      marker="D", s=28, facecolors="none", edgecolors=INK2, lw=0.7, zorder=4)
    a.set_xticks(range(len(pts)))
    a.set_xticklabels([p.replace("DH000", "DH") for p in pts])
    a.set_ylabel("T-cell fraction (%)")
    a.set_ylim(bottom=0)
    a.yaxis.grid(True, color=GRID, lw=0.4)
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=SPECIMEN_COLOR[s], ms=4, label=s) for s in ("N", "J", "T")
    ]
    handles.append(plt.Line2D([], [], marker="o", ls="", ms=5, markerfacecolor="none",
                              markeredgecolor="#4D4D4D", label="not detected"))
    handles.append(plt.Line2D([], [], marker="D", ls="", ms=5, markerfacecolor="none",
                              markeredgecolor=INK2, label="qcFit > 4"))
    a.legend(handles=handles, ncol=3, loc="upper center", fontsize=5, handletextpad=0.3,
             columnspacing=0.6, bbox_to_anchor=(0.5, 1.28), frameon=False)
else:
    placeholder(a, "Awaiting T cell ExTRECT")
panel_label(a, "a", x=-0.22)

b = fig.add_subplot(gs[1])
counts_path = I / "candidate_neoantigen_counts.tsv"
if counts_path.exists():
    C = pd.read_csv(counts_path, sep="\t").set_index("patient").reindex(pts)
    y = np.arange(len(pts))
    left = np.zeros(len(pts))
    stacks = [
        ("trunk", CLASS_COLOR["shared"], "Trunk"),
        ("J_branch", SPECIMEN_COLOR["J"], "J branch"),
        ("T_branch", SPECIMEN_COLOR["T"], "T branch"),
        ("unassigned", "#bdbdbd", "Unassigned"),
    ]
    for col, color, lab in stacks:
        vals = C[col].fillna(0).to_numpy()
        if vals.sum() == 0 and col == "unassigned":
            continue
        b.barh(y, vals, left=left, color=color, height=0.62, edgecolor="white", label=lab)
        left = left + vals
    b.set_yticks(y)
    b.set_yticklabels([p.replace("DH000", "DH") for p in pts], fontsize=6)
    b.invert_yaxis()
    b.set_xlabel("Candidate neoantigens (no RNA)")
    b.xaxis.grid(True, color=GRID, lw=0.4)
    b.set_axisbelow(True)
    b.legend(fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4,
             columnspacing=0.6, handlelength=1, frameon=False)
else:
    placeholder(b, "Awaiting HLA typing + pVACseq")
panel_label(b, "b", x=-0.3)

c = fig.add_subplot(gs[2])
if len(LOH):
    order = sorted(LOH["sample_id"].unique())
    genes = ["HLA-A", "HLA-B", "HLA-C"]
    cmap = {"LOH": "#4D4D4D", "allelic imbalance": "#d8a0b8", "balanced": "#efeeea",
            "homozygous (not informative)": "#c9c8c2", "no segment": "#ffffff"}
    for i, sid in enumerate(order):
        for j, g in enumerate(genes):
            st = LOH[(LOH["sample_id"] == sid) & (LOH["gene"] == g)]["state"]
            st = st.iloc[0] if len(st) else "no segment"
            c.add_patch(plt.Rectangle((j, i), 0.92, 0.88, color=cmap.get(st, "#ffffff"), lw=0))
    c.set_xlim(0, 3)
    c.set_ylim(len(order), 0)
    c.set_xticks([0.46, 1.46, 2.46])
    c.set_xticklabels(["A", "B", "C"])
    c.set_yticks([i + 0.44 for i in range(len(order))])
    c.set_yticklabels([s.replace("DH000", "DH") for s in order], fontsize=5.5)
    c.tick_params(length=0)
    for sp in c.spines.values():
        sp.set_visible(False)
    c.set_title("HLA class I", fontsize=6.5)
    c.legend(handles=[Patch(color=v, label=k.split(" (")[0]) for k, v in cmap.items() if k != "no segment"],
             fontsize=5, loc="upper left", bbox_to_anchor=(1.02, 1))
else:
    placeholder(c, "Awaiting FACETS\n(HLA LOH)")
panel_label(c, "c", x=-0.35)
save(fig, "FigS_immune", kind="supplementary", inputs=[I / "immune_summary.tsv", I / "candidate_neoantigen_counts.tsv"],
     note="All T-cell points are the unadjusted TCRA fraction. Open circles are below the TCRA detection limit, not a measured fraction of 0. qcFit > 4 is the low-confidence flag already used for this assay. DH00024 junction and tumour are unchanged by copy-number adjustment. Whiskers are the TcellExTRECT interval.")
