#!/usr/bin/env python3
"""Figure 4 — Cancer cell fraction (CCF) of each variant in tumour vs junction, per patient.

One panel per patient (small multiples, shared axes). Points coloured by spatial class;
driver-gene variants labelled. Error bars = 95 % CI (Wilson on VAF, propagated).
Falls back to VAF axes for a patient without purity/CN (stated on the panel).
Input: 09_ccf_clonality/ccf.tsv (or 06_variants/somatic_filtered.tsv)
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import CLASS_COLOR, CLASS_NAME, GRID, INK2, PROJECT_DIR, WIDTH, panel_label, save, setup  # noqa: E402

setup()
cf = PROJECT_DIR / "09_ccf_clonality" / "ccf.tsv"
src = cf if cf.exists() else PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
D = pd.read_csv(src, sep="\t")
patients = sorted(D["patient"].unique())
order = ["shared", "T_private", "J_private", "T_unresolved", "J_unresolved"]

fig, axes = plt.subplots(1, len(patients), figsize=(WIDTH["double"], 58 / 25.4), sharex=True, sharey=True,
                         gridspec_kw={"wspace": 0.12})
axes = [axes] if len(patients) == 1 else list(axes)
for k, (ax, p) in enumerate(zip(axes, patients)):
    sub = D[D["patient"] == p]
    use_ccf = "T_ccf_raw" in sub and sub[["T_ccf_raw", "J_ccf_raw"]].notna().all(axis=1).mean() > 0.5
    xcol, ycol = ("T_ccf_raw", "J_ccf_raw") if use_ccf else ("T_vaf", "J_vaf")
    ax.plot([0, 1], [0, 1], color=GRID, lw=0.6, zorder=0)
    for c in order:
        s = sub[sub["spatial_class"] == c]
        if s.empty:
            continue
        sx = s[xcol].astype(float)
        sy = s[ycol].astype(float)
        over = (sx > 1) | (sy > 1)
        if use_ccf:
            keep = s.loc[~over.to_numpy()]
            if len(keep):
                ax.errorbar(keep[xcol], keep[ycol],
                            xerr=[(keep[xcol] - keep["T_ccf_lo"]).clip(lower=0), (keep["T_ccf_hi"] - keep[xcol]).clip(lower=0)],
                            yerr=[(keep[ycol] - keep["J_ccf_lo"]).clip(lower=0), (keep["J_ccf_hi"] - keep[ycol]).clip(lower=0)],
                            fmt="none", ecolor=CLASS_COLOR[c], alpha=0.12, lw=0.3, zorder=1)
        ax.scatter(sx[~over].fillna(0).clip(lower=0), sy[~over].fillna(0).clip(lower=0),
                   s=9, color=CLASS_COLOR[c], edgecolor="white", linewidth=0.4, zorder=2)
        if over.any():
            ax.scatter(sx[over].clip(upper=1.02), sy[over].clip(upper=1.02),
                       s=16, marker="^", color=CLASS_COLOR[c], edgecolor="white", linewidth=0.4, zorder=3)
    force = sub["gene"].isin(["PIK3CA", "TP53", "BRCA2", "KMT2C"]) & sub["is_coding"].isin([True, "True", "true"])
    flagged = sub["is_driver_gene"].isin([True, "True", "true"]) & sub["is_coding"].isin([True, "True", "true"])
    lab = sub[force | flagged].copy()
    kmt = lab[lab["gene"] == "KMT2C"]
    lab = pd.concat([lab[lab["gene"] != "KMT2C"].drop_duplicates("gene"), kmt], ignore_index=True)
    lab["_m"] = lab[[xcol, ycol]].fillna(0).max(axis=1)
    lab = lab.sort_values("_m", ascending=False)
    for k_lab, r in enumerate(lab.itertuples()):
        x, y = getattr(r, xcol), getattr(r, ycol)
        x = 0 if pd.isna(x) else min(float(x), 1.02)
        y = 0 if pd.isna(y) else min(float(y), 1.02)
        right = x > 0.6
        text = r.gene
        if r.gene == "KMT2C" and pd.notna(getattr(r, "hgvsp_short", None)) and str(r.hgvsp_short):
            text = f"KMT2C {r.hgvsp_short}"
        if y > 0.72:
            xy = (-16 if right else 12, -16 - 6 * (k_lab % 2))
        else:
            xy = (-16 if right else 12, 8)
        print(f"[fig4] {r.patient} {text} {getattr(r, 'var_id', '')} x={x:.3f} y={y:.3f}")
        ax.annotate(text, (x, y),
                    xytext=xy,
                    textcoords="offset points", ha="right" if right else "left",
                    fontsize=5.5, fontstyle="italic", color="#1a1a1a",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#52514e"),
                    annotation_clip=True,
                    bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(-0.03, 1.05)
    ax.set_ylim(-0.03, 1.05)
    title = p + ("" if use_ccf else "  (VAF; no CN)")
    if p == "DH00024":
        title = "DH00024  not evaluable"
    ax.set_title(title, fontsize=6.5, color="#1a1a1a" if p != "DH00024" else INK2)
    ax.set_xlabel("Tumour CCF" if use_ccf else "Tumour VAF")
    if k == 0:
        ax.set_ylabel("Junction CCF" if use_ccf else "Junction VAF")
        panel_label(ax, "a", x=-0.35)
    ax.set_aspect("equal")

n_over_t = int((D["T_ccf_raw"] > 1).sum()) if "T_ccf_raw" in D else 0
n_over_j = int((D["J_ccf_raw"] > 1).sum()) if "J_ccf_raw" in D else 0
max_t = float(D["T_ccf_raw"].max()) if "T_ccf_raw" in D else float("nan")
max_j = float(D["J_ccf_raw"].max()) if "J_ccf_raw" in D else float("nan")
fig.subplots_adjust(bottom=0.28)
fig.text(0.5, 0.02,
         "Points are the uncapped estimate. A triangle on the axis edge is above 1 "
         f"({n_over_t} tumour, {n_over_j} junction; maxima {max_t:.2f} and {max_j:.2f}). "
         "Whiskers are the stored Wilson interval and remain capped at 1. "
         "A value above 1 is a mismatch among purity, copy number and VAF, not a cell fraction.",
         ha="center", fontsize=5.5, color=INK2)
handles = [Line2D([], [], marker="o", ls="", color=CLASS_COLOR[c], ms=4, label=CLASS_NAME[c])
           for c in order[:4]]
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.08))
save(fig, "Fig4_CCF_T_vs_J", inputs=[src],
     note="Points are uncapped T_ccf_raw and J_ccf_raw. Triangles mark estimates above 1, drawn at the axis edge. Stored Wilson intervals stay capped at 1.")
