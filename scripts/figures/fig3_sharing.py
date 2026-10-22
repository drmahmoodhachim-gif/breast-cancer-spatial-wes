#!/usr/bin/env python3
"""Figure 3 — Spatial sharing of somatic variants between tumour (T) and junction (J).

a  Counts per patient: shared / T-private / J-private / unresolved (insufficient power).
b  Same as proportions, with tier-1 (Mutect2 + Strelka2) fraction marked.
Private = present in one, ABSENT (>= 80 % power to detect half the other sample's VAF) in the other.
Input: 06_variants/somatic_filtered.tsv
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import CLASS_COLOR, CLASS_NAME, GRID, PROJECT_DIR, WIDTH, panel_label, save, setup  # noqa: E402

setup()
vf = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
V = pd.read_csv(vf, sep="\t")
V["cls"] = V["spatial_class"].replace({"J_unresolved": "unresolved", "T_unresolved": "unresolved"})
order = ["shared", "T_private", "J_private", "unresolved"]
color = {**CLASS_COLOR, "unresolved": "#c9c8c2"}
name = {**CLASS_NAME, "unresolved": "Unresolved (low power)"}
tab = V.groupby(["patient", "cls"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
patients = list(tab.index)

fig, axes = plt.subplots(1, 2, figsize=(WIDTH["double"], 62 / 25.4), gridspec_kw={"wspace": 0.45})
resolved = ["shared", "T_private", "J_private"]
for ax, prop, letter in ((axes[0], False, "a"), (axes[1], True, "b")):
    if prop:
        den = tab[resolved].sum(axis=1).replace(0, np.nan)
        data = tab[resolved].div(den, axis=0) * 100
        data["unresolved"] = 0.0
        if "DH00024" in data.index:
            data.loc["DH00024", order] = 0.0
        data = data[order]
    else:
        data = tab
    left = np.zeros(len(patients))
    for c in order:
        vals = data[c].to_numpy()
        hatches = []
        for i, patient in enumerate(patients):
            hatches.append("////" if (prop and patient == "DH00024") or c == "unresolved" else None)
        for i, (l, v) in enumerate(zip(left, vals)):
            ax.barh(i, v, left=l, color=color[c], height=0.62, edgecolor="white", linewidth=1.0,
                    hatch=hatches[i])
        for i, (l, v) in enumerate(zip(left, vals)):
            if (prop and v >= 9) or (not prop and v >= 0.09 * data.sum(axis=1).max()):
                txt = f"{v:.0f}%" if prop else f"{int(v)}"
                ax.text(l + v / 2, i, txt, ha="center", va="center", fontsize=5.5,
                        color="white" if c in ("shared", "T_private") else "#1a1a1a")
        left += vals
    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels(patients)
    ax.invert_yaxis()
    if prop:
        ax.set_xlim(0, 100)
        if "DH00024" in patients:
            ax.text(50, patients.index("DH00024"), "not evaluable", va="center", ha="center",
                    fontsize=5.5, color="#1a1a1a",
                    bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.9))
    ax.set_xlabel("% of resolved variants" if prop else "Somatic variants (n)")
    ax.xaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    panel_label(ax, letter, x=-0.28)

fig.subplots_adjust(bottom=0.32)
fig.legend(handles=[Patch(facecolor=color[c], label=name[c], hatch="////" if c == "unresolved" else None,
                          edgecolor="white") for c in order] +
          [Patch(facecolor="#e4e3df", edgecolor="white", hatch="////", label="DH00024 not evaluable")],
           loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.01))
save(fig, "Fig3_TJ_sharing", inputs=[vf])
