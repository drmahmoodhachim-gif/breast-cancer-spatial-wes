#!/usr/bin/env python3
"""Figure 5 — Allele-specific copy number in junction and tumour (FACETS).

a  Genome-wide heatmap, one row per T/J sample (grouped by patient): log2(total CN / ploidy),
   blue = loss, red = gain, grey midpoint = neutral. Black tick under a segment = LOH (minor CN 0).
b  FACETS purity and the reconciled purity used for CCF, with ploidy.
Inputs: 07_cnv_purity/facets/*/*.segments.tsv, 07_cnv_purity/purity_ploidy.tsv,
        07_cnv_purity/purity_reconciled.tsv
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import GRID, INK2, PROJECT_DIR, SPECIMEN_COLOR, WIDTH, panel_label, save, setup  # noqa: E402

setup()
HG38 = {"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
        "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
        "19": 58617616, "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895}
offset, acc = {}, 0
for c, L in HG38.items():
    offset[c] = acc
    acc += L
GENOME = acc

files = list((PROJECT_DIR / "07_cnv_purity" / "facets").glob("*/*.segments.tsv"))
ppf = PROJECT_DIR / "07_cnv_purity" / "purity_ploidy.tsv"
recf = PROJECT_DIR / "07_cnv_purity" / "purity_reconciled.tsv"
S = pd.concat([pd.read_csv(f, sep="\t") for f in files])
S = S.rename(columns={"tcn.em": "tcn", "lcn.em": "lcn"})
S["chrom"] = S["chrom"].astype(str).str.replace("chr", "", regex=False)
PP = pd.read_csv(ppf, sep="\t").set_index("sample_id")
REC = pd.read_csv(recf, sep="\t").set_index("sample") if recf.exists() else None
samples = sorted(S["sample_id"].unique(), key=lambda s: (s.split("_")[0], {"J": 0, "T": 1}[s.split("_")[1]]))

cmap = LinearSegmentedColormap.from_list("cn", ["#2a78d6", "#efeeea", "#e34948"])
norm = Normalize(-1.2, 1.2)
fig = plt.figure(figsize=(WIDTH["double"], (22 + 6 * len(samples)) / 25.4))
gs = fig.add_gridspec(1, 2, width_ratios=[6, 1.3], wspace=0.06)
ax, axp = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

for i, sid in enumerate(samples):
    ploidy = PP["ploidy"].get(sid, 2.0)
    for r in S[S["sample_id"] == sid].itertuples():
        if r.chrom not in offset:
            continue
        x0 = offset[r.chrom] + r.start
        w = r.end - r.start
        val = np.log2(max(r.tcn, 0.1) / ploidy)
        ax.add_patch(plt.Rectangle((x0, i - 0.38), w, 0.62, color=cmap(norm(val)), lw=0))
        if pd.notna(r.lcn) and r.lcn == 0:
            ax.add_patch(plt.Rectangle((x0, i + 0.27), w, 0.12, color="#1a1a1a", lw=0))
for c in HG38:
    ax.axvline(offset[c], color="white", lw=0.8)
ax.set_xlim(0, GENOME)
ax.set_ylim(len(samples) - 0.5, -0.5)
ax.set_xticks([offset[c] + HG38[c] / 2 for c in HG38])
ax.set_xticklabels([c if c not in {"19", "21"} else "" for c in HG38], fontsize=5.5)
ax.set_yticks(range(len(samples)))
ax.set_yticklabels(samples)
for t, sid in zip(ax.get_yticklabels(), samples):
    t.set_color(SPECIMEN_COLOR[sid.split("_")[1]])
for i in range(2, len(samples), 2):
    ax.axhline(i - 0.5, color=INK2, lw=0.4)
ax.tick_params(length=0)
for s in ("left", "bottom"):
    ax.spines[s].set_visible(False)
ax.set_xlabel("Chromosome")
panel_label(ax, "a", x=-0.09, y=1.01)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation="horizontal",
                  fraction=0.04, pad=0.12 + 0.02 * (10 / max(len(samples), 1)), aspect=40)
cb.set_label("log2(total CN / ploidy)   ·   black bar = LOH", fontsize=6)
cb.outline.set_visible(False)

facets_p = [REC["facets_purity"].get(s, np.nan) if REC is not None else PP["purity"].get(s, np.nan) for s in samples]
used_p = [REC["used_purity"].get(s, np.nan) if REC is not None else facets_p[i] for i, s in enumerate(samples)]
y = np.arange(len(samples))
axp.barh(y - 0.16, facets_p, color="#c9c8c2", height=0.3, label="FACETS")
axp.barh(y + 0.16, used_p, color=[SPECIMEN_COLOR[s.split("_")[1]] for s in samples], height=0.3, label="Reconciled")
for i, s in enumerate(samples):
    shown = used_p[i] if used_p[i] == used_p[i] else facets_p[i]
    lbl = f"{PP['ploidy'].get(s, np.nan):.1f}n" if shown == shown else "not estimable (flat)"
    axp.text(min(shown if shown == shown else 0, 1) + 0.03, i, lbl,
             va="center", fontsize=5.5, color=INK2)
axp.legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=2, fontsize=5, frameon=False, borderaxespad=0)
axp.set_ylim(len(samples) - 0.5, -0.5)
axp.set_yticks([])
axp.set_xlim(0, 1.25)
axp.set_xticks([0, 0.5, 1])
axp.set_xlabel("Purity (ploidy)")
axp.xaxis.grid(True, color=GRID, lw=0.4)
axp.set_axisbelow(True)
panel_label(axp, "b", x=-0.1, y=1.01)
save(fig, "Fig5_CNV_purity", inputs=files + [ppf, recf])
