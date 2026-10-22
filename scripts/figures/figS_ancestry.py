#!/usr/bin/env python3
"""Supplementary ancestry: peddy PC1 vs PC2.

1000 Genomes superpopulations are faint clouds. Study samples are labelled by the
inferred group only. Heterozygosity and relatedness plots stay in 16_ancestry/.
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import INK2, PROJECT_DIR, WIDTH, save, setup  # noqa: E402

setup()
ped = pd.read_csv(PROJECT_DIR / "16_ancestry" / "peddy.peddy.ped", sep="\t")
bg_path = PROJECT_DIR / "16_ancestry" / "peddy.background_pca.json"
raw = json.loads(bg_path.read_text())
clouds = []
if isinstance(raw, dict):
    items = raw.items()
elif isinstance(raw, list) and raw and isinstance(raw[0], dict):
    items = []
    for row in raw:
        pop = row.get("super_pop") or row.get("population") or row.get("ancestry") or "1000G"
        items.append((pop, row))
else:
    items = []
by_pop = {}
for pop, row in items:
    if isinstance(row, dict) and "PC1" in row and "PC2" in row:
        by_pop.setdefault(str(pop), []).append((row["PC1"], row["PC2"]))
    elif isinstance(row, (list, tuple)) and len(row) >= 2:
        by_pop.setdefault(str(pop), []).append((row[0], row[1]))

fig, ax = plt.subplots(figsize=(WIDTH["single"], WIDTH["single"]))
pop_col = {
    "AFR": "#6B4C3B", "EUR": "#8A8A82", "EAS": "#3D6B66", "SAS": "#6B6B3A",
    "AMR": "#6B5478", "1000G": "#bdbdbd",
}
for pop, pts in by_pop.items():
    xs, ys = zip(*pts)
    ax.scatter(xs, ys, s=6, alpha=0.35, linewidths=0, color=pop_col.get(pop, "#9a9a94"), label=pop, zorder=1)
ax.scatter(ped["PC1"], ped["PC2"], s=22, c="#1a1a1a", zorder=3, label="This study")
labelled = set()
for group, sub in ped.groupby("ancestry-prediction"):
    x, y = float(sub["PC1"].mean()), float(sub["PC2"].mean())
    text = f"{group} (n = {len(sub)})" if len(sub) > 1 else str(group)
    ax.annotate(text, (x, y), textcoords="offset points", xytext=(6, 4), fontsize=6, color=INK2)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.legend(loc="best", fontsize=5.5, frameon=False, markerscale=1.4)
save(fig, "FigS_ancestry", kind="supplementary", inputs=[bg_path, PROJECT_DIR / "16_ancestry" / "peddy.peddy.ped"],
     note="D046. Labels are the inferred superpopulation only. Het and relatedness plots are not copied here.")
