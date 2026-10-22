"""Shared publication style for all figures.

* Sizes in mm for journal columns (single 89 mm, 1.5 col 120 mm, double 183 mm).
* Sans-serif 7 pt body / 8 pt panel letters; editable text in PDF/SVG (fonttype 42).
* Every figure saved as vector PDF + SVG + 600 dpi PNG (TIFF optional) with a provenance file.
* Colours: one fixed colour per specimen everywhere in the paper (validated colour-blind safe,
  all-pairs CVD dE >= 14.9). Identity is never colour-alone: specimens are also labelled.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402,F401

MM = 1 / 25.4
WIDTH = {"single": 89 * MM, "onehalf": 120 * MM, "double": 183 * MM}

SPECIMEN_COLOR = {"N": "#0072B2", "J": "#E69F00", "T": "#882255"}
SPECIMEN_NAME = {"N": "Adjacent normal", "J": "Junction", "T": "Tumour"}
CLASS_COLOR = {"shared": "#3d3d3a", "T_private": SPECIMEN_COLOR["T"], "J_private": SPECIMEN_COLOR["J"],
               "T_unresolved": "#c9c8c2", "J_unresolved": "#c9c8c2"}
CLASS_NAME = {"shared": "Shared T+J", "T_private": "T-enriched", "J_private": "J-enriched",
              "T_unresolved": "Unresolved (low power)", "J_unresolved": "Unresolved (low power)"}
# Oncoplot consequence classes (validated all-pairs; teal relies on cell gaps + legend)
CONSEQ_COLOR = {"Missense": "#2a78d6", "Truncating": "#eb6834", "In-frame / other": "#1baf7a"}
INK, INK2, GRID = "#1a1a1a", "#52514e", "#e4e3df"


def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 600, "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })


def panel_label(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def save(fig, name: str, kind: str = "main", inputs=None, note: str = ""):
    """Save PDF/SVG/PNG into 11_figures/<kind>/ and write provenance."""
    out = PROJECT_DIR / "11_figures" / kind
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("pdf", "svg", "png"):
        p = out / f"{name}.{ext}"
        fig.savefig(p)
        paths.append(p)
    write_provenance(out / f"{name}.pdf", inputs or [], {"note": note})
    plt.close(fig)
    print(f"[figure] {kind}/{name} (.pdf .svg .png)")
    return paths


def placeholder(ax, msg: str):
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, color=INK2, fontsize=6, wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
