"""Catalogue comparison against the open cBioPortal TCGA breast mutation list.

Study: brca_tcga_pan_can_atlas_2018
Profile: brca_tcga_pan_can_atlas_2018_mutations
Denominator: sample list brca_tcga_pan_can_atlas_2018_sequenced
Retrieved from https://www.cbioportal.org/api on 16 October 2026.
This is a bulk tumour catalogue. It has no junction or adjacent-normal sample.
"""
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEMP = Path(r"C:\Users\chime\AppData\Local\Temp")
OUT = Path(__file__).resolve().parent
seq = json.loads((TEMP / "cbio_seq.json").read_text(encoding="utf-8"))
sequenced = set(seq["sampleIds"])
n_seq = len(sequenced)
rows = json.loads((TEMP / "cbio_genes.json").read_text(encoding="utf-8"))
rows = [r for r in rows if r.get("sampleId") in sequenced]

INK = "#1a1a1a"
GRID = "#e4e3df"
MISS = "#2a78d6"
TRUNC = "#eb6834"
ZERO = "#bdbdbd"

def gene(r):
    return (r.get("gene") or {}).get("hugoGeneSymbol")

by_samples = defaultdict(set)
classes = defaultdict(Counter)
protein_samples = defaultdict(set)
for r in rows:
    g = gene(r)
    sid = r["sampleId"]
    by_samples[g].add(sid)
    classes[g][r.get("mutationType")] += 1
    protein_samples[(g, r.get("proteinChange") or "")].add(sid)

# Exact allele search by genomic position (hg19 coordinates are what cBioPortal stores
# for this TCGA profile; the cohort table is hg38. Positions are checked as protein change.)
cohort = [
    ("DH00019", "PIK3CA", "p.Glu545Lys", "E545K", "shared", "Hotspot missense"),
    ("DH00023", "PIK3CA", "p.His1047Arg", "H1047R", "shared", "Hotspot missense"),
    ("DH00022", "TP53", "p.Arg248Leu", "R248L", "shared", "Same codon as a hotspot; this substitution was absent"),
    ("DH00023", "KMT2C", "p.Gln2197Ter", "Q2197*", "shared, J-lineage", "Truncating class is recurrent; this allele was absent"),
    ("DH00023", "KMT2C", "p.Leu4553Arg", "L4553R", "T-enriched", "Missense in FYRN; this allele was absent"),
    ("DH00019", "TBX3", "p.Thr210AsnfsTer30", "T210Nfs*30", "shared", "TBX3 frameshift class is present; this allele was absent"),
    ("DH00019", "NF1", "splice donor c.7286_7321+5del", "", "shared", "NF1 splice-site mutations are present; this deletion was not matched by protein change"),
    ("DH00022", "BRCA2", "p.Ser795Leu", "S795L", "shared", "Variant of uncertain significance; this allele was absent"),
]

lines = [
    "patient\tgene\tprotein_change\tspatial_class\ttcga_tumours_with_exact_change\ttcga_tumours_with_any_mutation\ttcga_sequenced\tnote"
]
for patient, g, protein, short, spatial, note in cohort:
    exact = len(protein_samples[(g, short)]) if short else 0
    lines.append(
        f"{patient}\t{g}\t{protein}\t{spatial}\t{exact}\t{len(by_samples[g])}\t{n_seq}\t{note}"
    )
(OUT / "tables" / "table_cbioportal_tcga_brca.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")

kmt_order = [
    ("Missense", classes["KMT2C"]["Missense_Mutation"]),
    ("Nonsense", classes["KMT2C"]["Nonsense_Mutation"]),
    ("Frameshift", classes["KMT2C"]["Frame_Shift_Del"] + classes["KMT2C"]["Frame_Shift_Ins"]),
    ("Splice", classes["KMT2C"]["Splice_Site"] + classes["KMT2C"]["Splice_Region"]),
    ("In-frame", classes["KMT2C"]["In_Frame_Del"] + classes["KMT2C"]["In_Frame_Ins"]),
]
kmt_lines = ["mutation_class\tn_mutations\tn_tumours_with_any_KMT2C\tn_sequenced"]
for name, n in kmt_order:
    kmt_lines.append(f"{name}\t{n}\t{len(by_samples['KMT2C'])}\t{n_seq}")
(OUT / "tables" / "table_cbioportal_kmt2c_classes.tsv").write_text("\n".join(kmt_lines) + "\n", encoding="utf-8")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 7,
    "axes.linewidth": 0.4,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), gridspec_kw={"width_ratios": [1.35, 1]})

plot_rows = [
    ("PIK3CA p.His1047Arg", 129, MISS),
    ("PIK3CA p.Glu545Lys", 65, MISS),
    ("TP53 p.Arg248Trp", 5, MISS),
    ("TP53 p.Arg248Gln", 4, MISS),
    ("TP53 p.Arg248Leu", 0, ZERO),
    ("KMT2C p.Gln2197Ter", 0, ZERO),
    ("KMT2C p.Leu4553Arg", 0, ZERO),
    ("TBX3 p.Thr210AsnfsTer30", 0, ZERO),
    ("BRCA2 p.Ser795Leu", 0, ZERO),
]
axes[0].barh(
    range(len(plot_rows)),
    [n for _, n, _ in plot_rows][::-1],
    color=[c for _, _, c in plot_rows][::-1],
    height=0.62,
)
axes[0].set_yticks(range(len(plot_rows)))
axes[0].set_yticklabels([lab for lab, _, _ in plot_rows][::-1])
axes[0].set_xlabel(f"Tumours in TCGA breast (n = {n_seq})")
axes[0].set_xlim(0, 160)
axes[0].xaxis.grid(True, color=GRID, lw=0.4)
axes[0].set_axisbelow(True)
axes[0].text(0.0, 1.04, "a", transform=axes[0].transAxes, fontweight="bold", fontsize=8, va="bottom")
for i, (_, n, _) in enumerate(plot_rows[::-1]):
    axes[0].text(n + 2, i, str(n), va="center", fontsize=6, color=INK)

axes[1].bar(
    range(len(kmt_order)),
    [n for _, n in kmt_order],
    color=[MISS, TRUNC, TRUNC, TRUNC, "#8E99A4"],
    width=0.72,
)
axes[1].set_xticks(range(len(kmt_order)))
axes[1].set_xticklabels([name for name, _ in kmt_order], rotation=30, ha="right")
axes[1].set_ylabel("KMT2C mutations")
axes[1].set_ylim(0, 55)
axes[1].yaxis.grid(True, color=GRID, lw=0.4)
axes[1].set_axisbelow(True)
axes[1].text(0.0, 1.04, "b", transform=axes[1].transAxes, fontweight="bold", fontsize=8, va="bottom")
for i, (_, n) in enumerate(kmt_order):
    axes[1].text(i, n + 1.2, str(n), ha="center", fontsize=6, color=INK)

fig.tight_layout()
fig.savefig(OUT / "figures" / "FigS_catalogue_tcga.pdf")
fig.savefig(OUT / "figures" / "FigS_catalogue_tcga.png", dpi=200)
print("sequenced", n_seq)
print("PIK3CA", len(by_samples["PIK3CA"]), "TP53", len(by_samples["TP53"]), "KMT2C", len(by_samples["KMT2C"]))
print("wrote tables and figure")
