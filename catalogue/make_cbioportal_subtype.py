"""PAM50 subtype of TCGA breast tumours with or without a mutation in each gene.

Uses the open cBioPortal files already retrieved for brca_tcga_pan_can_atlas_2018.
Patients are the unit. A patient is mutated if any sequenced sample carries the call.
Subtype is the patient-level SUBTYPE attribute (BRCA_LumA, BRCA_LumB, BRCA_Basal, BRCA_Her2, BRCA_Normal).
Luminal here means LumA or LumB. This is an association in bulk tumours, not a predicted consequence.
"""
import json
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

TEMP = Path(r"C:\Users\chime\AppData\Local\Temp")
OUT = Path(__file__).resolve().parent
sequenced = set(json.loads((TEMP / "cbio_seq.json").read_text(encoding="utf-8"))["sampleIds"])
mutations = json.loads((TEMP / "cbio_genes.json").read_text(encoding="utf-8"))
subtypes = {
    r["patientId"]: r["value"]
    for r in json.loads((TEMP / "cbio_subtype.json").read_text(encoding="utf-8"))
}
ORDER = ["BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal", "BRCA_Normal"]
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Splice_Region"}

def patient_of(sample_id):
    parts = sample_id.split("-")
    return "-".join(parts[:3])

sequenced_patients = {patient_of(s) for s in sequenced}
typed = {p for p in sequenced_patients if p in subtypes}

gene_mut = defaultdict(set)
kmt_trunc = set()
kmt_missense = set()
for r in mutations:
    if r.get("sampleId") not in sequenced:
        continue
    g = (r.get("gene") or {}).get("hugoGeneSymbol")
    p = patient_of(r["sampleId"])
    gene_mut[g].add(p)
    if g == "KMT2C":
        if r.get("mutationType") in TRUNC:
            kmt_trunc.add(p)
        elif r.get("mutationType") == "Missense_Mutation":
            kmt_missense.add(p)

def counts(mut_patients):
    mut = typed & mut_patients
    wt = typed - mut_patients
    def bucket(ids):
        c = {s: 0 for s in ORDER}
        for p in ids:
            c[subtypes[p]] += 1
        return c
    return bucket(mut), bucket(wt), len(mut), len(wt)

def luminal_p(mut_c, wt_c):
    mut_lum = mut_c["BRCA_LumA"] + mut_c["BRCA_LumB"]
    wt_lum = wt_c["BRCA_LumA"] + wt_c["BRCA_LumB"]
    mut_n = sum(mut_c.values())
    wt_n = sum(wt_c.values())
    table = [[mut_lum, mut_n - mut_lum], [wt_lum, wt_n - wt_lum]]
    oddsratio, p = fisher_exact(table)
    return mut_lum, mut_n, wt_lum, wt_n, oddsratio, p

groups = [
    ("PIK3CA", "any mutation", gene_mut["PIK3CA"]),
    ("TP53", "any mutation", gene_mut["TP53"]),
    ("KMT2C", "any mutation", gene_mut["KMT2C"]),
    ("KMT2C", "truncating or splice", kmt_trunc),
    ("KMT2C", "missense and not truncating", kmt_missense - kmt_trunc),
    ("TBX3", "any mutation", gene_mut["TBX3"]),
    ("NF1", "any mutation", gene_mut["NF1"]),
    ("BRCA2", "any mutation", gene_mut["BRCA2"]),
]

header = ["gene", "contrast", "group", *ORDER, "n", "luminal_n", "luminal_fraction", "fisher_or_luminal_vs_rest", "fisher_p"]
lines = ["\t".join(header)]
summary = []
for gene, contrast, mut_ids in groups:
    mut_c, wt_c, n_mut, n_wt = counts(mut_ids)
    mut_lum, mut_n, wt_lum, wt_n, oddsratio, p = luminal_p(mut_c, wt_c)
    for group, c, n, lum in (
        ("mutated", mut_c, n_mut, mut_lum),
        ("no mutation in this contrast", wt_c, n_wt, wt_lum),
    ):
        frac = f"{lum / n:.3f}" if n else ""
        or_s = f"{oddsratio:.3f}" if group == "mutated" else ""
        p_s = f"{p:.4g}" if group == "mutated" else ""
        lines.append("\t".join([
            gene, contrast, group,
            *[str(c[s]) for s in ORDER],
            str(n), str(lum), frac, or_s, p_s,
        ]))
    summary.append((gene, contrast, n_mut, mut_lum / n_mut, n_wt, wt_lum / n_wt, oddsratio, p))

path = OUT / "tables" / "table_cbioportal_subtype.tsv"
path.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("typed_and_sequenced", len(typed))
for row in summary:
    print("\t".join(str(x if not isinstance(x, float) else f"{x:.4g}") for x in row))
print("wrote", path)
