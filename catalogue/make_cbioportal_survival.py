"""Mutated versus wild-type comparisons in open TCGA breast cancers.

Study brca_tcga_pan_can_atlas_2018. Patients with mutation data are the unit.
A patient is mutated if any sequenced sample carries the call.
Tests were specified before looking at the p-values: age, stage, fraction of the
genome altered, mutation count, nonsynonymous TMB, and log-rank tests of overall,
disease-specific, and progression-free survival. An overall-survival log-rank
restricted to LumA or LumB is the subtype check. Benjamini-Hochberg q-values are
computed across this table. This is not a model of the two KMT2C alleles that are
absent from the catalogue.
"""
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2, chi2_contingency, mannwhitneyu

TEMP = Path(r"C:\Users\chime\AppData\Local\Temp")
OUT = Path(__file__).resolve().parent
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Splice_Region"}
STAGE = {
    "STAGE I": "I", "STAGE IA": "I", "STAGE IB": "I",
    "STAGE II": "II", "STAGE IIA": "II", "STAGE IIB": "II",
    "STAGE III": "III", "STAGE IIIA": "III", "STAGE IIIB": "III", "STAGE IIIC": "III",
    "STAGE IV": "IV",
}


def patient_of(sample_id):
    return "-".join(sample_id.split("-")[:3])


def load_clinical(path, id_key):
    rows = json.loads(path.read_text(encoding="utf-8"))
    out = defaultdict(dict)
    for r in rows:
        out[r[id_key]][r["clinicalAttributeId"]] = r["value"]
    return out


def dead_flag(value, positive):
    if value is None:
        return None
    return 1 if value.startswith("1:") or value in positive else 0


def to_float(value):
    if value is None or value in {"", "NA", "[Not Available]"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


sequenced = set(json.loads((TEMP / "cbio_seq.json").read_text(encoding="utf-8"))["sampleIds"])
mutations = json.loads((TEMP / "cbio_genes.json").read_text(encoding="utf-8"))
patients = load_clinical(TEMP / "cbio_clin_patient.json", "patientId")
sample_clin = defaultdict(dict)
for name in ("cbio_fga.json", "cbio_mutcount.json", "cbio_tmb.json"):
    for r in json.loads((TEMP / name).read_text(encoding="utf-8")):
        sample_clin[r["sampleId"]][r["clinicalAttributeId"]] = r["value"]

sequenced_patients = sorted({patient_of(s) for s in sequenced})
primary = {}
for sid in sorted(sequenced):
    pid = patient_of(sid)
    if pid not in primary or sid.split("-")[3].startswith("01"):
        primary[pid] = sid

gene_mut = defaultdict(set)
kmt_trunc = set()
for r in mutations:
    if r.get("sampleId") not in sequenced:
        continue
    g = (r.get("gene") or {}).get("hugoGeneSymbol")
    pid = patient_of(r["sampleId"])
    gene_mut[g].add(pid)
    if g == "KMT2C" and r.get("mutationType") in TRUNC:
        kmt_trunc.add(pid)

contrasts = [
    ("PIK3CA", "any mutation", gene_mut["PIK3CA"]),
    ("TP53", "any mutation", gene_mut["TP53"]),
    ("KMT2C", "any mutation", gene_mut["KMT2C"]),
    ("KMT2C", "truncating or splice", kmt_trunc),
    ("TBX3", "any mutation", gene_mut["TBX3"]),
    ("NF1", "any mutation", gene_mut["NF1"]),
    ("BRCA2", "any mutation", gene_mut["BRCA2"]),
]


def split(mut_ids):
    base = [p for p in sequenced_patients]
    return [p for p in base if p in mut_ids], [p for p in base if p not in mut_ids]


def numeric_pair(ids, getter):
    return np.array([v for p in ids if (v := getter(p)) is not None], float)


def mw(mut, wt, getter):
    a, b = numeric_pair(mut, getter), numeric_pair(wt, getter)
    if len(a) < 5 or len(b) < 5:
        return None
    res = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_mut": len(a), "n_wt": len(b),
        "median_mut": float(np.median(a)), "median_wt": float(np.median(b)),
        "p": float(res.pvalue), "stat": f"median {np.median(a):.3g} vs {np.median(b):.3g}",
    }


def stage_of(pid):
    return STAGE.get(patients[pid].get("AJCC_PATHOLOGIC_TUMOR_STAGE", ""), None)


def stage_test(mut, wt):
    levels = ["I", "II", "III", "IV"]
    table = []
    for ids in (mut, wt):
        c = {k: 0 for k in levels}
        for p in ids:
            s = stage_of(p)
            if s:
                c[s] += 1
        table.append([c[k] for k in levels])
    arr = np.array(table)
    if (arr < 1).any() or arr.sum() < 20:
        return None
    chi, p, _, _ = chi2_contingency(arr)
    return {"n_mut": int(arr[0].sum()), "n_wt": int(arr[1].sum()), "median_mut": "", "median_wt": "",
            "p": float(p), "stat": "stage I/II/III/IV " + ",".join(map(str, arr[0])) + " vs " + ",".join(map(str, arr[1]))}


def survival_arrays(ids, month_key, status_key):
    t, e = [], []
    for p in ids:
        month = to_float(patients[p].get(month_key))
        flag = dead_flag(patients[p].get(status_key), set())
        if month is None or flag is None:
            continue
        t.append(month)
        e.append(flag)
    return np.array(t, float), np.array(e, int)


def logrank(t1, e1, t2, e2):
    times = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))
    if len(times) == 0:
        return None
    o_minus_e = 0.0
    var = 0.0
    for t in times:
        n1 = np.sum(t1 >= t)
        n2 = np.sum(t2 >= t)
        d1 = np.sum((t1 == t) & (e1 == 1))
        d2 = np.sum((t2 == t) & (e2 == 1))
        n, d = n1 + n2, d1 + d2
        if n <= 1 or d == 0:
            continue
        expected = n1 * d / n
        o_minus_e += d1 - expected
        var += n1 * n2 * d * (n - d) / (n * n * (n - 1))
    if var <= 0:
        return None
    stat = o_minus_e ** 2 / var
    p = float(chi2.sf(stat, 1))
    return p, int(e1.sum()), int(e2.sum()), len(t1), len(t2)


def km_curve(times, events):
    order = np.argsort(times)
    times, events = times[order], events[order]
    n = len(times)
    surv = 1.0
    xs, ys = [0.0], [1.0]
    at_risk = n
    i = 0
    uniq = np.unique(times)
    for t in uniq:
        d = c = 0
        while i < n and times[i] == t:
            if events[i] == 1:
                d += 1
            else:
                c += 1
            i += 1
        if at_risk and d:
            surv *= 1 - d / at_risk
        xs.append(float(t))
        ys.append(surv)
        at_risk -= d + c
    return np.array(xs), np.array(ys)


rows_out = []
curves = {}
tests = []


def add(gene, contrast, variable, result, curve=None):
    if result is None:
        return
    tests.append((gene, contrast, variable, result))
    if curve is not None:
        curves[(gene, contrast, variable)] = curve


for gene, contrast, mut_ids in contrasts:
    mut, wt = split(mut_ids)
    add(gene, contrast, "age_years", mw(mut, wt, lambda p: to_float(patients[p].get("AGE"))))
    add(gene, contrast, "stage", stage_test(mut, wt))
    add(gene, contrast, "fraction_genome_altered", mw(
        mut, wt, lambda p: to_float(sample_clin[primary[p]].get("FRACTION_GENOME_ALTERED"))))
    add(gene, contrast, "mutation_count", mw(
        mut, wt, lambda p: to_float(sample_clin[primary[p]].get("MUTATION_COUNT"))))
    add(gene, contrast, "tmb_nonsynonymous", mw(
        mut, wt, lambda p: to_float(sample_clin[primary[p]].get("TMB_NONSYNONYMOUS"))))
    for label, mk, sk in (
        ("overall_survival", "OS_MONTHS", "OS_STATUS"),
        ("disease_specific_survival", "DSS_MONTHS", "DSS_STATUS"),
        ("progression_free", "PFS_MONTHS", "PFS_STATUS"),
    ):
        t1, e1 = survival_arrays(mut, mk, sk)
        t2, e2 = survival_arrays(wt, mk, sk)
        lr = logrank(t1, e1, t2, e2)
        if lr is None:
            continue
        p, ev1, ev2, n1, n2 = lr
        add(gene, contrast, label, {
            "n_mut": n1, "n_wt": n2, "median_mut": ev1, "median_wt": ev2,
            "p": p, "stat": f"events {ev1}/{n1} vs {ev2}/{n2}",
        }, (t1, e1, t2, e2))
    lum = {p for p in sequenced_patients if patients[p].get("SUBTYPE") in {"BRCA_LumA", "BRCA_LumB"}}
    t1, e1 = survival_arrays([p for p in mut if p in lum], "OS_MONTHS", "OS_STATUS")
    t2, e2 = survival_arrays([p for p in wt if p in lum], "OS_MONTHS", "OS_STATUS")
    lr = logrank(t1, e1, t2, e2)
    if lr is not None:
        p, ev1, ev2, n1, n2 = lr
        add(gene, contrast, "overall_survival_luminal_only", {
            "n_mut": n1, "n_wt": n2, "median_mut": ev1, "median_wt": ev2,
            "p": p, "stat": f"events {ev1}/{n1} vs {ev2}/{n2}",
        }, (t1, e1, t2, e2))

ps = np.array([r[3]["p"] for r in tests], float)
order = np.argsort(ps)
q = np.empty(len(ps))
prev = 1.0
m = len(ps)
for rank in range(m, 0, -1):
    i = order[rank - 1]
    prev = min(prev, ps[i] * m / rank)
    q[i] = prev

lines = ["gene\tcontrast\tvariable\tn_mutated\tn_wild_type\tstatistic\tp\tq_bh"]
for (gene, contrast, variable, result), qq in zip(tests, q):
    lines.append("\t".join([
        gene, contrast, variable, str(result["n_mut"]), str(result["n_wt"]),
        result["stat"], f"{result['p']:.4g}", f"{qq:.4g}",
    ]))
(OUT / "tables" / "table_cbioportal_mut_vs_wt.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("tests", m)
for line in lines[1:]:
    p = float(line.split("\t")[-2])
    if p < 0.05:
        print(line)

# Survival figure: overall, all patients and luminal-only, three contrasts.
draw = [
    ("PIK3CA", "any mutation", "overall_survival", "PIK3CA, all"),
    ("TP53", "any mutation", "overall_survival", "TP53, all"),
    ("KMT2C", "truncating or splice", "overall_survival", "KMT2C truncating, all"),
    ("PIK3CA", "any mutation", "overall_survival_luminal_only", "PIK3CA, luminal"),
    ("TP53", "any mutation", "overall_survival_luminal_only", "TP53, luminal"),
    ("KMT2C", "truncating or splice", "overall_survival_luminal_only", "KMT2C truncating, luminal"),
]
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 7, "axes.linewidth": 0.4,
    "xtick.major.width": 0.4, "ytick.major.width": 0.4, "pdf.fonttype": 42,
})
fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.6), sharex=True, sharey=True)
lookup = {(g, c, v): r for (g, c, v), r in zip(
    [(a, b, c) for a, b, c, _ in ((t[0], t[1], t[2], None) for t in tests)],
    [t[3] for t in tests],
)}
# rebuild lookup properly
lookup = {}
for gene, contrast, variable, result in tests:
    lookup[(gene, contrast, variable)] = result

for ax, (gene, contrast, variable, title) in zip(axes.ravel(), draw):
    t1, e1, t2, e2 = curves[(gene, contrast, variable)]
    x1, y1 = km_curve(t1, e1)
    x2, y2 = km_curve(t2, e2)
    res = lookup[(gene, contrast, variable)]
    ax.step(x1, y1, where="post", color="#2a78d6", lw=1.0, label=f"Mutated {res['n_mut']}")
    ax.step(x2, y2, where="post", color="#4D4D4D", lw=1.0, label=f"Wild-type {res['n_wt']}")
    ax.set_title(title, fontsize=7)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(0, 180)
    ax.yaxis.grid(True, color="#e4e3df", lw=0.4)
    ax.set_axisbelow(True)
    ax.text(0.98, 0.08, f"log-rank p = {res['p']:.3g}\n{res['stat']}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.5, color="#1a1a1a")
    ax.legend(frameon=False, fontsize=5.5, loc="lower left")
for ax in axes[1, :]:
    ax.set_xlabel("Months")
axes[0, 0].set_ylabel("Overall survival")
axes[1, 0].set_ylabel("Overall survival")
fig.tight_layout()
fig.savefig(OUT / "figures" / "FigS_catalogue_survival.pdf")
fig.savefig(OUT / "figures" / "FigS_catalogue_survival.png", dpi=200)

# Clinical figure: fraction genome altered.
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.8), sharey=True)
clin_draw = [
    ("PIK3CA", "any mutation", "PIK3CA"),
    ("TP53", "any mutation", "TP53"),
    ("KMT2C", "truncating or splice", "KMT2C truncating"),
]
for ax, (gene, contrast, title) in zip(axes, clin_draw):
    mut_ids = dict(((g, c), ids) for g, c, ids in contrasts)[(gene, contrast)]
    mut, wt = split(mut_ids)
    a = numeric_pair(mut, lambda p: to_float(sample_clin[primary[p]].get("FRACTION_GENOME_ALTERED")))
    b = numeric_pair(wt, lambda p: to_float(sample_clin[primary[p]].get("FRACTION_GENOME_ALTERED")))
    bp = ax.boxplot([a, b], widths=0.55, patch_artist=True, showfliers=False,
                    medianprops={"color": "#1a1a1a", "lw": 0.8})
    for patch, color in zip(bp["boxes"], ("#2a78d6", "#d9d9d9")):
        patch.set_facecolor(color)
        patch.set_edgecolor("#1a1a1a")
        patch.set_linewidth(0.4)
    res = lookup[(gene, contrast, "fraction_genome_altered")]
    ax.set_xticklabels(["Mutated", "Wild-type"])
    ax.set_title(f"{title}\nq = {dict(zip([(g,c,v) for g,c,v,_ in tests], q))[(gene, contrast, 'fraction_genome_altered')]:.3g}", fontsize=7)
    ax.yaxis.grid(True, color="#e4e3df", lw=0.4)
    ax.set_axisbelow(True)
axes[0].set_ylabel("Fraction of genome altered")
fig.tight_layout()
fig.savefig(OUT / "figures" / "FigS_catalogue_clinical.pdf")
fig.savefig(OUT / "figures" / "FigS_catalogue_clinical.png", dpi=200)
print("figures written")
