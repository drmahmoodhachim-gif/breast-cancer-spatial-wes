#!/usr/bin/env python3
"""PyClone-VI clustering of T and J jointly, per patient — only where mutations suffice.

Patients with fewer than `pyclone_min_mutations` usable variants (with copy number in both
T and J) are skipped and reported as 'CCF only' — this is stated in the Methods/limitations.

Outputs: 09_ccf_clonality/pyclone/<patient>.input.tsv, <patient>.results.tsv,
         <patient>.clusters.tsv, 09_ccf_clonality/pyclone_status.tsv
"""
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PATIENTS, PROJECT_DIR, load_filters, write_provenance  # noqa: E402

F = load_filters()
D = PROJECT_DIR / "09_ccf_clonality"
ccf = pd.read_csv(D / "ccf.tsv", sep="\t")
have_pyclone = shutil.which("pyclone-vi") is not None
status = []
PP = pd.read_csv(PROJECT_DIR / "07_cnv_purity" / "purity_ploidy.tsv", sep="\t").set_index("sample_id")
_segs = pd.concat(pd.read_csv(f, sep="\t") for f in (PROJECT_DIR / "07_cnv_purity" / "facets").glob("*/*.segments.tsv"))
_segs["chrom"] = _segs["chrom"].astype(str).str.replace("chr", "", regex=False)
SEG = {sid: g for sid, g in _segs.groupby("sample_id")}


def local_cn(sid, chrom, pos):
    g = SEG.get(sid)
    if g is None:
        return np.nan
    hit = g[(g["chrom"] == str(chrom).replace("chr", "")) & (g["start"] <= pos) & (g["end"] >= pos)]
    if hit.empty:
        return np.nan
    return float(hit["tcn.em"].iloc[0])


def sample_stats(frame, sample_id):
    s = frame[frame["sample_id"] == sample_id].drop_duplicates("mutation_id")
    if s.empty:
        return np.nan, np.nan
    if s["cluster_id"].nunique() == 1:
        return float(s["cellular_prevalence"].iloc[0]), float(s["cellular_prevalence_std"].iloc[0])
    return float(s["cellular_prevalence"].mean()), float(s["cellular_prevalence"].std(ddof=0))


def frac_gained(patient, mutation_ids):
    """Fraction of mutations on segments with total CN above rounded ploidy (tumour)."""
    sid = f"{patient}_T"
    ploidy = PP["ploidy"].get(sid, np.nan)
    if pd.isna(ploidy):
        return np.nan
    cutoff = int(round(float(ploidy)))
    sub = ccf[ccf["var_id"].isin(mutation_ids)]
    flags = []
    for r in sub.itertuples(index=False):
        tcn = local_cn(sid, r.chrom, r.pos)
        if pd.notna(tcn):
            flags.append(tcn > cutoff)
    return float(np.mean(flags)) if flags else np.nan


def write_clusters(patient, res_path):
    R = pd.read_csv(res_path, sep="\t")
    n_mut = R.groupby("cluster_id")["mutation_id"].nunique()
    small = set(n_mut[n_mut < 5].index)
    rows = []
    for cid in sorted(c for c in n_mut.index if c not in small):
        sub = R[R["cluster_id"] == cid]
        prev_j, sd_j = sample_stats(sub, f"{patient}_J")
        prev_t, sd_t = sample_stats(sub, f"{patient}_T")
        rec = {"cluster": int(cid), "n_mut": int(n_mut[cid]),
               "prev_J": prev_j, "prev_T": prev_t, "sd_J": sd_j, "sd_T": sd_t}
        if patient == "DH00019":
            rec["frac_in_gained_segment"] = frac_gained(patient, set(sub["mutation_id"]))
        for key in ("prev_J", "prev_T", "sd_J", "sd_T", "frac_in_gained_segment"):
            if key in rec and pd.notna(rec[key]):
                rec[key] = round(float(rec[key]), 4)
        rows.append(rec)
    if small:
        sub = R[R["cluster_id"].isin(small)]
        prev_j, sd_j = sample_stats(sub, f"{patient}_J")
        prev_t, sd_t = sample_stats(sub, f"{patient}_T")
        rec = {"cluster": "unassigned", "n_mut": int(sub["mutation_id"].nunique()),
               "prev_J": prev_j, "prev_T": prev_t, "sd_J": sd_j, "sd_T": sd_t}
        if patient == "DH00019":
            rec["frac_in_gained_segment"] = frac_gained(patient, set(sub["mutation_id"]))
        for key in ("prev_J", "prev_T", "sd_J", "sd_T", "frac_in_gained_segment"):
            if key in rec and pd.notna(rec[key]):
                rec[key] = round(float(rec[key]), 4)
        rows.append(rec)
    out = res_path.with_name(f"{patient}.clusters.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    write_provenance(out, [res_path], {"unassigned_below": 5})
    return out


for p in PATIENTS:
    sub = ccf[(ccf["patient"] == p)].copy()
    usable = sub.dropna(subset=["T_tcn", "J_tcn", "T_purity", "J_purity"])
    n = len(usable)
    if n < F["pyclone_min_mutations"]:
        status.append({"patient": p, "usable_mutations": n, "status": "skipped_CCF_only"})
        continue
    rows = []
    for r in usable.itertuples(index=False):
        for s in ("T", "J"):
            tcn = int(round(getattr(r, f"{s}_tcn")))
            lcn = getattr(r, f"{s}_lcn")
            minor = int(lcn) if pd.notna(lcn) else 0
            major = max(1, tcn - minor)
            rows.append({"mutation_id": r.var_id, "sample_id": f"{p}_{s}",
                         "ref_counts": int(getattr(r, f"{s}_depth") - getattr(r, f"{s}_alt")),
                         "alt_counts": int(getattr(r, f"{s}_alt")), "major_cn": major, "minor_cn": minor,
                         "normal_cn": 2, "tumour_content": float(getattr(r, f"{s}_purity"))})
    inp = D / "pyclone" / f"{p}.input.tsv"
    pd.DataFrame(rows).to_csv(inp, sep="\t", index=False)
    if not have_pyclone:
        status.append({"patient": p, "usable_mutations": n, "status": "input_written_pyclone_not_installed"})
        continue
    h5, res = D / "pyclone" / f"{p}.h5", D / "pyclone" / f"{p}.results.tsv"
    subprocess.run(["pyclone-vi", "fit", "-i", str(inp), "-o", str(h5), "-c", str(F["pyclone_num_clusters"]),
                    "-d", "beta-binomial", "-r", str(F["pyclone_restarts"]), "--seed", "1234"], check=True)
    subprocess.run(["pyclone-vi", "write-results-file", "-i", str(h5), "-o", str(res)], check=True)
    write_provenance(res, [inp], {k: F[k] for k in F if k.startswith("pyclone")})
    write_clusters(p, res)
    status.append({"patient": p, "usable_mutations": n, "status": "clustered"})

pd.DataFrame(status).to_csv(D / "pyclone_status.tsv", sep="\t", index=False)
print(pd.DataFrame(status).to_string(index=False))


def sum_rule_rows():
    """Children of the trunk cannot sum to more than the trunk prevalence in a sample."""
    rows = []
    for path in sorted((D / "pyclone").glob("*.clusters.tsv")):
        patient = path.name.split(".")[0]
        cl = pd.read_csv(path, sep="\t")
        cl = cl[cl["cluster"].astype(str) != "unassigned"]
        if cl.empty:
            continue
        for specimen, col in (("J", "prev_J"), ("T", "prev_T")):
            trunk_i = cl[col].idxmax()
            trunk_prev = float(cl.loc[trunk_i, col])
            sum_children = float(cl.drop(index=trunk_i)[col].sum())
            rows.append({
                "patient": patient,
                "sample": f"{patient}_{specimen}",
                "trunk_cluster": cl.loc[trunk_i, "cluster"],
                "trunk_prev": round(trunk_prev, 4),
                "sum_children": round(sum_children, 4),
                "sum_rule_ok": sum_children <= trunk_prev + 0.05,
            })
    return pd.DataFrame(rows)


sum_path = D / "pyclone" / "sum_rule_check.tsv"
sum_df = sum_rule_rows()
sum_df.to_csv(sum_path, sep="\t", index=False)
write_provenance(sum_path, list((D / "pyclone").glob("*.clusters.tsv")),
                 {"rule": "sum of non-trunk cluster prevalences <= trunk prevalence + 0.05; unassigned excluded"})
print(sum_df.to_string(index=False))
