#!/usr/bin/env python3
"""TMB per sample and a restricted COSMIC signature fit.

TMB = coding somatic variants present in the sample / recorded target size (Mb).
The pad100 BED is byte-identical to the unpadded BED, so this denominator is not enlarged by the 100 bp pad.

Signatures use every PASS SNV present in the sample, including non-coding sites in the
±100 bp flanks, not the coding subset. Each T and J sample is fitted, plus one pooled
T+J catalogue per patient (each variant once). The fit is limited to COSMIC v3.4
SBS1, SBS2, SBS3, SBS4, SBS5 and SBS13. APOBEC is SBS2 + SBS13.
n_snv < 50 is flagged low_count and is still fitted.

Outputs:
  12_tables/main/tmb.tsv
  10_tmb_signatures/input_vcfs/*.vcf
  10_tmb_signatures/signature_exposures.tsv
"""
import inspect
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

SBS_KEEP = ["SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS13"]

V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
V = V[V["filter"].astype(str).eq("PASS")].copy()
size_f = PROJECT_DIR / "03_resources" / "target_size.tsv"
mb = float(pd.read_csv(size_f, sep="\t").set_index("metric").loc["padded_target_mb", "value"]) if size_f.exists() else float("nan")

sig_root = PROJECT_DIR / "10_tmb_signatures"
vdir = sig_root / "input_vcfs"
vdir.mkdir(parents=True, exist_ok=True)


def present_snvs(frame):
    snv = frame[frame["variant_type"] == "SNV"]
    return snv.drop_duplicates(["chrom", "pos", "ref", "alt"])


def write_vcf(path, snv):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in snv.sort_values(["chrom", "pos"]).itertuples():
            fh.write(f"{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\tPASS\t.\n")


catalogues = []
tmb_rows = []
for patient, g in V.groupby("patient"):
    by_sample = {}
    for specimen in ("T", "J"):
        pres = g[g[f"{specimen}_present"]]
        snv = present_snvs(pres)
        by_sample[specimen] = snv
        sample_id = f"{patient}_{specimen}"
        write_vcf(vdir / f"{sample_id}.vcf", snv)
        catalogues.append((sample_id, patient, specimen, snv))
        tmb_rows.append({
            "sample_id": sample_id, "patient": patient, "specimen": specimen,
            "n_variants": int(pres.shape[0]), "n_coding": int(pres["is_coding"].sum()),
            "n_snv": int(snv.shape[0]), "target_mb": mb,
            "TMB_coding_per_mb": round(pres["is_coding"].sum() / mb, 2) if mb == mb else None,
            "signature_fit": "ok" if len(snv) >= 50 else "low_count",
        })
    pooled = present_snvs(pd.concat(by_sample.values(), ignore_index=True))
    sample_id = f"{patient}_pooled"
    write_vcf(vdir / f"{sample_id}.vcf", pooled)
    catalogues.append((sample_id, patient, "pooled", pooled))

tmb = pd.DataFrame(tmb_rows)
tmb_out = PROJECT_DIR / "12_tables" / "main" / "tmb.tsv"
tmb.to_csv(tmb_out, sep="\t", index=False)
write_provenance(tmb_out, [PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"], {"target_mb": mb})
print(tmb.to_string(index=False))
if os.environ.get("TMB_ONLY") == "1":
    sys.exit(0)


def cosmic_sbs_table():
    import SigProfilerAssignment
    root = Path(SigProfilerAssignment.__file__).resolve().parent
    hits = sorted(root.rglob("*SBS*GRCh38*.txt"))
    hits = [p for p in hits if "3.4" in p.name or "v3.4" in p.name or "3.4" in str(p)]
    if not hits:
        hits = sorted(root.rglob("*SBS*GRCh38*.txt"))
    if not hits:
        raise FileNotFoundError("COSMIC SBS GRCh38 table not found inside SigProfilerAssignment")
    # Prefer the genome (not already exome-renormalised) table so exome=True is applied once.
    genome = [p for p in hits if "exome" not in p.name.lower()]
    return genome[0] if genome else hits[0]


def restricted_database(src: Path, keep, out: Path) -> Path:
    db = pd.read_csv(src, sep="\t")
    key = db.columns[0]
    missing = [s for s in keep if s not in db.columns]
    if missing:
        raise SystemExit(f"signatures missing from {src}: {missing}")
    db[[key, *keep]].to_csv(out, sep="\t", index=False)
    return out


try:
    from SigProfilerAssignment import Analyzer as Analyze
except ImportError:
    print("[signatures] SigProfilerAssignment is not installed in the analysis environment")
    sys.exit(1)

src_db = cosmic_sbs_table()
params = inspect.signature(Analyze.cosmic_fit).parameters
SBS_ALT = ["SBS1", "SBS2", "SBS5", "SBS13"]
meta = pd.DataFrame([
    {"sample": sid, "patient": patient, "n_snv": int(snv.shape[0])}
    for sid, patient, _specimen, snv in catalogues
])


def run_fit(keep, db_name, fit_dir: Path):
    db_path = restricted_database(src_db, keep, sig_root / db_name)
    print(f"[signatures] reference {src_db}")
    print(f"[signatures] restricted database {db_path}")
    kwargs = {
        "samples": str(vdir),
        "output": str(fit_dir),
        "input_type": "vcf",
        "context_type": "96",
        "genome_build": "GRCh38",
        "collapse_to_SBS96": True,
        "exome": True,
        "make_plots": False,
        "signature_database": str(db_path),
        "add_background_signatures": False,
    }
    if "cosmic_version" in params:
        kwargs["cosmic_version"] = 3.4
    Analyze.cosmic_fit(**{k: v for k, v in kwargs.items() if k in params or k in ("samples", "output")})
    activity_files = sorted(Path(fit_dir).rglob("*Activities*.txt"))
    if not activity_files:
        raise SystemExit(f"no Activities file under {fit_dir}")
    act = pd.read_csv(activity_files[0], sep="\t")
    sample_col = "Samples" if "Samples" in act.columns else act.columns[0]
    act = act.rename(columns={sample_col: "sample"})
    for sbs in keep:
        if sbs not in act.columns:
            act[sbs] = 0.0
    act["APOBEC"] = act["SBS2"].astype(float) + act["SBS13"].astype(float)
    exp = meta.merge(act[["sample", *keep, "APOBEC"]], on="sample", how="left")
    for sbs in [*keep, "APOBEC"]:
        exp[sbs] = exp[sbs].fillna(0)
    exp["low_count"] = exp["n_snv"] < 50
    return exp, db_path, activity_files[0]


def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom else np.nan


def _cosine_from_stats(fit_dir: Path, sample: str):
    """Cosine similarity of the observed spectrum versus the signature reconstruction."""
    stats = list(Path(fit_dir).rglob("*Samples_Stats.txt"))
    if not stats:
        return np.nan
    tab = pd.read_csv(stats[0], sep="\t")
    name = "Sample Names" if "Sample Names" in tab.columns else tab.columns[0]
    col = "Cosine Similarity" if "Cosine Similarity" in tab.columns else None
    if col is None:
        return np.nan
    hit = tab[tab[name].astype(str) == sample]
    return float(hit[col].iloc[0]) if len(hit) else np.nan


primary, db_path, act_file = run_fit(SBS_KEEP, "COSMIC_v3.4_SBS1_2_3_4_5_13.tsv", sig_root / "output")
primary = primary[["sample", "patient", "n_snv", *SBS_KEEP, "APOBEC", "low_count"]]
exp_out = sig_root / "signature_exposures.tsv"
primary.to_csv(exp_out, sep="\t", index=False)
write_provenance(exp_out, [PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", db_path],
                 {"cosmic": "v3.4", "signatures": ",".join(SBS_KEEP), "activities": str(act_file)})
print(primary.to_string(index=False))
print(f"[signatures] wrote {exp_out}")

alt, alt_db, alt_act = run_fit(SBS_ALT, "COSMIC_v3.4_SBS1_2_5_13.tsv", sig_root / "output_no_sbs3")
# Cosine similarity of the two exposure vectors (SBS3 is 0 in the refit).
vec_cols = SBS_KEEP
left = primary.set_index("sample")
right = alt.set_index("sample")
alt["cosine_similarity"] = [
    cosine([left.loc[s, c] if c in left.columns else 0 for c in vec_cols],
           [right.loc[s, c] if c in right.columns else 0 for c in vec_cols])
    for s in alt["sample"]
]
alt["cosine_observed_no_sbs3"] = alt["sample"].map(
    lambda s: _cosine_from_stats(sig_root / "output_no_sbs3", s))
alt["cosine_observed"] = alt["sample"].map(
    lambda s: _cosine_from_stats(sig_root / "output", s))
alt = alt[["sample", "patient", "n_snv", *SBS_ALT, "APOBEC", "low_count",
           "cosine_similarity", "cosine_observed", "cosine_observed_no_sbs3"]]
alt_out = sig_root / "signature_exposures_alt.tsv"
alt.to_csv(alt_out, sep="\t", index=False)
write_provenance(alt_out, [exp_out, alt_db],
                 {"cosmic": "v3.4", "signatures": ",".join(SBS_ALT),
                  "cosine": "exposure vector vs the SBS1/2/3/4/5/13 fit", "activities": str(alt_act)})
print(alt.to_string(index=False))
print(f"[signatures] wrote {alt_out}")
