#!/usr/bin/env python3
"""Spike a known fraction of tumour reads into the adjacent normal and re-estimate TiN.

samtools view -s keeps that fraction of templates from the recalibrated tumour BAM.
Those reads are merged into the recalibrated normal BAM. Alleles are recounted at
the patient's PASS sites with the same filters as the look-back. The observed slope
is the read-weighted through-origin TiN on shared sites. The expected slope adds
fraction f of the original tumour counts to the original normal counts:

    (N_alt + f * T_alt) / ((N_depth + f * T_depth) * T_vaf)

The CI is the 2.5–97.5 percentile of 1000 site bootstraps of the observed slope.
"""
import subprocess
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import linregress

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_env, load_filters, write_provenance  # noqa: E402

env = load_env()
F = load_filters()
RES = Path(env["SAREK_OUTDIR"])
SAMTOOLS = env.get("SAMTOOLS") or shutil.which("samtools") or "samtools"
MINMQ, MINBQ = F["lookback_min_mapq"], F["lookback_min_baseq"]
PATIENTS = ["DH00018", "DH00023"]
# 0 keeps no tumour templates (-s 1.000) and still rewrites the normal by the same merge.
FRACTIONS = [0.0, 0.005, 0.01, 0.02, 0.05]
CORE = ["patient", "spiked_fraction", "n_shared", "baseline_slope",
        "expected_slope", "observed_slope", "ci_lo", "ci_hi"]
N_BOOT = 1000
WORK = PROJECT_DIR / "08_lookback_N" / "dilution_bams"
WORK.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(1234)


def bam_of(sid):
    return RES / "preprocessing" / "recalibrated" / sid / f"{sid}.recal.bam"


def count_site(bam, chrom, pos, ref, alt):
    is_snv = len(ref) == 1 and len(alt) == 1
    indel_len = len(alt) - len(ref)
    depth = alt_n = 0
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter",
                          min_base_quality=0, ignore_overlaps=True, max_depth=100000):
        for pr in col.pileups:
            aln = pr.alignment
            if (aln.is_duplicate or aln.is_secondary or aln.is_supplementary or aln.is_qcfail
                    or aln.mapping_quality < MINMQ or pr.is_refskip or pr.is_del):
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            if is_snv:
                bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 40
                if bq < MINBQ:
                    continue
                base = aln.query_sequence[qpos].upper()
                is_alt = base == alt.upper()
                is_ref = base == ref.upper()
            else:
                is_alt = pr.indel == indel_len and indel_len != 0
                is_ref = pr.indel == 0
            if is_alt or is_ref:
                depth += 1
                if is_alt:
                    alt_n += 1
            else:
                depth += 1
    return depth, alt_n


def slope(depth, alt, t_vaf):
    den = np.sum(depth * t_vaf)
    return float(np.sum(alt) / den) if den > 0 else np.nan


def bootstrap(depth, alt, t_vaf):
    n = len(depth)
    if n < 5:
        return np.nan, np.nan
    vals = []
    for _ in range(N_BOOT):
        ix = rng.integers(0, n, size=n)
        vals.append(slope(depth[ix], alt[ix], t_vaf[ix]))
    vals = [v for v in vals if np.isfinite(v)]
    if not vals:
        return np.nan, np.nan
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
out = PROJECT_DIR / "08_lookback_N" / "dilution_test.tsv"
prev = pd.read_csv(out, sep="\t") if out.exists() else pd.DataFrame(columns=CORE)
rows = []
for patient in PATIENTS:
    g = V[(V["patient"] == patient) & (V["spatial_class"] == "shared") & (V["T_vaf"] > 0)].copy()
    g = g.reset_index(drop=True)
    base = slope(g["N_depth"].to_numpy(float), g["N_alt"].to_numpy(float), g["T_vaf"].to_numpy(float))
    n_bam = bam_of(f"{patient}_N")
    t_bam = bam_of(f"{patient}_T")
    done = prev[prev["patient"] == patient] if len(prev) else prev
    for frac in FRACTIONS:
        hit = done[np.isclose(done["spiked_fraction"].astype(float), frac)] if len(done) else done
        if len(hit):
            rec = {c: hit.iloc[0][c] for c in CORE}
            rows.append(rec)
            print(f"[dilution] {patient} {frac:.3f} reused", flush=True)
            continue
        tag = f"{frac * 100:g}".replace(".", "p")
        spike = WORK / f"{patient}_T_frac{tag}.bam"
        merged = WORK / f"{patient}_N_plus_{tag}.bam"
        # -s seed.fraction: 1.005 keeps 0.5%, 1.050 keeps 5%. A zero fraction disables subsampling.
        subsample = f"{1 + frac:.3f}"
        if frac == 0.0:
            print(f"[dilution] {patient} spike 0.000 (header-only tumour BAM, no reads)", flush=True)
        else:
            print(f"[dilution] {patient} spike {frac:.3f} (-s {subsample})", flush=True)
        if frac == 0.0:
            # samtools treats a zero fractional part as "subsampling off" and would keep every read.
            header = spike.with_suffix(".sam")
            subprocess.check_call([SAMTOOLS, "view", "-H", "-o", str(header), str(t_bam)])
            subprocess.check_call([SAMTOOLS, "view", "-b", "-o", str(spike), str(header)])
            header.unlink(missing_ok=True)
            nreads = int(subprocess.check_output([SAMTOOLS, "view", "-c", str(spike)], text=True))
            if nreads != 0:
                raise SystemExit(f"[dilution] 0% spike kept {nreads} reads; refusing to merge")
        else:
            subprocess.check_call([SAMTOOLS, "view", "-s", subsample, "-b", "-@", "2",
                                   "-o", str(spike), str(t_bam)])
        subprocess.check_call([SAMTOOLS, "merge", "-f", "-@", "2", str(merged), str(n_bam), str(spike)])
        subprocess.check_call([SAMTOOLS, "index", str(merged)])
        depth, alt = [], []
        with pysam.AlignmentFile(str(merged)) as bam:
            for r in g.itertuples():
                d, a = count_site(bam, r.chrom, int(r.pos), r.ref, r.alt)
                depth.append(d)
                alt.append(a)
        depth = np.asarray(depth, float)
        alt = np.asarray(alt, float)
        t_vaf = g["T_vaf"].to_numpy(float)
        obs = slope(depth, alt, t_vaf)
        lo, hi = bootstrap(depth, alt, t_vaf)
        exp_depth = g["N_depth"].to_numpy(float) + frac * g["T_depth"].to_numpy(float)
        exp_alt = g["N_alt"].to_numpy(float) + frac * g["T_alt"].to_numpy(float)
        exp = slope(exp_depth, exp_alt, t_vaf)
        rows.append({
            "patient": patient, "spiked_fraction": frac, "n_shared": len(g),
            "baseline_slope": base, "expected_slope": exp, "observed_slope": obs,
            "ci_lo": lo, "ci_hi": hi,
        })
        print(f"[dilution] {patient} {frac:.3f} expected {exp:.4f} observed {obs:.4f} CI {lo:.4f}-{hi:.4f}", flush=True)
        merged.unlink(missing_ok=True)
        Path(str(merged) + ".bai").unlink(missing_ok=True)
        spike.unlink(missing_ok=True)

res = pd.DataFrame(rows)
zero = (res.loc[np.isclose(res["spiked_fraction"].astype(float), 0.0), ["patient", "observed_slope"]]
        .drop_duplicates("patient")
        .rename(columns={"observed_slope": "pipeline_baseline"}))
res = res.merge(zero, on="patient", how="left")
res["observed_minus_baseline"] = res["observed_slope"] - res["pipeline_baseline"]
res["expected_increment"] = res["expected_slope"] - res["baseline_slope"]
res["recovered"] = res["observed_slope"] - res["baseline_slope"]
fit_rows = []
for patient, g in res.groupby("patient"):
    fit = linregress(g["expected_increment"].to_numpy(float), g["recovered"].to_numpy(float))
    fit_rows.append((patient, float(fit.slope), float(fit.intercept), float(fit.rvalue ** 2)))
fit_df = pd.DataFrame(fit_rows, columns=["patient", "fit_slope", "fit_intercept", "fit_r2"])
res = res.merge(fit_df, on="patient", how="left")
cols = ["patient", "spiked_fraction", "n_shared", "baseline_slope", "pipeline_baseline",
        "expected_slope", "observed_slope", "observed_minus_baseline", "expected_increment",
        "recovered", "ci_lo", "ci_hi", "fit_slope", "fit_intercept", "fit_r2"]
res = res[cols].sort_values(["patient", "spiked_fraction"])
res.to_csv(out, sep="\t", index=False)
write_provenance(out, [PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"],
                 {"fractions": FRACTIONS, "patients": PATIENTS, "bootstrap": N_BOOT,
                  "decision": "D047", "pipeline_baseline": "observed slope of the 0% re-merge",
                  "recovered": "observed_slope - baseline_slope",
                  "expected_increment": "expected_slope - baseline_slope"})
print(res.to_string(index=False))
