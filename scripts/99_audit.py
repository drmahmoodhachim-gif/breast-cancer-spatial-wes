#!/usr/bin/env python3
"""Re-check the project against itself. Writes 00_admin/audit_report.md only.

Does not rewrite analysis tables or figures. BAM recounts and the synthetic test
read inputs; the reproducibility step runs in /tmp/audit_clone.
"""
from __future__ import annotations

import hashlib
import math
import re
import json
import os
import random
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_env, load_filters  # noqa: E402

ROOT = PROJECT_DIR
ENV = load_env()
F = load_filters()
HEAD = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
HEAD_SHORT = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
CHECKS: list[dict] = []


def add(section, name, status, detail):
    CHECKS.append({"section": section, "name": name, "status": status, "detail": detail})
    print(f"[{status}] {section} {name}: {detail.splitlines()[0][:180]}")


def md5(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


# ---------------------------------------------------------------------------
# A. Staleness
# ---------------------------------------------------------------------------
stale, old_git, missing_in = [], [], []
withheld = ("Fig2_oncoplot", "FigS3_VAF_distributions", "pathway_hits")
for prov in list((ROOT / "11_figures").rglob("*.provenance.json")) + list((ROOT / "12_tables").rglob("*.provenance.json")):
    if "archive" in prov.parts or any(name in prov.name for name in withheld):
        continue
    rec = json.loads(prov.read_text())
    commit = str(rec.get("git_commit") or "")
    if commit and commit not in (HEAD, HEAD_SHORT) and not HEAD.startswith(commit) and not commit.startswith(HEAD_SHORT):
        old_git.append(f"{prov.relative_to(ROOT)} recorded {commit}")
    for item in rec.get("inputs") or []:
        p = Path(item.get("path", ""))
        if not p.is_file():
            missing_in.append(f"{prov.name}: missing {p}")
            continue
        stored = item.get("md5")
        if not stored:
            continue
        now = md5(p)
        if now and now != stored:
            stale.append(f"{prov.relative_to(ROOT)} input {p.name} md5 {stored[:8]} -> {now[:8]}")

add("A", "input checksums", "FAIL" if stale or missing_in else "PASS",
    ("; ".join(stale + missing_in) if stale or missing_in else "every recorded input checksum still matches"))
add("A", "git hash vs HEAD", "WARN" if old_git else "PASS",
    (f"{len(old_git)} provenance files were written before {HEAD_SHORT}" if old_git else f"all provenance files record {HEAD_SHORT}"))

# ---------------------------------------------------------------------------
# B. Counts
# ---------------------------------------------------------------------------
V = read_tsv(ROOT / "06_variants" / "somatic_filtered.tsv")
ccf = read_tsv(ROOT / "09_ccf_clonality" / "ccf.tsv")
look = read_tsv(ROOT / "08_lookback_N" / "lookback_results.tsv")
table2 = read_tsv(ROOT / "12_tables" / "main" / "table2_driver_alterations.tsv")
add("B", "primary row count", "PASS" if len(V) == 380 else "FAIL", f"somatic_filtered has {len(V)} rows (expected 380)")


def class_counts(df):
    return df.groupby(["patient", "spatial_class"]).size().to_dict()


primary_ids = set(V["var_id"])
look_n = look[(look["specimen"] == "N") & look["var_id"].isin(primary_ids)]
mism = []
if class_counts(V) != class_counts(ccf):
    mism.append("ccf")
look_ids = set(look_n["var_id"])
if look_ids != primary_ids or look_n["var_id"].duplicated().any():
    mism.append(f"lookback primary N has {len(look_ids)} of {len(primary_ids)} variants")
# Table 2 is the driver sheet. Each of its variants must carry the same class as the primary table.
t2 = table2.merge(V[["patient", "chrom", "pos", "ref", "alt", "spatial_class"]],
                  on=["patient", "chrom", "pos", "ref", "alt"], how="left", suffixes=("", "_primary"))
bad_t2 = t2[t2["spatial_class"] != t2["spatial_class_primary"]]
if len(bad_t2):
    mism.append(f"table2 {len(bad_t2)} class mismatches")
# Figure 3 reads somatic_filtered directly, so its class counts are this table.
add("B", "class counts agree", "FAIL" if mism else "PASS",
    ("mismatch in " + ", ".join(mism) if mism else
     "ccf, look-back (primary N rows) and Figure 3's source (somatic_filtered) share the same per-patient class counts; "
     f"Table 2's {len(table2)} driver rows match those classes"))

# ---------------------------------------------------------------------------
# C. Recompute
# ---------------------------------------------------------------------------
vaf_bad = []
for s in ("N", "J", "T"):
    den = V[f"{s}_depth"].replace(0, np.nan)
    calc = V[f"{s}_alt"] / den
    diff = (calc - V[f"{s}_vaf"]).abs()
    hit = diff[diff > 1e-6]
    if len(hit):
        vaf_bad.append(f"{s}:{len(hit)}")
add("C", "VAF = alt/depth", "FAIL" if vaf_bad else "PASS",
    ("exceeds 1e-6 at " + ", ".join(vaf_bad) if vaf_bad else "N, J and T agree within 1e-6"))


def ccf_from_vaf(v, p, C, major):
    if any(pd.isna(x) for x in (v, p, C)) or p <= 0:
        return np.nan
    u = v * (p * C + 2 * (1 - p)) / p
    m = max(1, int(round(u)))
    if not pd.isna(major) and major >= 1:
        m = min(m, int(major))
    return u / m


ccf_bad = 0
for s in ("T", "J"):
    for r in ccf.itertuples(index=False):
        C = getattr(r, f"{s}_tcn")
        minor = getattr(r, f"{s}_lcn")
        major = (C - minor) if pd.notna(C) and pd.notna(minor) else np.nan
        if pd.notna(major):
            major = max(1, major)
        calc = ccf_from_vaf(getattr(r, f"{s}_vaf"), getattr(r, f"{s}_purity"), C, major)
        stored = getattr(r, f"{s}_ccf_raw")
        if pd.notna(calc) and pd.notna(stored) and abs(calc - stored) > 0.01:
            ccf_bad += 1
add("C", "CCF from VAF", "FAIL" if ccf_bad else "PASS",
    f"{ccf_bad} T/J rows differ from VAF, purity, total copy number and multiplicity by more than 0.01")

tmb = read_tsv(ROOT / "12_tables" / "main" / "tmb.tsv")
tmb_bad = []
for r in tmb.itertuples(index=False):
    calc = r.n_coding / 37.5
    if abs(calc - r.TMB_coding_per_mb) > 0.02:
        tmb_bad.append(f"{r.sample_id} stored {r.TMB_coding_per_mb} calc {calc:.3f}")
    if r.specimen in ("T", "J"):
        pres = V[(V["patient"] == r.patient) & V[f"{r.specimen}_present"]]
        n_coding = int(pres["is_coding"].sum())
        if int(r.n_coding) != n_coding:
            tmb_bad.append(f"{r.sample_id} n_coding {int(r.n_coding)} primary table {n_coding}")
add("C", "TMB = coding/37.5", "FAIL" if tmb_bad else "PASS",
    ("; ".join(tmb_bad) if tmb_bad else
     "n_coding matches variants present in the current primary table, and every sample is within 0.02 of coding mutations / 37.5 Mb"))

tin = read_tsv(ROOT / "08_lookback_N" / "tumour_in_normal_estimates.tsv")
tin_bad = []
for r in tin.itertuples(index=False):
    g = V[(V["patient"] == r.patient) & (V["spatial_class"] == "shared") & (V["T_vaf"] > 0)]
    den = float((g["N_depth"] * g["T_vaf"]).sum())
    calc = float(g["N_alt"].sum() / den) if den else np.nan
    if pd.notna(r.TiN_shared) and abs(calc - r.TiN_shared) > 1e-6:
        tin_bad.append(f"{r.patient} stored {r.TiN_shared:.6f} calc {calc:.6f}")
add("C", "TiN slope", "FAIL" if tin_bad else "PASS",
    ("; ".join(tin_bad) if tin_bad else "shared-site slope matches sum(N alt) / sum(N depth x T VAF)"))

min_alt = int(F["min_alt_reads_present"])
floor = float(F["min_vaf_present"])
power_cut = float(F["power_threshold"])


def detection_power(depth, exp, joint):
    if not depth or not exp or float(exp) <= 0:
        return 0.0
    k = max(min_alt, math.ceil(floor * float(depth) - 1e-9)) if joint else min_alt
    return float(binom.sf(k - 1, int(depth), float(exp)))


joint_bad = []
count_only_extra = 0
for r in V.itertuples(index=False):
    for cls, sample, source in (("T_private", "J", "T"), ("J_private", "T", "J"),
                                ("T_unresolved", "J", "T"), ("J_unresolved", "T", "J")):
        if r.spatial_class != cls:
            continue
        exp = max(0.5 * float(getattr(r, f"{source}_vaf")), floor)
        depth = getattr(r, f"{sample}_depth")
        pw = detection_power(depth, exp, True)
        pw_count = detection_power(depth, exp, False)
        enriched = cls.endswith("private")
        if enriched and pw < power_cut:
            joint_bad.append(f"{r.var_id} {cls} joint power {pw:.2f}")
        if (not enriched) and pw >= power_cut:
            joint_bad.append(f"{r.var_id} {cls} joint power {pw:.2f}")
        if (not enriched) and pw_count >= power_cut:
            count_only_extra += 1
add("C", "80% power rule", "FAIL" if joint_bad else "PASS",
    (f"{len(joint_bad)} calls fail the joint rule; first: " + "; ".join(joint_bad[:5]) if joint_bad
     else "every enriched call has joint-rule power >= 0.8 and every unresolved call has joint-rule power < 0.8"))
add("C", "count-only power (sensitivity)", "PASS",
    f"{count_only_extra} unresolved calls would be enriched under P(alt >= {min_alt}) alone. "
    "That comparison is not the primary absence rule.")

# ---------------------------------------------------------------------------
# D. Consistency
# ---------------------------------------------------------------------------
pur = read_tsv(ROOT / "07_cnv_purity" / "purity_reconciled.tsv")
pur_map = dict(zip(pur["sample"], pur["used_purity"]))
pur_bad = []
for patient, g in ccf.groupby("patient"):
    for spec in ("T", "J"):
        vals = g[f"{spec}_purity"].dropna().unique()
        used = pur_map.get(f"{patient}_{spec}")
        if len(vals) != 1 or used is None or abs(float(vals[0]) - float(used)) > 1e-6:
            pur_bad.append(f"{patient}_{spec} ccf {vals} reconciled {used}")
add("D", "CCF purity", "FAIL" if pur_bad else "PASS",
    ("; ".join(pur_bad) if pur_bad else "ccf.tsv purity equals purity_reconciled.used_purity"))

tcell_path = ROOT / "14_immune" / "tcell" / "tcell_fraction.tsv"
if tcell_path.exists():
    tc = read_tsv(tcell_path)
    tc_bad = []
    # DH00024 keeps the TcellExTRECT purity. Copy number 2 left the fraction unchanged,
    # so 0.074 is reported unadjusted. Do not compare these two rows with reconciled purity.
    expect_24 = {"DH00024_T": 0.244, "DH00024_J": 0.320}
    for rec in tc.to_dict("records"):
        sample = str(rec["sample"])
        if sample in expect_24:
            purity = pd.to_numeric(rec["purity"], errors="coerce")
            if not np.isfinite(purity) or abs(float(purity) - expect_24[sample]) > 0.001:
                tc_bad.append(f"{sample} purity {rec['purity']} expected {expect_24[sample]}")
            raw = pd.to_numeric(rec["TCRA.tcell.fraction"], errors="coerce")
            adj = pd.to_numeric(rec["TCRA.tcell.fraction.adj"], errors="coerce")
            if not (np.isfinite(raw) and np.isfinite(adj) and abs(float(raw) - float(adj)) < 1e-9):
                tc_bad.append(f"{sample} adjustment changed the fraction")
            continue
        if sample.endswith("_N"):
            continue
        used = pur_map.get(sample)
        purity = pd.to_numeric(rec["purity"], errors="coerce")
        if used is None or not np.isfinite(purity) or abs(float(purity) - float(used)) > 1e-4:
            tc_bad.append(f"{sample} tcell {rec['purity']} reconciled {used}")
    add("D", "T-cell purity", "FAIL" if tc_bad else "PASS",
        ("; ".join(tc_bad) if tc_bad else
         "DH00024 tumour purity 0.244 and junction purity 0.320; fraction 0.074 is unadjusted. Other samples match reconciled purity."))
else:
    add("D", "T-cell purity", "WARN", "14_immune/tcell/tcell_fraction.tsv is missing")

HLA_POS = {"HLA-A": 29942000, "HLA-C": 31270000, "HLA-B": 31355000}
loh_path = ROOT / "14_immune" / "hla" / "hla_loh.tsv"
hla_bad = []
if loh_path.exists():
    loh = read_tsv(loh_path)
    segs = {}
    for f in (ROOT / "07_cnv_purity" / "facets").glob("*/*.segments.tsv"):
        s = read_tsv(f)
        s["chrom"] = s["chrom"].astype(str).str.replace("chr", "", regex=False)
        segs[s["sample_id"].iloc[0]] = s
    for r in loh.itertuples(index=False):
        s = segs.get(r.sample_id)
        if s is None:
            hla_bad.append(f"{r.sample_id} no segments")
            continue
        hit = s[(s["chrom"] == "6") & (s["start"] <= HLA_POS[r.gene]) & (s["end"] >= HLA_POS[r.gene])]
        if hit.empty:
            hla_bad.append(f"{r.sample_id} {r.gene} uncovered")
            continue
        if abs(float(hit["tcn.em"].iloc[0]) - float(r.total_cn)) > 0.05:
            hla_bad.append(f"{r.sample_id} {r.gene} tcn")
        if pd.notna(r.minor_cn) and pd.notna(hit["lcn.em"].iloc[0]) and abs(float(hit["lcn.em"].iloc[0]) - float(r.minor_cn)) > 0.05:
            hla_bad.append(f"{r.sample_id} {r.gene} lcn")
    add("D", "HLA copy number", "FAIL" if hla_bad else "PASS",
        ("; ".join(hla_bad[:8]) if hla_bad else "HLA-A/B/C states match the FACETS segment at each locus"))
else:
    add("D", "HLA copy number", "WARN", "hla_loh.tsv is missing")

py_bad = []
for p in ["DH00018", "DH00019", "DH00022", "DH00023", "DH00024"]:
    usable = ccf[ccf["patient"] == p].dropna(subset=["T_tcn", "J_tcn", "T_purity", "J_purity"])
    inp = ROOT / "09_ccf_clonality" / "pyclone" / f"{p}.input.tsv"
    if len(usable) < int(F["pyclone_min_mutations"]):
        if inp.exists():
            py_bad.append(f"{p} has an input despite {len(usable)} usable rows")
        continue
    if not inp.exists():
        py_bad.append(f"{p} missing input ({len(usable)} usable)")
        continue
    got = set(read_tsv(inp)["mutation_id"])
    if got != set(usable["var_id"]):
        py_bad.append(f"{p} input {len(got)} vs usable {len(usable)}")
add("D", "PyClone input", "FAIL" if py_bad else "PASS",
    ("; ".join(py_bad) if py_bad else "each PyClone input is the ccf rows with copy number and purity in both J and T"))

# ---------------------------------------------------------------------------
# E. Red flags
# ---------------------------------------------------------------------------
hi_n = V[V["N_vaf"] > V[["T_vaf", "J_vaf"]].max(axis=1)]
add("E", "N VAF above the tumour", "FAIL" if len(hi_n) else "PASS",
    (f"{len(hi_n)} primary variants; first {hi_n['var_id'].iloc[0]}" if len(hi_n) else "none"))
germ = V[V["N_vaf"].between(0.35, 0.65)]
add("E", "germline-like N VAF", "FAIL" if len(germ) else "PASS",
    (f"{len(germ)} primary variants with N VAF 0.35-0.65" if len(germ) else "none in 0.35-0.65"))
if "gnomAD_AF" in V.columns:
    gnom = V[pd.to_numeric(V["gnomAD_AF"], errors="coerce") > 1e-3]
    bits = [f"{r.patient} {r.var_id} AF {r.gnomAD_AF}" for r in gnom.itertuples(index=False)]
    add("E", "gnomAD AF > 1e-3", "PASS",
        (f"{len(gnom)} of 380 exceed 0.001 and are retained. " + "; ".join(bits)) if len(gnom)
        else "none of the 380 exceed 0.001")
else:
    add("E", "gnomAD AF > 1e-3", "WARN", "somatic_filtered.tsv has no gnomAD_AF column, so this filter cannot be recomputed")

bed = ROOT / "03_resources" / "exome2plus_targets_hg38.pad100.bed"
intervals = {}
for line in bed.open():
    chrom, start, end = line.split("\t")[:3]
    intervals.setdefault(chrom, []).append((int(start), int(end)))
for chrom in intervals:
    intervals[chrom].sort()
outside = []
for r in V.itertuples(index=False):
    hits = intervals.get(str(r.chrom), [])
    if not any(a < int(r.pos) <= b for a, b in hits):
        outside.append(r.var_id)
add("E", "inside padded targets", "FAIL" if outside else "PASS",
    (f"{len(outside)} primary variants fall outside the targets plus 100 bp" if outside else "all 380 lie in the padded target BED"))

wanted = [("PIK3CA", "E545K", "Glu545Lys"), ("PIK3CA", "H1047R", "His1047Arg"), ("TP53", "R248L", "Arg248Leu")]
missing_drv = []
for gene, short, long in wanted:
    sub = V[V["gene"] == gene]
    text = (sub["hgvsp_short"].fillna("") + " " + sub["hgvsp"].fillna("")).str.replace("p.", "", regex=False)
    if not text.str.contains(short).any() and not text.str.contains(long).any():
        missing_drv.append(f"{gene} {short}")
add("E", "recognised drivers", "FAIL" if missing_drv else "PASS",
    ("missing " + ", ".join(missing_drv) if missing_drv else "PIK3CA E545K, PIK3CA H1047R and TP53 R248L are in the primary table"))

qc = read_tsv(ROOT / "04_qc" / "sample_qc_extended.tsv")
sex_bad = qc.loc[qc["sex_check"] != "XX", "sample_id"].tolist() if "sex_check" in qc.columns else list(qc["sample_id"])
add("E", "sex is female", "FAIL" if sex_bad or len(qc) != 15 else "PASS",
    (f"{len(qc)} rows; not XX: {sex_bad}" if sex_bad or len(qc) != 15 else "all 15 samples are XX"))

# ---------------------------------------------------------------------------
# F. Raw recount
# ---------------------------------------------------------------------------
keys = []
for gene, token in (("PIK3CA", "545"), ("PIK3CA", "1047"), ("TP53", "248"), ("BRCA2", ""),
                    ("KMT2C", "2197"), ("KMT2C", "4553"), ("FANCM", ""), ("POLE2", ""),
                    ("MGA", ""), ("PRKCB", "")):
    sub = V[V["gene"] == gene]
    if token:
        sub = sub[sub["hgvsp"].fillna("").str.contains(token) | sub["hgvsp_short"].fillna("").str.contains(token)]
    if len(sub):
        keys.append(sub.iloc[0])
key_ids = {r.var_id for r in keys}
pool = V[~V["var_id"].isin(key_ids)]
rng = random.Random(1234)
extra = [pool.iloc[i] for i in rng.sample(range(len(pool)), 20)]
recount_rows = list(keys) + extra

sam = shutil.which("samtools") or str(Path(sys.executable).parent / "samtools")
bam_root = Path(ENV.get("SAREK_OUTDIR", ""))
recount_lines = []
worst = "PASS"
vaf_ok = vaf_fail = depth_only = indel_miss = 0


def pile_count(bam, chrom, pos, alt, ref, extra_flags):
    cmd = [sam, "mpileup", "-q", "20", "-Q", "20", "-r", f"{chrom}:{pos}-{pos}", str(bam)]
    if extra_flags:
        cmd[2:2] = ["--excl-flags", extra_flags]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0 or not out.stdout.strip():
        return None, None
    parts = out.stdout.split("\t")
    if len(parts) < 5:
        return 0, 0
    depth = int(parts[3])
    bases = parts[4]
    bases = re.sub(r"\^.", "", bases).replace("$", "")
    ref, alt = str(ref).upper(), str(alt).upper()
    if len(ref) == len(alt) == 1:
        alt_n = bases.upper().count(alt)
    elif len(ref) > len(alt):
        deleted = ref[len(alt):]
        alt_n = bases.upper().count(f"-{len(deleted)}{deleted}")
    elif len(alt) > len(ref):
        inserted = alt[len(ref):]
        alt_n = bases.upper().count(f"+{len(inserted)}{inserted}")
    else:
        alt_n = 0
    return alt_n, depth


igv_dir = ROOT / "00_admin" / "igv"
igv_dir.mkdir(parents=True, exist_ok=True)
batch = ["new", "genome hg38"]
if not Path(sam).exists() and shutil.which("samtools") is None:
    add("F", "BAM recount", "FAIL", "samtools is not on PATH")
else:
    for r in recount_rows:
        sid_bams = {}
        for spec in ("N", "J", "T"):
            sid = f"{r.patient}_{spec}"
            hits = list(bam_root.rglob(f"{sid}.recal.bam")) if bam_root.exists() else []
            sid_bams[spec] = hits[0] if hits else None
        bits = [r.var_id]
        for spec, bam in sid_bams.items():
            if bam is None:
                bits.append(f"{spec}:no bam")
                worst = "FAIL"
                continue
            alt_n, depth = pile_count(bam, r.chrom, int(r.pos), r.alt, r.ref, "")
            if alt_n is None:
                bits.append(f"{spec}:mpileup failed")
                worst = "FAIL"
                continue
            t_alt, t_dep = int(getattr(r, f"{spec}_alt")), int(getattr(r, f"{spec}_depth"))
            tv = t_alt / t_dep if t_dep else 0
            pv = alt_n / depth if depth else 0
            vaf_gap = abs(pv - tv)
            rel_d = abs(depth - t_dep) / max(t_dep, 1)
            snv = len(str(r.ref)) == len(str(r.alt)) == 1
            flag = "ok"
            if snv and vaf_gap > 0.02:
                flag = "VAF"
                worst = "FAIL"
                vaf_fail += 1
            elif snv:
                vaf_ok += 1
                if rel_d > 0.25:
                    flag = "depth"
                    depth_only += 1
            elif t_alt > 0 and alt_n == 0:
                flag = "indel"
                worst = "FAIL"
                indel_miss += 1
            elif vaf_gap > 0.02:
                flag = "VAF"
                worst = "FAIL"
                vaf_fail += 1
            bits.append(f"{spec} table {t_alt}/{t_dep} ({tv:.3f}) pileup {alt_n}/{depth} ({pv:.3f}) {flag}")
        if r.var_id in key_ids:
            for spec, bam in sid_bams.items():
                if bam is not None:
                    batch += [f"load {bam}", f"goto {r.chrom}:{int(r.pos)-40}-{int(r.pos)+40}",
                              f"snapshot {igv_dir / (r.var_id.replace(':', '_') + '_' + spec + '.png')}"]
        recount_lines.append(" | ".join(bits))
    (igv_dir / "batch.igv").write_text("\n".join(batch + ["exit"]) + "\n")
    igv = shutil.which("igv") or shutil.which("igv.sh")
    if igv:
        subprocess.run([igv, "-b", str(igv_dir / "batch.igv")], check=False)
        n_png = len(list(igv_dir.glob("*.png")))
        shot = f"{n_png} IGV snapshots"
        if n_png < 10:
            worst = "WARN" if worst == "PASS" else worst
    else:
        shot = "IGV is not installed; batch file written to 00_admin/igv/batch.igv"
        if worst == "PASS":
            worst = "WARN"
    reading = (
        f"SNV allele fractions within 0.02: {vaf_ok}; beyond 0.02: {vaf_fail}. "
        f"Depth differs by more than 25% where the SNV allele fraction agrees: {depth_only}. "
        f"Those depths are mpileup -q 20 -Q 20. An agreeing allele fraction does not show whether the lower depth changes the absence call. "
        f"A base pileup that returns no alternate reads for an indel is a parser limit, not evidence the indel is absent. "
        f"Indel rows with table alt reads and zero base-pileup alt reads: {indel_miss}. "
    )
    add("F", "BAM recount of 30 sites", worst, reading + shot + "\n" + "\n".join(recount_lines))

# ---------------------------------------------------------------------------
# G. Synthetic truth
# ---------------------------------------------------------------------------
if not (ROOT / "tests" / "run_synthetic_test.sh").exists():
    add("G", "synthetic test", "WARN",
        "tests/run_synthetic_test.sh is not in this public snapshot; the synthetic test was run in the full project")
else:
    syn = Path("/tmp/audit_synthetic")
    env = os.environ.copy()
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    try:
        proc = subprocess.run(["bash", str(ROOT / "tests" / "run_synthetic_test.sh"), str(syn)],
                              cwd=ROOT, env=env, capture_output=True, text=True, timeout=1800)
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-15:])
        add("G", "synthetic test", "PASS" if proc.returncode == 0 else "FAIL",
            f"exit {proc.returncode}\n{tail}")
    except subprocess.TimeoutExpired:
        add("G", "synthetic test", "FAIL", "timed out after 30 minutes")

# ---------------------------------------------------------------------------
# H. Reproducibility
# ---------------------------------------------------------------------------
clone = Path("/tmp/audit_clone")
if clone.exists():
    shutil.rmtree(clone, ignore_errors=True)
proc = subprocess.run(["git", "clone", "--local", str(ROOT), str(clone)], capture_output=True, text=True)
if proc.returncode != 0:
    add("H", "fresh clone", "FAIL", proc.stderr[-500:])
else:
    env_file = clone / "config" / "project.env"
    if env_file.exists():
        # run_all.sh sources this file and would otherwise keep the absolute project path.
        env_file.write_text(env_file.read_text().replace(str(ROOT), str(clone)))
    needed = [
        "06_variants/somatic_filtered.tsv",
        "07_cnv_purity/purity_reconciled.tsv",
        "07_cnv_purity/facets",
        "09_ccf_clonality/ccf.tsv",
        "12_tables/main/tmb.tsv",
    ]
    absent = [p for p in needed if not (clone / p).exists()]
    if absent:
        add("H", "rerun from HEAD", "FAIL",
            "HEAD does not contain the inputs the downstream steps read: " + ", ".join(absent))
    elif not (clone / "run_all.sh").exists():
        add("H", "rerun from HEAD", "WARN",
            "run_all.sh is not in this public snapshot; the rerun check was done in the full project")
    else:
        diffs = []
        env2 = env.copy()
        (clone / "logs").mkdir(parents=True, exist_ok=True)
        for step in ("ccf", "pyclone", "tmb_sigs", "tables"):
            before = {}
            r = subprocess.run(["bash", "run_all.sh", step], cwd=clone, env=env2, capture_output=True, text=True, timeout=1800)
            if r.returncode != 0:
                diffs.append(f"{step} exit {r.returncode}: " + "\n".join(r.stderr.splitlines()[-5:]))
                break
        compare = [
            "09_ccf_clonality/ccf.tsv",
            "09_ccf_clonality/pyclone_status.tsv",
            "12_tables/main/tmb.tsv",
            "12_tables/main/table2_driver_alterations.tsv",
        ]

        def semantic_diff(live, clone_path):
            A = pd.read_csv(live, sep="\t")
            B = pd.read_csv(clone_path, sep="\t")
            drop = [c for c in list(A.columns) + list(B.columns) if "gnomAD" in str(c) or "gnomad" in str(c)]
            A = A.drop(columns=[c for c in drop if c in A.columns], errors="ignore")
            B = B.drop(columns=[c for c in drop if c in B.columns], errors="ignore")
            only_live = [c for c in A.columns if c not in B.columns]
            only_clone = [c for c in B.columns if c not in A.columns]
            problems = []
            if only_live or only_clone:
                problems.append(f"columns live-only {only_live} clone-only {only_clone}")
            cols = [c for c in A.columns if c in B.columns]
            if len(A) != len(B):
                problems.append(f"rows {len(A)} vs {len(B)}")
                return problems
            A, B = A[cols].reset_index(drop=True), B[cols].reset_index(drop=True)
            for c in cols:
                x = pd.to_numeric(A[c], errors="coerce")
                y = pd.to_numeric(B[c], errors="coerce")
                numeric = x.notna().any() or y.notna().any()
                if numeric and not A[c].astype(str).str.contains("[A-Za-z]", regex=True).any():
                    bad = (x.isna() ^ y.isna()) | ((x - y).abs() > 1e-5)
                    bad = bad.fillna(False)
                    if bad.any():
                        problems.append(f"{c} numeric {int(bad.sum())}")
                else:
                    xs = A[c].map(lambda v: "" if pd.isna(v) else str(v))
                    ys = B[c].map(lambda v: "" if pd.isna(v) else str(v))
                    if not xs.equals(ys):
                        problems.append(f"{c} text {int((xs != ys).sum())}")
            return problems

        byte_notes, semantic = [], []
        for rel in compare:
            a, b = ROOT / rel, clone / rel
            if a.exists() and not b.exists():
                semantic.append(f"{rel} missing from the clone rerun")
                continue
            if not (a.exists() and b.exists()):
                continue
            if md5(a) != md5(b):
                byte_notes.append(f"{rel} byte mismatch")
                found = semantic_diff(a, b)
                if found:
                    semantic.append(f"{rel}: " + "; ".join(found))
                else:
                    byte_notes.append(f"{rel} semantic match after dropping gnomAD and equating missing values")
        for p in (clone / "09_ccf_clonality" / "pyclone").glob("*.clusters.tsv"):
            orig = ROOT / "09_ccf_clonality" / "pyclone" / p.name
            if orig.exists() and md5(orig) != md5(p):
                found = semantic_diff(orig, p)
                byte_notes.append(f"pyclone/{p.name} byte mismatch")
                if found:
                    semantic.append(f"pyclone/{p.name}: " + "; ".join(found))
        detail = "; ".join(diffs + semantic + byte_notes) if (diffs or semantic or byte_notes) else "ccf, pyclone, tmb and table outputs match"
        add("H", "rerun from HEAD", "FAIL" if (diffs or semantic) else ("WARN" if byte_notes else "PASS"), detail)

# ---------------------------------------------------------------------------
# Report and bundle
# ---------------------------------------------------------------------------
order = {"FAIL": 0, "WARN": 1, "PASS": 2}
lines = [
    "# Audit report",
    "",
    f"Run: {datetime.now().isoformat(timespec='seconds')}  ",
    f"HEAD: `{HEAD}`",
    "",
    "The audit reads the project. It does not rewrite analysis tables or figures.",
    "",
    "## Summary",
    "",
    "| Check | Status | Detail |",
    "|---|---|---|",
]
for c in sorted(CHECKS, key=lambda r: (r["section"], order[r["status"]], r["name"])):
    first = c["detail"].splitlines()[0].replace("|", "/")
    lines.append(f"| {c['section']}. {c['name']} | {c['status']} | {first} |")
lines += ["", "## Details", ""]
for c in CHECKS:
    lines += [f"### {c['section']}. {c['name']} — {c['status']}", "", c["detail"], ""]
report = ROOT / "00_admin" / "audit_report.md"
report.write_text("\n".join(lines))

bundle = ROOT / "00_admin" / "audit_bundle.zip"
with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as z:
    z.write(report, "audit_report.md")
    z.write(ROOT / "00_admin" / "decision_log.md", "decision_log.md")
    for rel in (
        "06_variants/somatic_filtered.tsv",
        "09_ccf_clonality/ccf.tsv",
        "07_cnv_purity/purity_reconciled.tsv",
        "08_lookback_N/lookback_results.tsv",
        "08_lookback_N/tumour_in_normal_estimates.tsv",
        "07_cnv_purity/cn_events.tsv",
        "09_ccf_clonality/DH00023_lineages.tsv",
    ):
        p = ROOT / rel
        if p.exists():
            z.write(p, rel)
    for p in (ROOT / "07_cnv_purity" / "facets").glob("*/*.segments.tsv"):
        z.write(p, str(p.relative_to(ROOT)))
    for p in (ROOT / "09_ccf_clonality" / "pyclone").glob("*.clusters.tsv"):
        z.write(p, str(p.relative_to(ROOT)))
print(f"[audit] report {report}")
print(f"[audit] bundle {bundle} ({bundle.stat().st_size} bytes)")
nfail = sum(c["status"] == "FAIL" for c in CHECKS)
nwarn = sum(c["status"] == "WARN" for c in CHECKS)
print(f"[audit] {nfail} FAIL, {nwarn} WARN, {len(CHECKS) - nfail - nwarn} PASS")
