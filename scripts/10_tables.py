#!/usr/bin/env python3
"""Assemble manuscript tables (TSV + one Excel workbook for submission).

Main
  Table 1  per-sample summary: clinical, pathology cellularity, sequencing QC, purity/ploidy, TMB
  Table 2  driver-gene alterations with N/J/T VAF, T/J CCF, spatial class, N look-back class
Supplementary
  S1 FASTQ pre-flight   S2 sequencing QC   S3 all filtered somatic variants
  S4 look-back results  S5 purity/ploidy   S6 thresholds used (filters.yaml)
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, load_filters, write_provenance  # noqa: E402

P = PROJECT_DIR
T = P / "12_tables"


def rd(path):
    path = P / path
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


clin, path_ = rd("config/clinical_metadata.tsv"), rd("config/section_pathology.tsv")
qc, pp, tmb = rd("04_qc/qc_summary.tsv"), rd("07_cnv_purity/purity_ploidy.tsv"), rd("12_tables/main/tmb.tsv")
ccf, lb = rd("09_ccf_clonality/ccf.tsv"), rd("08_lookback_N/lookback_results.tsv")
var = ccf if len(ccf) else rd("06_variants/somatic_filtered.tsv")

# ---- Table 1 ------------------------------------------------------------------
t1 = path_.copy()
for df, key in ((qc.drop(columns=["patient", "specimen"], errors="ignore"), "sample_id"), (pp, "sample_id"),
                (tmb[["sample_id", "n_variants", "TMB_coding_per_mb"]] if len(tmb) else pd.DataFrame(), "sample_id")):
    if len(df):
        t1 = t1.merge(df, on=key, how="left")
if len(clin):
    t1 = t1.merge(clin, on="patient", how="left")
t1.to_csv(T / "main" / "table1_sample_summary.tsv", sep="\t", index=False)

# ---- Table 2 ------------------------------------------------------------------
t2 = pd.DataFrame()
if len(var):
    d = var[var["is_driver_gene"] & var["is_coding"]].copy()
    if len(lb):
        n_lb = lb[lb["specimen"] == "N"][["var_id", "evidence_class", "q_value"]].rename(
            columns={"evidence_class": "N_lookback_class", "q_value": "N_lookback_q"})
        d = d.merge(n_lb, on="var_id", how="left")
    cols = [c for c in ["patient", "gene", "hgvsp", "hgvsc", "consequence", "chrom", "pos", "ref", "alt",
                        "N_vaf", "J_vaf", "T_vaf", "N_alt", "N_depth", "J_ccf", "J_ccf_lo", "J_ccf_hi",
                        "T_ccf", "T_ccf_lo", "T_ccf_hi", "spatial_class", "tier", "N_lookback_class",
                        "N_lookback_q"] if c in d.columns]
    t2 = d[cols].sort_values(["patient", "gene"])
    t2.to_csv(T / "main" / "table2_driver_alterations.tsv", sep="\t", index=False)

supp = {
    "S1_fastq_preflight": rd("02_fastq_check/fastq_headcheck.tsv"),
    "S2_sequencing_qc": qc,
    "S3_somatic_variants": var,
    "S4_lookback_N": lb,
    "S5_purity_ploidy": pp,
    "S6_thresholds": pd.DataFrame([(k, str(v)) for k, v in load_filters().items()], columns=["parameter", "value"]),
}
for name, df in supp.items():
    if len(df):
        df.to_csv(T / "supplementary" / f"{name}.tsv", sep="\t", index=False)

try:
    xlsx = T / "tables_for_submission.xlsx"
    with pd.ExcelWriter(xlsx) as xw:
        t1.to_excel(xw, sheet_name="Table1", index=False)
        if len(t2):
            t2.to_excel(xw, sheet_name="Table2_drivers", index=False)
        for name, df in supp.items():
            if len(df):
                df.to_excel(xw, sheet_name=name[:31], index=False)
    write_provenance(xlsx)
    print(f"[tables] wrote {xlsx}")
except ImportError:
    print("[tables] openpyxl not installed; TSVs written only")
