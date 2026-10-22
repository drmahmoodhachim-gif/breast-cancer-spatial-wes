#!/usr/bin/env python3
"""D066. Review of the nine gnomAD AF>0.001 variants. None are dropped."""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance

V = pd.read_csv(PROJECT_DIR / "06_variants" / "somatic_filtered.tsv", sep="\t")
both = pd.read_csv(PROJECT_DIR / "06_variants" / "spatial_class_both_rules.tsv", sep="\t")
high = V[V["gnomAD_AF"].fillna(0) > 0.001].merge(
    both[["var_id", "spatial_class_reads_only"]], on="var_id", how="left"
)
keys = ["chrom", "pos", "ref", "alt"]
n_pat = V.groupby(keys)["patient"].nunique()
rows = []
for r in high.itertuples(index=False):
    chrom = str(r.chrom).replace("chr", "")
    rows.append({
        "patient": r.patient,
        "var_id": r.var_id,
        "gene": "" if pd.isna(r.gene) else r.gene,
        "consequence": r.consequence,
        "gnomAD_AF": r.gnomAD_AF,
        "N_alt": r.N_alt, "N_depth": r.N_depth, "N_vaf": r.N_vaf,
        "J_alt": r.J_alt, "J_depth": r.J_depth, "J_vaf": r.J_vaf,
        "T_alt": r.T_alt, "T_depth": r.T_depth, "T_vaf": r.T_vaf,
        "strelka_T": r.strelka_T, "strelka_J": r.strelka_J,
        "segdup": r.segdup, "mappability": r.mappability, "low_mappability": r.low_mappability,
        "frac_mapq_lt50_N": r.frac_mapq_lt50_N,
        "frac_mapq_lt50_J": r.frac_mapq_lt50_J,
        "frac_mapq_lt50_T": r.frac_mapq_lt50_T,
        "spatial_class_joint": r.spatial_class,
        "spatial_class_reads_only": r.spatial_class_reads_only,
        "n_patients_with_allele": int(n_pat.loc[(r.chrom, r.pos, r.ref, r.alt)]),
        "recurrent_chr9_64407964": chrom == "9" and int(r.pos) == 64407964,
        "kept_in_primary": True,
    })
out = PROJECT_DIR / "12_tables" / "supplementary" / "gnomad_af_review.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
write_provenance(out, inputs=[PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"], params={"decision": "D066"})
print(pd.DataFrame(rows)[["patient", "gene", "gnomAD_AF", "N_alt", "J_alt", "T_alt", "strelka_T", "strelka_J", "segdup", "low_mappability", "n_patients_with_allele"]].to_string(index=False))
