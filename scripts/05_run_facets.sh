#!/usr/bin/env bash
# FACETS allele-specific copy number, purity and ploidy for T and J (each vs matched N).
# Requires: snp-pileup + R package facets  (conda env: envs/facets.yml)
# Caveat logged in decision_log: if N contains tumour cells, purity in T/J can be underestimated.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$PROJECT_DIR/config/project.env"
OUT="$PROJECT_DIR/07_cnv_purity/facets"
[[ -s "$FACETS_SNP_VCF" ]] || { echo "[ERROR] common SNP VCF missing: $FACETS_SNP_VCF"; exit 1; }

TARGET_SNPS="$PROJECT_DIR/03_resources/dbsnp_common_hg38.targets.vcf.gz"
if [[ ! -s "$TARGET_SNPS" ]]; then
  bcftools view -R "$TARGET_BED_PADDED" -v snps -m2 -M2 "$FACETS_SNP_VCF" -Ou \
    | bcftools sort -Oz -o "$TARGET_SNPS"
  tabix -f -p vcf "$TARGET_SNPS"
  echo "[facets] target SNPs: $(bcftools view -H "$TARGET_SNPS" | wc -l)"
fi
FACETS_SNP_VCF="$TARGET_SNPS"

findbam() {  # final recalibrated BAM (or CRAM) for a sample
  find "$SAREK_OUTDIR/preprocessing" \( -name "$1.recal.bam" -o -name "$1.recal.cram" \) | head -n1
}

for p in DH00018 DH00019 DH00022 DH00023 DH00024; do
  nbam="$(findbam "${p}_N")"
  for s in T J; do
    sid="${p}_${s}"; tbam="$(findbam "$sid")"
    [[ -n "$nbam" && -n "$tbam" ]] || { echo "[WARN] BAM missing for $sid"; continue; }
    mkdir -p "$OUT/$sid"
    pile="$OUT/$sid/$sid.snp_pileup.csv.gz"
    if [[ ! -s "$pile" ]]; then
      snp-pileup -g -q15 -Q20 -P100 -r25,0 "$FACETS_SNP_VCF" "$pile" "$nbam" "$tbam"
    fi
    Rscript "$PROJECT_DIR/scripts/05_facets.R" "$pile" "$OUT/$sid" "$sid" "$FACETS_CVAL"
  done
done

# Combine per-sample purity/ploidy into one table
python3 - "$OUT" "$PROJECT_DIR/07_cnv_purity/purity_ploidy.tsv" <<'PY'
import sys, pathlib, pandas as pd
d = pathlib.Path(sys.argv[1]); rows = [pd.read_csv(f, sep="\t") for f in d.glob("*/*.purity_ploidy.tsv")]
pd.concat(rows).sort_values("sample_id").to_csv(sys.argv[2], sep="\t", index=False) if rows else print("no results")
PY
