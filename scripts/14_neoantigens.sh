#!/usr/bin/env bash
# Neoantigen prediction (HLA class I) with pVACseq, per patient.
#  1. Take the joint Mutect2 VCF (T+J vs N), keep PASS calls that are in 06_variants/somatic_filtered.tsv
#  2. Re-annotate with VEP + Frameshift/Wildtype plugins (required by pVACseq), same cache release as sarek
#  3. pVACseq with the patient's OptiType class I alleles (typed from N), MHCflurry, 8–11-mers
#     VAF/coverage filters are switched OFF here: region assignment (T, J, shared) is done afterwards
#     from our own spatial classes, so junction-private variants are not discarded for lacking T reads.
# Output: 14_immune/neoantigens/<patient>/MHC_Class_I/*.all_epitopes.tsv and *.filtered.tsv
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$PROJECT_DIR/config/project.env"
export PATH="$(dirname "$PY"):${PATH}"
OUT="$PROJECT_DIR/14_immune/neoantigens"; mkdir -p "$OUT"
TYPES="$PROJECT_DIR/14_immune/hla/hla_types.tsv"
[[ -s "$TYPES" ]] || { echo "[ERROR] run HLA typing first ($TYPES)"; exit 1; }
[[ "$VEP_CACHE_DIR" != "__SET_ME__" && -d "$VEP_CACHE_DIR" ]] || { echo "[ERROR] set VEP_CACHE_DIR in config/project.env"; exit 1; }
[[ "$REFERENCE_FASTA" != "__SET_ME__" && -s "$REFERENCE_FASTA" ]] || { echo "[ERROR] set REFERENCE_FASTA in config/project.env"; exit 1; }

PLUG="${VEP_PLUGIN_DIR:-$OUT/vep_plugins}"
if [[ ! -s "$PLUG/Wildtype.pm" ]]; then
  mkdir -p "$PLUG"
  docker run --rm -v "$PLUG":/plug "$PVACTOOLS_IMAGE" pvacseq install_vep_plugin /plug
fi

"$PY" - "$PROJECT_DIR" <<'PY' > "$OUT/loci_by_patient.tsv"
import sys, pandas as pd
v = pd.read_csv(f"{sys.argv[1]}/06_variants/somatic_filtered.tsv", sep="\t")
v[["patient", "chrom", "pos"]].to_csv(sys.stdout, sep="\t", index=False, header=False)
PY

for p in DH00018 DH00019 DH00022 DH00023 DH00024; do
  vcf=$(find "$SAREK_OUTDIR/variant_calling" -path "*mutect2*" -name "${p}.mutect2.filtered.vcf.gz" | head -n1)
  [[ -n "$vcf" ]] || { echo "[WARN] no joint Mutect2 VCF for $p"; continue; }
  d="$OUT/$p"; mkdir -p "$d"
  awk -v p="$p" '$1==p {print $2"\t"$3}' "$OUT/loci_by_patient.tsv" > "$d/loci.tsv"
  [[ -s "$d/loci.tsv" ]] || { echo "[WARN] no filtered variants for $p"; continue; }
  bcftools view -f PASS -T "$d/loci.tsv" "$vcf" -Oz -o "$d/${p}.pass.vcf.gz"; tabix -f -p vcf "$d/${p}.pass.vcf.gz"

  docker run --rm -u "$(id -u):$(id -g)" -v "$d":/work -v "$VEP_CACHE_DIR":/cache -v "$PLUG":/plugins \
    -v "$(dirname "$REFERENCE_FASTA")":/ref "$VEP_IMAGE" vep \
    --input_file "/work/${p}.pass.vcf.gz" --output_file "/work/${p}.vep.vcf" --format vcf --vcf \
    --symbol --terms SO --tsl --mane_select --canonical --biotype --hgvs --pick \
    --fasta "/ref/$(basename "$REFERENCE_FASTA")" --offline --cache --dir_cache /cache \
    --cache_version "$VEP_CACHE_VERSION" --assembly GRCh38 --plugin Frameshift --plugin Wildtype --dir_plugins /plugins --force_overwrite
  bgzip -f "$d/${p}.vep.vcf"; tabix -f -p vcf "$d/${p}.vep.vcf.gz"

  tum=$(bcftools query -l "$d/${p}.vep.vcf.gz" | grep -E "${p}[_-]T" | head -n1)
  nor=$(bcftools query -l "$d/${p}.vep.vcf.gz" | grep -E "${p}[_-]N" | head -n1)
  alleles=$(awk -F'\t' -v s="${p}_N" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["sample_id"]==s{
      printf "HLA-%s,HLA-%s,HLA-%s,HLA-%s,HLA-%s,HLA-%s", $h["A1"],$h["A2"],$h["B1"],$h["B2"],$h["C1"],$h["C2"]}' "$TYPES" \
      | tr ',' '\n' | sort -u | paste -sd, -)
  [[ -n "$alleles" ]] || { echo "[WARN] no HLA alleles for ${p}_N"; continue; }
  echo "[neoag] $p tumour=$tum normal=$nor alleles=$alleles"

  docker run --rm -u "$(id -u):$(id -g)" -v "$d":/work "$PVACTOOLS_IMAGE" pvacseq run \
    "/work/${p}.vep.vcf.gz" "$tum" "$alleles" $NEOAG_ALGORITHMS /work/pvacseq \
    --normal-sample-name "$nor" -e1 8,9,10,11 --binding-threshold "$NEOAG_IC50_CUTOFF" \
    --tdna-vaf 0 --tdna-cov 0 --normal-vaf 1 --normal-cov 0 --pass-only -t "${THREADS:-4}"
done
echo "[neoag] done -> $OUT/<patient>/pvacseq/MHC_Class_I/"
