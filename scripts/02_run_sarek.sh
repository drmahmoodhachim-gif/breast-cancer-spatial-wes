#!/usr/bin/env bash
# Run nf-core/sarek: fastp -> BWA-MEM2 -> MarkDuplicates -> BQSR -> Mutect2 (joint T+J vs N),
# Strelka2 (T vs N, J vs N), CNVkit, NGSCheckMate, VEP.
# Re-running resumes from cache (-resume). Logs + pipeline_info are kept for the Methods.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$PROJECT_DIR/config/project.env"

[[ -s "$TARGET_BED_PADDED" ]] || { echo "[ERROR] Run step 'bed' first ($TARGET_BED_PADDED missing)"; exit 1; }
command -v nextflow >/dev/null || { echo "[ERROR] nextflow not found (see envs/README)"; exit 1; }

mkdir -p "$SAREK_OUTDIR" "$NXF_WORK_DIR"
cd "$PROJECT_DIR/05_sarek"

nextflow run nf-core/sarek \
  -r "$SAREK_VERSION" \
  -profile "$NXF_PROFILE" \
  -params-file "$PROJECT_DIR/config/sarek_params.yaml" \
  -c "$PROJECT_DIR/config/sarek_custom.config" \
  -work-dir "$NXF_WORK_DIR" \
  --outdir "$SAREK_OUTDIR" \
  -resume \
  "$@"

# Freeze exact software versions for the Methods section
find "$SAREK_OUTDIR/pipeline_info" -name "*software*versions*.yml" -exec cp {} "$PROJECT_DIR/13_paper/software_versions_sarek.yml" \; 2>/dev/null || true
echo "[sarek] finished. MultiQC report: $SAREK_OUTDIR/multiqc/multiqc_report.html"
