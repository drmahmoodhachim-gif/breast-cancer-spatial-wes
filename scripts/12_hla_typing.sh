#!/usr/bin/env bash
# HLA class I typing with OptiType (nf-core/hlatyping), memory- and disk-light version.
#
# Instead of converting whole BAMs to FASTQ (slow, ~36 GB requests), we first extract only the reads
# OptiType can use: the MHC region on chr6 (28.5–33.5 Mb), any HLA-* alt contigs present in the
# reference, and unmapped reads. D037: those reads come from the markduplicates BAM, which still
# holds off-target and unmapped reads. This sarek run wrote <sample>.md.bam, not CRAM.
# Output: 14_immune/hla/hla_types.tsv (+ concordance warnings T/J vs N). pVACseq is not started here.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$PROJECT_DIR/config/project.env"
OUT="$PROJECT_DIR/14_immune/hla"; FQ="$OUT/fastq_md"; mkdir -p "$OUT" "$FQ"
SHEET="$OUT/hlatyping_samplesheet_md.csv"
THREADS_X="${THREADS_HLA:-4}"
SAMTOOLS="${SAMTOOLS:-samtools}"
# Keep the recalibrated-BAM OptiType table for the comparison.
if [[ -s "$OUT/hla_types.tsv" && ! -s "$OUT/hla_types_recal.tsv" ]]; then
  cp "$OUT/hla_types.tsv" "$OUT/hla_types_recal.tsv"
fi

echo "sample,fastq_1,fastq_2,seq_type" > "$SHEET"
while IFS= read -r bam; do
  sid="$(basename "$bam" .md.bam)"
  r1="$FQ/${sid}_R1.fastq.gz"; r2="$FQ/${sid}_R2.fastq.gz"
  if [[ ! -s "$r1" ]]; then
    bed="$FQ/${sid}.regions.bed"
    echo -e "chr6\t28500000\t33500000" > "$bed"
    "$SAMTOOLS" idxstats "$bam" | awk '$1 ~ /^HLA[-_]/ {print $1"\t0\t"$2}' >> "$bed"
    echo "[hla] $sid: extracting MHC + HLA contigs + unmapped reads from markduplicates BAM"
    { "$SAMTOOLS" view -u -L "$bed" -@ "$THREADS_X" "$bam"; } > "$FQ/${sid}.mhc.bam"
    "$SAMTOOLS" view -u -f 4 -@ "$THREADS_X" "$bam" > "$FQ/${sid}.unmapped.bam"
    "$SAMTOOLS" cat "$FQ/${sid}.mhc.bam" "$FQ/${sid}.unmapped.bam" \
      | "$SAMTOOLS" collate -u -O -@ "$THREADS_X" - \
      | "$SAMTOOLS" fastq -@ "$THREADS_X" -F 0x900 -1 "$r1" -2 "$r2" -0 /dev/null -s /dev/null -n -
    rm -f "$FQ/${sid}.mhc.bam" "$FQ/${sid}.unmapped.bam"
  fi
  echo "${sid},${r1},${r2},dna" >> "$SHEET"
done < <(find "$SAREK_OUTDIR/preprocessing/markduplicates" -name "DH*_?.md.bam" | sort)
n=$(( $(wc -l < "$SHEET") - 1 )); echo "[hla] samples: $n"
[[ $n -gt 0 ]] || { echo "[ERROR] no markduplicates BAMs found"; exit 1; }

cd "$OUT"
export NXF_VER="${NXF_HLA_VER:-25.04.2}"
nextflow run nf-core/hlatyping -r "$HLATYPING_VERSION" -profile "$NXF_PROFILE" \
  -c "$PROJECT_DIR/config/hlatyping_custom.config" \
  -work-dir "$OUT/work_md" \
  --input "$SHEET" --outdir "$OUT/results_md" --genome GRCh38

PY="${PY:-python3}"
"$PY" - "$OUT/results_md" "$OUT/hla_types.tsv" "$OUT/hla_types_recal.tsv" "$OUT/hla_discordance_recal_vs_md.tsv" <<'PY'
import sys, glob, os, re, pandas as pd
rows = []
for f in glob.glob(os.path.join(sys.argv[1], "**", "*_result.tsv"), recursive=True):
    sid = re.search(r"(DH\d{5}_[NJT])", f)
    if not sid: continue
    d = pd.read_csv(f, sep="\t").iloc[0]
    rows.append({"sample_id": sid.group(1), "patient": sid.group(1)[:7],
                 **{k: d.get(k) for k in ["A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"]}})
if not rows: sys.exit("[hla] no OptiType results found")
alleles = ["A1", "A2", "B1", "B2", "C1", "C2"]
t = pd.DataFrame(rows).sort_values("sample_id")
t.to_csv(sys.argv[2], sep="\t", index=False)
print(t.to_string(index=False))
for p, g in t.groupby("patient"):
    n = g[g.sample_id.str.endswith("_N")]
    if n.empty: continue
    ref = set(n.iloc[0][alleles])
    for r in g.itertuples():
        diff = ref ^ set([r.A1, r.A2, r.B1, r.B2, r.C1, r.C2])
        if diff: print(f"[hla] WARNING {r.sample_id} differs from {p}_N: {sorted(map(str, diff))} (swap? HLA LOH?)")
old = pd.read_csv(sys.argv[3], sep="\t")
cmp_rows = []
for r in t.itertuples():
    prev = old[old.sample_id == r.sample_id]
    if prev.empty: continue
    p = prev.iloc[0]
    changed = [c for c in alleles if str(getattr(r, c)) != str(p[c])]
    old_n = old[(old.patient == r.patient) & old.sample_id.str.endswith("_N")].iloc[0]
    new_n = t[(t.patient == r.patient) & t.sample_id.str.endswith("_N")].iloc[0]
    is_n = str(r.sample_id).endswith("_N")
    old_diff = [] if is_n else sorted(set(map(str, old_n[alleles])) ^ set(map(str, p[alleles])))
    new_vals = [str(getattr(r, c)) for c in alleles]
    new_diff = [] if is_n else sorted(set(map(str, new_n[alleles])) ^ set(new_vals))
    resolved = [a for a in old_diff if a not in new_diff]
    cmp_rows.append({
        "sample_id": r.sample_id,
        "alleles_changed": ",".join(changed) if changed else "",
        "discordant_vs_N_recal": ",".join(old_diff),
        "discordant_vs_N_md": ",".join(new_diff),
        "resolved_alleles": ",".join(resolved),
    })
out = pd.DataFrame(cmp_rows)
out.to_csv(sys.argv[4], sep="\t", index=False)
print("[hla] allele changes versus the recalibrated-BAM run:")
print(out.to_string(index=False))
PY
