#!/usr/bin/env python3
"""Allele-aware, left-aligned recount of audited indels and the named drivers.

Does not edit the variant table. IGV is used only if the igv executable exists.
Otherwise the snapshots are matplotlib pileups and are labelled as such.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[1]


def _env_value(key):
    env_file = ROOT / "config" / "project.env"
    if env_file.exists():
        for line in env_file.read_text().splitlines():
            if line.startswith(key + "="):
                return line.split("=", 1)[1].strip().strip('"')
    return ""


BAM_ROOT = Path(_env_value("SAREK_OUTDIR") or ROOT / "sarek_results") / "preprocessing" / "recalibrated"
SAMTOOLS = shutil.which("samtools") or "samtools"
OUT = ROOT / "00_admin" / "indel_allele_recount.tsv"
SNAP = ROOT / "00_admin" / "pileup_snapshots"
SNAP.mkdir(parents=True, exist_ok=True)


def fasta():
    for line in (ROOT / "config" / "project.env").read_text().splitlines():
        if line.startswith("REFERENCE_FASTA="):
            return line.split("=", 1)[1].strip().strip('"')
    return ""


def bam_path(sample):
    return BAM_ROOT / sample / f"{sample}.recal.bam"


def call_allele(read, pos0, ref, alt, min_bq):
    ref = ref.upper()
    alt = alt.upper()
    end = pos0 + len(ref)
    seq = []
    covered = [False] * len(ref)
    started = False
    for qpos, rpos in read.get_aligned_pairs():
        if rpos is not None and rpos >= end:
            break
        if rpos is not None and rpos < pos0:
            continue
        if qpos is None:
            if rpos is not None and pos0 <= rpos < end:
                covered[rpos - pos0] = True
            continue
        base = read.query_sequence[qpos] if read.query_sequence else None
        bq = read.query_qualities[qpos] if read.query_qualities is not None else 60
        if base is None:
            return None
        if rpos is None:
            if started:
                if bq < min_bq:
                    return None
                seq.append(base)
            continue
        if pos0 <= rpos < end:
            started = True
            covered[rpos - pos0] = True
            if bq < min_bq:
                return None
            seq.append(base)
    if not all(covered):
        return None
    observed = "".join(seq).upper()
    if observed == ref:
        return "ref"
    if observed == alt:
        return "alt"
    return "other"


def recount(sample, chrom, pos, ref, alt, min_mapq, min_bq):
    path = bam_path(sample)
    if not path.exists():
        return {"error": "bam missing", "alt": 0, "ref": 0, "other": 0, "depth": 0}
    chrom = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
    pos0 = int(pos) - 1
    alt_n = ref_n = other_n = 0
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(chrom, pos0, pos0 + max(len(ref), 1)):
            if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if read.mapping_quality < min_mapq:
                continue
            allele = call_allele(read, pos0, ref, alt, min_bq)
            if allele == "alt":
                alt_n += 1
            elif allele == "ref":
                ref_n += 1
            elif allele == "other":
                other_n += 1
    return {"alt": alt_n, "ref": ref_n, "other": other_n, "depth": alt_n + ref_n + other_n, "error": ""}


def snapshot(sample, chrom, pos, ref, alt, gene):
    path = bam_path(sample)
    if not path.exists():
        return ""
    chrom = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
    pos0 = int(pos) - 1
    rows = []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(chrom, max(0, pos0 - 5), pos0 + len(ref) + 5):
            if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            allele = call_allele(read, pos0, ref, alt, 0)
            if allele is None:
                continue
            rows.append(allele)
            if len(rows) >= 40:
                break
    if not rows:
        return ""
    fig, ax = plt.subplots(figsize=(4.2, 2.4))
    colors = {"ref": "#4C78A8", "alt": "#E45756", "other": "#B0B0B0"}
    for i, allele in enumerate(rows):
        ax.barh(i, 1, color=colors[allele], height=0.8)
    ax.set_yticks([])
    ax.set_xlim(0, 1)
    ax.set_xlabel("spanning reads, first 40")
    ax.set_title(f"{gene} {sample} {chrom}:{pos}\nnot IGV", fontsize=8)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=k) for k, c in colors.items()], fontsize=6, loc="upper right")
    out = SNAP / f"{gene}_{sample}_{chrom}_{pos}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return str(out)


def main():
    follow = pd.read_csv(ROOT / "00_admin" / "bam_recount_followup.tsv", sep="\t")
    V = pd.read_csv(ROOT / "06_variants" / "somatic_filtered.tsv", sep="\t")
    sites = []
    for site in follow["site"].unique():
        hit = V[V["var_id"] == site]
        if hit.empty:
            continue
        r = hit.iloc[0]
        sites.append((r.patient, r.chrom, int(r.pos), r.ref, r.alt, r.gene if pd.notna(r.gene) else site, "audited"))
    named = V[V["hgvsc"].isin(["c.628dup", "c.7286_7321+5del"]) | (
        (V["gene"] == "KMT2C") & (V["patient"] == "DH00023")
    ) | (
        (V["gene"] == "PIK3CA") & (V["hgvsp_short"].isin(["p.E545K", "p.H1047R"]))
    ) | (
        (V["gene"] == "TP53") & (V["patient"] == "DH00022")
    )]
    for r in named.itertuples(index=False):
        sites.append((r.patient, r.chrom, int(r.pos), r.ref, r.alt, r.gene if pd.notna(r.gene) else "", "driver"))
    # unique
    seen = set()
    uniq = []
    for s in sites:
        key = (s[0], s[1], s[2], s[3], s[4])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(s)
    rows = []
    igv = shutil.which("igv") or shutil.which("igv.sh")
    for patient, chrom, pos, ref, alt, gene, kind in uniq:
        specimens = ("N", "J", "T") if kind == "driver" or True else ("J",)
        # audited indels: recount every specimen that has a follow-up row, plus N/J/T for drivers
        wanted = ["N", "J", "T"] if kind == "driver" else sorted(set(
            follow.loc[follow["site"] == f"{patient}:{chrom}:{pos}:{ref}>{alt}", "specimen"]
        ) | {"N", "J", "T"})
        for spec in wanted:
            sample = f"{patient}_{spec}"
            for label, mq, bq in (("no_quality_filter", 0, 0), ("mapq20_baseq20", 20, 20)):
                rec = recount(sample, chrom, pos, ref, alt, mq, bq)
                table = V[(V["patient"] == patient) & (V["chrom"] == chrom) & (V["pos"] == pos)]
                rows.append({
                    "patient": patient,
                    "specimen": spec,
                    "gene": gene,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "kind": kind,
                    "filter": label,
                    "table_alt": int(table.iloc[0][f"{spec}_alt"]) if len(table) else "",
                    "table_depth": int(table.iloc[0][f"{spec}_depth"]) if len(table) else "",
                    "recount_alt": rec["alt"],
                    "recount_ref": rec["ref"],
                    "recount_other": rec["other"],
                    "recount_depth": rec["depth"],
                    "matches_table_alt": rec["alt"] == int(table.iloc[0][f"{spec}_alt"]) if len(table) else False,
                    "error": rec["error"],
                    "igv_available": bool(igv),
                })
        if kind == "driver" or gene == "chr9_64407964":
            for spec in ("J", "T"):
                snapshot(f"{patient}_{spec}", chrom, pos, ref, alt, f"{gene}_{str(ref)[:8]}")
    # Recurrent gnomAD site. Review only; it stays in the primary table.
    for patient in ("DH00022", "DH00024"):
        hit = V[(V["patient"] == patient) & (V["pos"] == 64407964)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        for spec in ("N", "J", "T"):
            sample = f"{patient}_{spec}"
            rec = recount(sample, r.chrom, int(r.pos), r.ref, r.alt, 0, 0)
            rows.append({
                "patient": patient, "specimen": spec, "gene": r.gene if pd.notna(r.gene) else "chr9:64407964",
                "chrom": r.chrom, "pos": int(r.pos), "ref": r.ref, "alt": r.alt,
                "kind": "gnomad_recurrent", "filter": "no_quality_filter",
                "table_alt": int(r[f"{spec}_alt"]), "table_depth": int(r[f"{spec}_depth"]),
                "recount_alt": rec["alt"], "recount_ref": rec["ref"], "recount_other": rec["other"],
                "recount_depth": rec["depth"], "matches_table_alt": rec["alt"] == int(r[f"{spec}_alt"]),
                "error": rec["error"], "igv_available": bool(igv),
            })
            snapshot(sample, r.chrom, int(r.pos), r.ref, r.alt, "chr9_64407964")
    # The one SNV, unfiltered mpileup, kept separate from the allele counter.
    ref_fa = fasta()
    snv_bam = bam_path("DH00018_J")
    if ref_fa and snv_bam.exists():
        cmd = [SAMTOOLS, "mpileup", "-f", ref_fa, "-r", "chr11:85664329-85664329", str(snv_bam)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        (ROOT / "00_admin" / "dh00018_chr11_85664329_mpileup.txt").write_text(proc.stdout + "\n" + proc.stderr)
        rows.append({
            "patient": "DH00018", "specimen": "J", "gene": "", "chrom": "chr11", "pos": 85664329,
            "ref": "C", "alt": "A", "kind": "snv_mpileup_no_q", "filter": "mpileup without -q -Q",
            "table_alt": 23, "table_depth": 105, "recount_alt": "", "recount_ref": "",
            "recount_other": "", "recount_depth": "", "matches_table_alt": "",
            "error": proc.stdout.strip().replace("\t", " ")[:500], "igv_available": bool(igv),
        })
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    print("igv", bool(igv))
    print("rows", len(rows))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
