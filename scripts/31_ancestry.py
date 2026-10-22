#!/usr/bin/env python3
"""1000 Genomes superpopulation for each adjacent normal, via peddy."""
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import peddy


_ORIG_FROMSTRING = np.fromstring


def _numpy2_fromstring(string, dtype=float, count=-1, sep="", *, like=None):
    """peddy 0.4.8 still calls the removed binary mode of np.fromstring."""
    if isinstance(sep, str) and sep != "":
        return _ORIG_FROMSTRING(string, dtype=dtype, count=count, sep=sep)
    buf = string if isinstance(string, (bytes, bytearray, memoryview)) else memoryview(string)
    arr = np.frombuffer(buf, dtype=dtype)
    if count != -1:
        arr = arr[:count]
    return np.array(arr, copy=True)

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PATIENTS, PROJECT_DIR, load_env, write_provenance  # noqa: E402

env = load_env()
RES = Path(env["SAREK_OUTDIR"])
FASTA = env["REFERENCE_FASTA"]
BCFTOOLS = env.get("BCFTOOLS") or shutil.which("bcftools") or "bcftools"
SITES = Path(env.get("PEDDY_SITES") or Path(peddy.__file__).resolve().parent / "GRCH38.sites")
OUT = PROJECT_DIR / "16_ancestry"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    targets = OUT / "peddy_hg38.targets.tsv"
    with SITES.open() as src, targets.open("w") as dst:
        for line in src:
            chrom, pos, _ref, _alt = line.strip().split(":")
            dst.write(f"chr{chrom}\t{pos}\n")
    bams = OUT / "normals.bamlist"
    with bams.open("w") as fh:
        for patient in PATIENTS:
            sid = f"{patient}_N"
            fh.write(str(RES / "preprocessing" / "recalibrated" / sid / f"{sid}.recal.bam") + "\n")
    vcf = OUT / "normals.vcf.gz"
    if not vcf.exists() or vcf.stat().st_size < 1000:
        subprocess.check_call(
            f"{BCFTOOLS} mpileup -f {FASTA} -b {bams} -T {targets} -q 20 -Q 20 "
            f"-a FORMAT/AD,FORMAT/DP -Ou | {BCFTOOLS} call -m -Oz -o {vcf}",
            shell=True,
        )
        subprocess.check_call([BCFTOOLS, "index", "-f", str(vcf)])
    samples = subprocess.check_output([BCFTOOLS, "query", "-l", str(vcf)], text=True).split()
    rename = OUT / "sample_rename.txt"
    with rename.open("w") as fh:
        for name in samples:
            short = name
            for patient in PATIENTS:
                if name.endswith(f"_{patient}_N") or name == f"{patient}_N":
                    short = f"{patient}_N"
            fh.write(f"{name}\t{short}\n")
    called = OUT / "normals.named.vcf.gz"
    subprocess.check_call([BCFTOOLS, "reheader", "-s", str(rename), "-o", str(called), str(vcf)])
    subprocess.check_call([BCFTOOLS, "index", "-f", str(called)])
    ped = OUT / "normals.ped"
    with ped.open("w") as fh:
        fh.write("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\n")
        for patient in PATIENTS:
            fh.write(f"{patient}\t{patient}_N\t0\t0\t2\t-9\n")
    np.fromstring = _numpy2_fromstring
    from peddy.cli import peddy as peddy_cmd
    peddy_cmd.main(args=[
        "--sites", "hg38", "--plot", "--prefix", str(OUT / "peddy"), str(called), str(ped),
    ], standalone_mode=False)
    het = pd.read_csv(OUT / "peddy.het_check.csv")
    keep = [c for c in ("sample_id", "ancestry-prediction", "ancestry-prob", "PC1", "PC2", "PC3", "PC4") if c in het.columns]
    table = het[keep].copy()
    table.insert(0, "patient", table["sample_id"].str.replace("_N", "", regex=False))
    dest = OUT / "ancestry_superpopulation.tsv"
    table.to_csv(dest, sep="\t", index=False)
    write_provenance(dest, [], {"tool": "peddy 0.4.8", "sites": "peddy GRCH38.sites", "panel": "1000 Genomes", "samples": "adjacent normal only"})
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
