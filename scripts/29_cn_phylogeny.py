#!/usr/bin/env python3
"""FACETS breakpoint sharing and, when installed, a MEDICC2 copy-number tree."""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

FAC = PROJECT_DIR / "07_cnv_purity" / "facets"
OUT = PROJECT_DIR / "07_cnv_purity" / "cn_phylogeny"
OUT.mkdir(parents=True, exist_ok=True)
segs = pd.concat(pd.read_csv(f, sep="\t") for f in FAC.glob("*/*.segments.tsv"))
segs["chrom"] = segs["chrom"].astype(str).str.replace("chr", "", regex=False)
auto = {str(i) for i in range(1, 23)}
segs = segs[segs["chrom"].isin(auto)].copy()


def breakpoints(g):
    pts = []
    for chrom, cg in g.groupby("chrom"):
        cg = cg.sort_values("start")
        prev = None
        for _, r in cg.iterrows():
            tcn = None if pd.isna(r["tcn.em"]) else int(round(r["tcn.em"]))
            lcn = None if pd.isna(r["lcn.em"]) else int(round(r["lcn.em"]))
            state = (tcn, lcn)
            if prev is not None and state != prev:
                pts.append((str(chrom), int(r["start"])))
            prev = state
    return pts


rows = []
for patient in sorted(set(segs["sample_id"].str.replace(r"_[JT]$", "", regex=True))):
    t = segs[segs["sample_id"] == f"{patient}_T"]
    j = segs[segs["sample_id"] == f"{patient}_J"]
    tb, jb = breakpoints(t), breakpoints(j)
    near = 0
    for chrom, pos in tb:
        if any(c == chrom and abs(pos - jp) <= 1_000_000 for c, jp in jb):
            near += 1
    rows.append({
        "patient": patient,
        "n_T_breakpoints": len(tb),
        "n_J_breakpoints": len(jb),
        "n_T_within_1Mb_of_J": near,
        "fraction_T_within_1Mb_of_J": (near / len(tb) if tb else float("nan")),
    })
share = pd.DataFrame(rows)
share_path = OUT / "breakpoint_sharing.tsv"
share.to_csv(share_path, sep="\t", index=False)
write_provenance(share_path, [], {"window": "1 Mb", "breakpoints": "segment starts where allele-specific CN changes"})
print(share.round(3).to_string(index=False))

IMAGE = "quay.io/biocontainers/medicc2:1.1.2--py310h8ea774a_2"
medicc = shutil.which("medicc2")
docker = shutil.which("docker")
if medicc is None and docker is None:
    print("[cn] medicc2 is not on PATH and docker is unavailable; tree not fitted")
    (OUT / "medicc2_status.txt").write_text("medicc2 not installed and docker unavailable\n")
else:
    (OUT / "medicc2_status.txt").write_text(
        f"binary {medicc}\n" if medicc else f"docker image {IMAGE}\n"
    )
    for patient in share["patient"]:
        by_spec = {spec: segs[segs["sample_id"] == f"{patient}_{spec}"] for spec in ("J", "T")}

        def allele_at(spec, chrom, pos):
            g = by_spec[spec]
            hit = g[(g["chrom"] == chrom) & (g["start"] <= pos) & (g["end"] >= pos)]
            if hit.empty or pd.isna(hit.iloc[0]["tcn.em"]) or pd.isna(hit.iloc[0]["lcn.em"]):
                return None
            r = hit.iloc[0]
            minor = int(round(r["lcn.em"]))
            major = int(round(r["tcn.em"] - r["lcn.em"]))
            return minor, max(major, 0)

        pieces = []
        chroms = sorted(set(by_spec["J"]["chrom"]) | set(by_spec["T"]["chrom"]), key=int)
        for chrom in chroms:
            bounds = set()
            for spec in ("J", "T"):
                g = by_spec[spec][by_spec[spec]["chrom"] == chrom]
                bounds.update(int(x) for x in g["start"])
                bounds.update(int(x) for x in g["end"])
            bounds = sorted(bounds)
            for left, right in zip(bounds, bounds[1:]):
                if right <= left:
                    continue
                mid = (left + right) // 2
                cn_j, cn_t = allele_at("J", chrom, mid), allele_at("T", chrom, mid)
                if cn_j is None or cn_t is None:
                    continue
                for spec, cn in (("N", (1, 1)), ("J", cn_j), ("T", cn_t)):
                    pieces.append({"sample_id": f"{patient}_{spec}", "chrom": chrom,
                                   "start": left, "end": right, "cn_a": cn[0], "cn_b": cn[1]})
        inp = OUT / f"{patient}.medicc2_input.tsv"
        pd.DataFrame(pieces).to_csv(inp, sep="\t", index=False)
        dest = OUT / patient
        dest.mkdir(exist_ok=True)
        normal = f"{patient}_N"
        if medicc:
            cmd = [medicc, "-n", normal, "-j", "2", "--events", str(inp), str(dest)]
        else:
            cmd = ["docker", "run", "--rm", "-v", f"{OUT}:/data", IMAGE,
                   "medicc2", "-n", normal, "-j", "2", "--events",
                   f"/data/{inp.name}", f"/data/{patient}"]
        subprocess.check_call(cmd)
        print(f"[cn] medicc2 {patient}")
