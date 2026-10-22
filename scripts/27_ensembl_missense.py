#!/usr/bin/env python3
"""Fill AlphaMissense, REVEL and Pfam domain for primary missense variants from Ensembl VEP."""
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

VAR = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
URL = ("https://rest.ensembl.org/vep/human/region"
       "?AlphaMissense=1&REVEL=1&canonical=1&domains=1")
V = pd.read_csv(VAR, sep="\t")
miss = V[V["consequence"].fillna("").str.contains("missense_variant")].copy()
print(f"[rest] missense {len(miss)} of {len(V)}")


def region(row):
    chrom = str(row.chrom).replace("chr", "")
    start = int(row.pos)
    end = start + len(str(row.ref)) - 1
    return f"{chrom} {start} {end} {row.ref}/{row.alt} 1"


keys = {}
for r in miss.itertuples(index=False):
    keys.setdefault(region(r), (str(r.chrom), int(r.pos), r.ref, r.alt))

payloads = list(keys)
scored = {}


def post(batch, attempt=1):
    body = json.dumps({"variants": batch}).encode()
    req = urllib.request.Request(
        URL, data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": "WES_TJN_project/missense-scores"},
    )
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            return json.loads(resp.read().decode())
    except Exception as exc:
        if attempt >= 3:
            raise
        print(f"[rest] retry after {exc}")
        time.sleep(5 * attempt)
        return post(batch, attempt + 1)


def pfam(conseq):
    names = []
    for dom in conseq.get("domains") or []:
        text = " ".join(str(dom.get(k, "")) for k in ("name", "db"))
        if "Pfam" in text:
            names.append(dom.get("name") or text)
    return ",".join(names)


for i in range(0, len(payloads), 200):
    batch = payloads[i:i + 200]
    print(f"[rest] batch {i // 200 + 1} n={len(batch)}")
    data = post(batch)
    for item in data:
        hit = None
        for c in item.get("transcript_consequences") or []:
            if c.get("canonical") == 1 or str(c.get("canonical")) == "1":
                hit = c
                break
        hit = hit or (item.get("transcript_consequences") or [None])[0]
        if not hit:
            continue
        key = keys.get(item.get("input"))
        if key is None:
            continue
        am = hit.get("alphamissense") or {}
        if not isinstance(am, dict):
            am = {}
        scored[key] = {
            "am_pathogenicity": am.get("am_pathogenicity", hit.get("am_pathogenicity", "")),
            "am_class": am.get("am_class", hit.get("am_class", "")),
            "revel": hit.get("revel", ""),
            "pfam_domain": pfam(hit),
        }
    time.sleep(1)

print(f"[rest] scored {len(scored)}")
V["am_pathogenicity"] = pd.NA
V["am_class"] = pd.NA
V["revel"] = pd.NA
V["pfam_domain"] = pd.NA
for idx, r in V.iterrows():
    rec = scored.get((str(r.chrom), int(r.pos), r.ref, r.alt))
    if not rec:
        continue
    V.at[idx, "am_pathogenicity"] = rec["am_pathogenicity"]
    V.at[idx, "am_class"] = rec["am_class"]
    V.at[idx, "revel"] = rec["revel"]
    V.at[idx, "pfam_domain"] = rec["pfam_domain"]
V.to_csv(VAR, sep="\t", index=False)
write_provenance(VAR, [], {"scores": "Ensembl REST VEP AlphaMissense=1 REVEL=1 canonical=1 domains=1"})
filled = lambda col: V[col].notna() & V[col].astype(str).str.len().gt(0)
print(f"[rest] rows with REVEL {int(filled('revel').sum())} AlphaMissense {int(filled('am_pathogenicity').sum())}")
