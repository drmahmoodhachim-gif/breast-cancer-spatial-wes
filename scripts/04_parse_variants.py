#!/usr/bin/env python3
"""Parse joint Mutect2 (T+J vs N, VEP-annotated) VCFs + Strelka2 support into tidy tables.

For every patient:
  * reads the joint Mutect2 VCF (prefers the VEP-annotated copy)
  * extracts per-sample depth / alt reads / VAF for N, J and T
  * annotates Strelka2 PASS support per pair (confidence tier, NOT a requirement)
  * applies config/filters.yaml call filters
  * classifies each variant per specimen: present / absent (powered) / not_evaluable
  * assigns a spatial class: shared, T_private, J_private, or *_unresolved

Outputs (06_variants/):
  somatic_all.tsv        every Mutect2 record (incl. filtered) with counts
  somatic_filtered.tsv   variants passing call filters (the analysis set)
  variant_counts.tsv     per-patient counts by class and tier
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import (PATIENTS, PROJECT_DIR, driver_genes, load_env, load_filters,  # noqa: E402
                    parse_sample_id, write_provenance)

env = load_env()
F = load_filters()
RES = Path(env["SAREK_OUTDIR"])
OUT = PROJECT_DIR / "06_variants"
DRIVERS = set(driver_genes())
CODING = set(F["coding_consequences"])


def find_mutect2_vcf(patient: str) -> Path | None:
    cands = [p for p in RES.rglob(f"*{patient}*.vcf.gz")
             if "mutect2" in str(p).lower() and "_vs_" not in p.name and "filtered" in p.name.lower()]
    ann = [p for p in cands if "vep" in p.name.lower() or "ann" in p.name.lower()]
    return sorted(ann or cands)[0] if (ann or cands) else None


def strelka_pass_keys(patient: str, specimen: str) -> set:
    keys = set()
    for p in RES.rglob(f"*{patient}_{specimen}_vs_*somatic*.vcf.gz"):
        if "strelka" not in str(p).lower():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                if f[6] == "PASS":
                    for alt in f[4].split(","):
                        keys.add((f[0], int(f[1]), f[3], alt))
    return keys


def vep_allele(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return alt
    if ref[0] == alt[0]:
        return alt[1:] or "-"
    return alt


AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
    "Sec": "U", "Pyl": "O", "Xaa": "X",
}
PIK3CA_HOTSPOT = {"E542K", "E545K", "H1047R", "H1047L"}


def csq_get(csq: dict, name: str) -> str:
    """Case-sensitive lookup of a CSQ field named in the VCF header."""
    return str(csq.get(name, "") or "")


def hgvsp_short(hgvsp: str) -> str:
    """Drop the transcript prefix and write amino acids as one letter."""
    import re
    text = hgvsp.split(":")[-1].replace("%3D", "=")
    return re.sub(r"[A-Z][a-z]{2}", lambda m: AA3.get(m.group(0), m.group(0)), text)


def hotspot_label(gene: str, short: str) -> str:
    import re
    change = short[2:] if short.startswith("p.") else short
    if gene == "PIK3CA" and change in PIK3CA_HOTSPOT:
        return f"PIK3CA_{change}"
    if gene == "TP53":
        nums = [int(n) for n in re.findall(r"(\d+)", change)]
        if any(102 <= n <= 292 for n in nums):
            return "TP53_DNA_binding_domain"
    return ""


def parse_vcf(path: Path, patient: str) -> pd.DataFrame:
    csq_fields, samples, recs = None, [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = line.split("Format: ")[1].split('"')[0].split("|")
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, filt, info, fmt = f[0], int(f[1]), f[3], f[4].split(","), f[6], f[7], f[8].split(":")
            info_d = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in info.split(";"))
            sdata = {}
            for sname, sval in zip(samples, f[9:]):
                _, spec = parse_sample_id(sname)
                if spec is None:
                    continue
                d = dict(zip(fmt, sval.split(":")))
                ad = [int(x) if x not in (".", "") else 0 for x in d.get("AD", "0").split(",")]
                sdata[spec] = ad
            csq = [dict(zip(csq_fields, c.split("|"))) for c in info_d.get("CSQ", "").split(",")] if csq_fields and "CSQ" in info_d else []
            for i, alt in enumerate(alts):
                rec = {"patient": patient, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                       "var_id": f"{patient}:{chrom}:{pos}:{ref}>{alt}", "filter": filt,
                       "variant_type": "SNV" if len(ref) == len(alt) == 1 else ("MNV" if len(ref) == len(alt) else "INDEL")}
                for spec in ("N", "J", "T"):
                    ad = sdata.get(spec)
                    if ad is None:
                        rec.update({f"{spec}_depth": np.nan, f"{spec}_alt": np.nan, f"{spec}_vaf": np.nan})
                        continue
                    depth = ad[0] + (ad[i + 1] if i + 1 < len(ad) else 0)
                    a = ad[i + 1] if i + 1 < len(ad) else 0
                    rec.update({f"{spec}_depth": depth, f"{spec}_alt": a, f"{spec}_vaf": a / depth if depth else np.nan})
                sel = [c for c in csq if c.get("Allele") == vep_allele(ref, alt)] or csq
                can = [c for c in sel if c.get("CANONICAL") == "YES"] or sel
                c = can[0] if can else {}
                rec.update({
                    "gene": csq_get(c, "SYMBOL"),
                    "consequence": csq_get(c, "Consequence"),
                    "impact": csq_get(c, "IMPACT"),
                    "hgvsc": csq_get(c, "HGVSc").split(":")[-1].replace("%3D", "="),
                    "hgvsp": csq_get(c, "HGVSp").split(":")[-1].replace("%3D", "="),
                    "existing_variation": csq_get(c, "Existing_variation"),
                    "clinvar_sig": csq_get(c, "CLIN_SIG").replace("&", ","),
                    "protein_position": csq_get(c, "Protein_position"),
                    "amino_acids": csq_get(c, "Amino_acids"),
                    "protein_domain": csq_get(c, "DOMAINS").replace("&", ",").replace("%3A", ":"),
                })
                rec["hgvsp_short"] = hgvsp_short(rec["hgvsp"])
                rec["hotspot"] = hotspot_label(rec["gene"], rec["hgvsp_short"])
                recs.append(rec)
    return pd.DataFrame(recs)


def detection_power(depth, vaf_expected, min_alt):
    if not depth or np.isnan(depth) or np.isnan(vaf_expected) or vaf_expected <= 0:
        return 0.0
    return float(binom.sf(min_alt - 1, int(depth), vaf_expected))


def absence_threshold(depth, floor, min_alt):
    """Smallest alt count that meets both alt>=min_alt and VAF>=floor. D065."""
    if not depth or np.isnan(depth) or depth <= 0:
        return min_alt
    return max(min_alt, int(np.ceil(float(floor) * float(depth) - 1e-9)))


def classify(df: pd.DataFrame) -> pd.DataFrame:
    floor = F["min_vaf_present"]
    min_alt = F["min_alt_reads_present"]
    for s in ("T", "J", "N"):
        df[f"{s}_present"] = (df[f"{s}_alt"] >= min_alt) & (df[f"{s}_vaf"] >= floor)
    states = {}
    reads_only = {}
    for s, other in (("T", "J"), ("J", "T")):
        st, old = [], []
        for _, r in df.iterrows():
            if r[f"{s}_present"]:
                st.append("present")
                old.append("present")
                continue
            other_v = r[f"{other}_vaf"]
            other_v = 0.0 if np.isnan(other_v) else float(other_v)
            exp = max(0.5 * other_v, floor)
            depth = r[f"{s}_depth"]
            pw_reads = detection_power(depth, exp, min_alt)
            pw_joint = detection_power(depth, exp, absence_threshold(depth, floor, min_alt))
            old.append("absent" if pw_reads >= F["power_threshold"] else "not_evaluable")
            st.append("absent" if pw_joint >= F["power_threshold"] else "not_evaluable")
        states[s] = st
        reads_only[s] = old
    df["T_state"], df["J_state"] = states["T"], states["J"]
    df["T_state_reads_only"], df["J_state_reads_only"] = reads_only["T"], reads_only["J"]

    def cls_from(t, j):
        if t == "present" and j == "present":
            return "shared"
        if t == "present":
            return "T_private" if j == "absent" else "T_unresolved"
        if j == "present":
            return "J_private" if t == "absent" else "J_unresolved"
        return "neither"

    df["spatial_class"] = [cls_from(t, j) for t, j in zip(df["T_state"], df["J_state"])]
    df["spatial_class_reads_only"] = [cls_from(t, j) for t, j in zip(df["T_state_reads_only"], df["J_state_reads_only"])]
    return df


all_rows, inputs = [], []
for p in PATIENTS:
    vcf = find_mutect2_vcf(p)
    if vcf is None:
        print(f"[WARN] no joint Mutect2 VCF found for {p}")
        continue
    inputs.append(vcf)
    df = parse_vcf(vcf, p)
    for spec in ("T", "J"):
        keys = strelka_pass_keys(p, spec)
        df[f"strelka_{spec}"] = [(r.chrom, r.pos, r.ref, r.alt) in keys for r in df.itertuples()]
    all_rows.append(df)
    print(f"[variants] {p}: {len(df)} records from {vcf.name}")

if not all_rows:
    sys.exit("[ERROR] no Mutect2 VCFs parsed")

v = pd.concat(all_rows, ignore_index=True)
hgvs_rows = []
for p in PATIENTS:
    for name in (f"{p}.pass.vep.vcf.gz", f"{p}.rescued.vep.vcf.gz"):
        hgvs = PROJECT_DIR / "06_variants" / "vep_hgvs" / name
        if hgvs.exists():
            hgvs_rows.append(parse_vcf(hgvs, p))
            inputs.append(hgvs)
if hgvs_rows:
    h = pd.concat(hgvs_rows, ignore_index=True)
    h = h[["chrom", "pos", "ref", "alt", "hgvsc", "hgvsp", "hgvsp_short", "hotspot", "clinvar_sig",
           "protein_position", "amino_acids", "protein_domain"]].drop_duplicates(["chrom", "pos", "ref", "alt"])
    v = v.drop(columns=[c for c in h.columns if c not in ("chrom", "pos", "ref", "alt")], errors="ignore")
    v = v.merge(h, on=["chrom", "pos", "ref", "alt"], how="left")
v["revel"] = "NA"
v["alphamissense"] = "NA"
if "protein_domain" not in v.columns:
    v["protein_domain"] = ""
v["is_driver_gene"] = v["gene"].isin(DRIVERS)
v["is_coding"] = v["consequence"].fillna("").apply(lambda c: any(x in CODING for x in c.split("&")))
v.to_csv(OUT / "somatic_all.tsv", sep="\t", index=False)

def only_normal_artifact(val):
    parts = [x for x in str(val).split(";") if x and x != "."]
    return parts == ["normal_artifact"]


v["only_normal_artifact"] = v["filter"].map(only_normal_artifact)
pass_mask = v["filter"].eq("PASS") if F["require_filter_pass"] else pd.Series(True, index=v.index)
tumor_ok = pd.Series(False, index=v.index)
for s in ("T", "J"):
    tumor_ok = tumor_ok | ((v[f"{s}_depth"] >= F["min_depth_tumor_side"]) & (v[f"{s}_alt"] >= F["min_alt_reads_call"])
                            & (v[f"{s}_vaf"] >= F["min_vaf_call"]))
depth_ok = tumor_ok & (v["N_depth"] >= F["min_depth_normal"])
keep = v[(pass_mask | v["only_normal_artifact"]) & depth_ok].copy()
keep["call_set"] = np.where(keep["filter"].eq("PASS"), "PASS", "rescued_normal_artifact")
rescue_n = (keep["call_set"] == "rescued_normal_artifact").groupby(keep["patient"]).sum()
only_n = v["only_normal_artifact"].groupby(v["patient"]).sum()
pd.DataFrame({
    "only_normal_artifact": only_n,
    "rescued_passing_filters": rescue_n,
}).fillna(0).astype(int).rename_axis("patient").to_csv(OUT / "normal_artifact_rescue.tsv", sep="\t")
print("[variants] normal_artifact rescue")
print(pd.read_csv(OUT / "normal_artifact_rescue.tsv", sep="\t").to_string(index=False))
keep = classify(keep)
keep["tier"] = np.where((keep["T_present"] & keep["strelka_T"]) | (keep["J_present"] & keep["strelka_J"]),
                        "tier1_mutect2+strelka2", "tier2_mutect2_only")
primary = keep[keep["call_set"].eq("PASS")].copy()
rescued_raw = keep[keep["call_set"].eq("rescued_normal_artifact")].copy()
primary.to_csv(OUT / "somatic_filtered.tsv", sep="\t", index=False)
rescued_raw.to_csv(OUT / "rescued_normal_artifact_raw.tsv", sep="\t", index=False)
keep = primary

counts = keep.groupby(["patient", "spatial_class", "tier"]).size().unstack(fill_value=0)
counts.to_csv(OUT / "variant_counts.tsv", sep="\t")
for f in ("somatic_all.tsv", "somatic_filtered.tsv"):
    write_provenance(OUT / f, inputs, F)
print(keep.groupby(["patient", "spatial_class"]).size().unstack(fill_value=0))
