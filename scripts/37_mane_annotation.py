#!/usr/bin/env python3
"""Rebuild gene, consequence and is_coding from the VEP 116 MANE/canonical transcript.

Does not re-run VEP. Reads 06_variants/vep_hgvs/<patient>.pass.vep.vcf.gz.
The chosen transcript is MANE Select when that field is set, otherwise the
canonical transcript. Protein HGVS is taken from that same transcript.
"""
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, driver_genes, load_filters, write_provenance  # noqa: E402

CODING = set(load_filters()["coding_consequences"])
DRIVERS = set(driver_genes())
AA = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


def protein(hgvs: str) -> str:
    from urllib.parse import unquote
    text = "" if hgvs is None or (isinstance(hgvs, float)) else str(hgvs)
    if text in ("", "nan", "None"):
        return ""
    text = unquote(text)
    if ":p." in text:
        return "p." + text.split(":p.")[-1]
    return text if text.startswith("p.") else ""


def cdna(hgvs: str) -> str:
    from urllib.parse import unquote
    text = "" if hgvs is None else unquote(str(hgvs))
    if ":c." in text:
        return "c." + text.split(":c.")[-1]
    if ":n." in text:
        return "n." + text.split(":n.")[-1]
    return text


def short_protein(p: str) -> str:
    if not p.startswith("p."):
        return ""
    body = p[2:]
    for long, one in sorted(AA.items(), key=lambda kv: -len(kv[0])):
        body = body.replace(long, one)
    return "p." + body


def csq_fields(header_line: str) -> list[str]:
    fmt = header_line.split("Format: ")[1].split('">')[0]
    return fmt.split("|")


def load_vep(patient: str) -> dict:
    path = PROJECT_DIR / "06_variants" / "vep_hgvs" / f"{patient}.pass.vep.vcf.gz"
    fields = None
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = csq_fields(line)
                continue
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, _qual, _filt, info = line.rstrip("\n").split("\t")[:8]
            if "CSQ=" not in info:
                continue
            raw = info.split("CSQ=")[1].split(";")[0]
            blocks = []
            for piece in raw.split(","):
                parts = piece.split("|")
                blocks.append(dict(zip(fields, parts)))
            out[(chrom, int(pos), ref, alt)] = blocks
    return out


def choose(blocks: list[dict]) -> tuple[dict, str]:
    mane = [b for b in blocks if b.get("MANE_SELECT")]
    if mane:
        return mane[0], "MANE"
    canon = [b for b in blocks if b.get("CANONICAL") == "YES"]
    if canon:
        return canon[0], "canonical"
    return blocks[0], "only_transcript"


def blank(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value)
    return "" if text in ("nan", "None") else text


def main():
    src = PROJECT_DIR / "06_variants" / "somatic_filtered.tsv"
    v = pd.read_csv(src, sep="\t")
    if len(v) != 380:
        raise SystemExit(f"[mane] expected 380 rows, found {len(v)}")
    caches = {p: load_vep(p) for p in v["patient"].unique()}
    changes = []
    unmatched = []
    for i, r in v.iterrows():
        key = (str(r["chrom"]), int(r["pos"]), str(r["ref"]), str(r["alt"]))
        blocks = caches[r["patient"]].get(key)
        if blocks is None:
            alt_key = (key[0], key[1] - 1, key[2], key[3])
            blocks = caches[r["patient"]].get(alt_key)
        if not blocks:
            unmatched.append(r["var_id"])
            continue
        chosen, rule = choose(blocks)
        new_gene = chosen.get("SYMBOL") or ""
        new_csq = chosen.get("Consequence") or ""
        new_p = protein(chosen.get("HGVSp", ""))
        new_c = cdna(chosen.get("HGVSc", ""))
        terms = [t for t in new_csq.split("&") if t]
        new_coding = any(t in CODING for t in terms)
        old_gene = blank(r["gene"])
        old_csq = blank(r["consequence"])
        old_coding = bool(r["is_coding"]) if pd.notna(r["is_coding"]) else False
        old_p = protein(blank(r["hgvsp"]))
        changed = (new_gene != old_gene) or (new_csq != old_csq) or (new_coding != old_coding)
        v.at[i, "gene"] = new_gene
        v.at[i, "consequence"] = new_csq
        v.at[i, "impact"] = chosen.get("IMPACT") or blank(r["impact"])
        if new_c:
            v.at[i, "hgvsc"] = new_c
        if new_p:
            v.at[i, "hgvsp"] = new_p
            if new_p != old_p or blank(r["hgvsp_short"]) == "":
                v.at[i, "hgvsp_short"] = short_protein(new_p)
        if chosen.get("Protein_position"):
            v.at[i, "protein_position"] = chosen.get("Protein_position")
        if chosen.get("Amino_acids"):
            v.at[i, "amino_acids"] = chosen.get("Amino_acids")
        v.at[i, "is_coding"] = new_coding
        v.at[i, "is_driver_gene"] = new_gene in DRIVERS
        if changed:
            changes.append({
                "var_id": r["var_id"],
                "patient": r["patient"],
                "old_gene": old_gene,
                "new_gene": new_gene,
                "in_driver_genes": new_gene in DRIVERS,
                "old_consequence": old_csq,
                "new_consequence": new_csq,
                "old_is_coding": old_coding,
                "new_is_coding": new_coding,
                "old_hgvsp": old_p,
                "new_hgvsp": new_p,
                "transcript": chosen.get("Feature", ""),
                "transcript_rule": rule,
                "spatial_class": r["spatial_class"],
            })
    if unmatched:
        raise SystemExit("[mane] unmatched variants:\n" + "\n".join(unmatched))
    v.to_csv(src, sep="\t", index=False)
    write_provenance(src, [PROJECT_DIR / "06_variants" / "vep_hgvs" / f"{p}.pass.vep.vcf.gz" for p in caches],
                     {"rule": "MANE Select else canonical", "n_changed": len(changes)})
    ch = pd.DataFrame(changes)
    out = PROJECT_DIR / "06_variants" / "annotation_mane_changes.tsv"
    ch.to_csv(out, sep="\t", index=False)
    write_provenance(out, [src], {"n_changed": len(changes)})
    print(f"[mane] changed {len(ch)} of 380; coding now {int(v.is_coding.sum())}")
    if len(ch):
        print(ch.loc[ch["in_driver_genes"], ["var_id", "old_gene", "new_gene", "new_consequence"]].to_string(index=False))
        print("[mane] new genes", ", ".join(sorted(set(ch["new_gene"]) - {""})))

    ann = v.set_index("var_id")
    cols = ["gene", "consequence", "impact", "hgvsc", "hgvsp", "hgvsp_short",
            "protein_position", "amino_acids", "is_coding", "is_driver_gene"]
    for path in (
        PROJECT_DIR / "09_ccf_clonality" / "ccf.tsv",
        PROJECT_DIR / "09_ccf_clonality" / "DH00023_lineages.tsv",
    ):
        if not path.exists():
            continue
        t = pd.read_csv(path, sep="\t")
        if "var_id" not in t.columns:
            continue
        for col in cols:
            if col not in t.columns or col not in ann.columns:
                continue
            mapped = t["var_id"].map(ann[col])
            t[col] = mapped.where(mapped.notna(), t[col])
        t.to_csv(path, sep="\t", index=False)
        print(f"[mane] refreshed annotation columns in {path.name}")


if __name__ == "__main__":
    main()
