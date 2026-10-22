"""Shared helpers: project paths, settings, sample sheet, provenance records.

Every analysis script imports this so that paths and thresholds come from one place
(config/project.env and config/filters.yaml) and every output gets a provenance file.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import yaml

PROJECT_DIR = Path(__file__).resolve().parents[2]
CONFIG = PROJECT_DIR / "config"

PATIENTS = ["DH00018", "DH00019", "DH00022", "DH00023", "DH00024"]
SPECIMENS = ["N", "J", "T"]               # display order: normal -> junction -> tumour
TUMOR_SIDE = ["T", "J"]
SPECIMEN_LABEL = {"N": "Adjacent normal", "J": "Junction", "T": "Tumour"}


def load_env(path: Path = CONFIG / "project.env") -> dict:
    """Parse simple KEY="value" lines, expanding ${VAR} references."""
    env: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        val = val.split(" #")[0].strip().strip('"').strip("'")
        val = re.sub(r"\$\{(\w+)\}", lambda m: env.get(m.group(1), os.environ.get(m.group(1), "")), val)
        env[key.strip()] = val
    if "__PROJECT_DIR__" in env.get("PROJECT_DIR", "__PROJECT_DIR__"):
        env = {k: v.replace("__PROJECT_DIR__", str(PROJECT_DIR)) for k, v in env.items()}
        env["PROJECT_DIR"] = str(PROJECT_DIR)
    return env


def load_filters() -> dict:
    return yaml.safe_load((CONFIG / "filters.yaml").read_text())


def driver_genes() -> list[str]:
    lines = (CONFIG / "driver_genes.txt").read_text().splitlines()
    return [g.strip() for g in lines if g.strip() and not g.startswith("#")]


def sample_id(patient: str, specimen: str) -> str:
    return f"{patient}_{specimen}"


def parse_sample_id(name: str) -> tuple[str | None, str | None]:
    """Map any VCF/BAM sample name containing e.g. 'DH00018_T' or 'DH00018-T' to (patient, specimen)."""
    m = re.search(r"(DH\d{5})[_-]([NTJ])(?![A-Za-z0-9])", name)
    return (m.group(1), m.group(2)) if m else (None, None)


def git_hash() -> str:
    try:
        return subprocess.check_output(["git", "-C", str(PROJECT_DIR), "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "no-git"


def _md5(path: Path, limit: int = 200_000_000) -> str | None:
    if not path.is_file() or path.stat().st_size > limit:
        return None
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_provenance(output: Path, inputs: list[Path] | None = None, params: dict | None = None) -> None:
    """Write <output>.provenance.json next to every result: what made it, from what, when."""
    output = Path(output)
    rec = {
        "output": str(output.relative_to(PROJECT_DIR)) if output.is_relative_to(PROJECT_DIR) else str(output),
        "created": _dt.datetime.now().isoformat(timespec="seconds"),
        "script": Path(sys.argv[0]).name,
        "argv": sys.argv[1:],
        "git_commit": git_hash(),
        "python": sys.version.split()[0],
        "params": params or {},
        "inputs": [{"path": str(p), "md5": _md5(Path(p))} for p in (inputs or [])],
    }
    Path(str(output) + ".provenance.json").write_text(json.dumps(rec, indent=2))


def ensure_dir(p: Path) -> Path:
    Path(p).mkdir(parents=True, exist_ok=True)
    return Path(p)
