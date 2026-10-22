#!/usr/bin/env python3
"""D065. Primary absence rule is the full presence criterion.

Present: alt >= 3 and VAF >= 0.03 (unchanged).
Absent: P(alt >= max(3, ceil(0.03 * depth))) >= 0.80 at max(half the other VAF, 0.03).
The historical rule, P(alt >= 3) only, is written beside it and is not the presence test.

Does not re-read VCFs, BAMs, FACETS, or the dilution/rescue simulators.
Does not recompute cancer-cell-fraction numbers.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import PROJECT_DIR, write_provenance  # noqa: E402

FLOOR = 0.03
MIN_ALT = 3
POWER = 0.80
LO, HI = 0.6, 1.6
N_SIM = 10000


def k_of(depth: float) -> int:
    if not np.isfinite(depth) or depth <= 0:
        return MIN_ALT
    return max(MIN_ALT, int(np.ceil(FLOOR * float(depth) - 1e-9)))


def power_at(depth, expected_vaf, k) -> float:
    if not np.isfinite(depth) or depth < k or not np.isfinite(expected_vaf) or expected_vaf <= 0:
        return 0.0
    p = float(np.clip(expected_vaf, 1e-12, 1 - 1e-12))
    return float(binom.sf(k - 1, int(depth), p))


def assign_states(df: pd.DataFrame, joint: bool) -> tuple[list[str], list[str]]:
    t_state, j_state = [], []
    for r in df.itertuples(index=False):
        pair = {}
        for side, other in (("T", "J"), ("J", "T")):
            if bool(getattr(r, f"{side}_present")):
                pair[side] = "present"
                continue
            other_v = getattr(r, f"{other}_vaf")
            other_v = 0.0 if pd.isna(other_v) else float(other_v)
            expected = max(0.5 * other_v, FLOOR)
            depth = getattr(r, f"{side}_depth")
            k = k_of(depth) if joint else MIN_ALT
            pair[side] = "absent" if power_at(depth, expected, k) >= POWER else "not_evaluable"
        t_state.append(pair["T"])
        j_state.append(pair["J"])
    return t_state, j_state


def spatial_class(t_state: str, j_state: str) -> str:
    if t_state == "present" and j_state == "present":
        return "shared"
    if t_state == "present":
        return "T_private" if j_state == "absent" else "T_unresolved"
    if j_state == "present":
        return "J_private" if t_state == "absent" else "J_unresolved"
    return "neither"


def summarise(df: pd.DataFrame, col: str) -> pd.DataFrame:
    rows = []
    for patient, g in df.groupby("patient", sort=True):
        c = g[col].value_counts()
        shared = int(c.get("shared", 0))
        t_en = int(c.get("T_private", 0))
        j_en = int(c.get("J_private", 0))
        unresolved = int(c.get("T_unresolved", 0) + c.get("J_unresolved", 0) + c.get("neither", 0))
        resolved = shared + t_en + j_en
        rows.append({
            "patient": patient,
            "rule": col,
            "n": len(g),
            "shared": shared,
            "T_enriched": t_en,
            "J_enriched": j_en,
            "unresolved": unresolved,
            "resolved": resolved,
            "shared_of_resolved": shared / resolved if resolved else np.nan,
        })
    return pd.DataFrame(rows)


def ratio_lineage(spatial: str, t_vaf, j_vaf) -> str:
    both = pd.notna(t_vaf) and pd.notna(j_vaf) and t_vaf > 0 and j_vaf > 0
    ratio = float(t_vaf) / float(j_vaf) if both else np.nan
    if spatial == "shared":
        if both and LO <= ratio <= HI:
            return "trunk"
        if both and ratio < LO:
            return "J-lineage"
        return "unassigned"
    if spatial == "J_private":
        return "J-lineage"
    if spatial == "T_private":
        return "T-lineage"
    return "unassigned"


def presence_lineage(t_present, j_present, t_vaf, j_vaf) -> str:
    """Uses presence only. Does not use the power call."""
    both = bool(t_present) and bool(j_present)
    if both:
        return ratio_lineage("shared", t_vaf, j_vaf)
    if bool(t_present) and not bool(j_present):
        return "T-lineage"
    if bool(j_present) and not bool(t_present):
        return "J-lineage"
    return "unassigned"


def poisbinom_tails(ps, obs: int) -> dict:
    dp = np.array([1.0])
    for p in ps:
        p = float(np.clip(p, 0.0, 1.0))
        nxt = np.zeros(len(dp) + 1)
        nxt[:-1] += dp * (1.0 - p)
        nxt[1:] += dp * p
        dp = nxt
    lo = float(dp[: obs + 1].sum())
    hi = float(dp[obs:].sum())
    return {
        "p_lower_leq_observed": lo,
        "p_upper_geq_observed": hi,
        "expected_positive": float(np.sum(ps)),
    }


def site_probabilities(depth, vaf, slope, error) -> np.ndarray:
    mu = np.clip(slope * np.asarray(vaf, float) + error, 0.0, 1.0)
    n = np.asarray(depth, int)
    return 1.0 - np.power(1.0 - mu, n)


def betabinom_tests(depth, vaf, obs_reads: int, obs_sites: int, slope, conc, seed: int) -> dict:
    depth = np.asarray(depth, int)
    vaf = np.asarray(vaf, float)
    if len(depth) == 0:
        return {"p_reads": np.nan, "p_sites": np.nan, "lo": np.nan, "hi": np.nan, "site_expectation": np.nan}
    mu = np.clip(slope * vaf, 1e-8, 1 - 1e-8)
    alpha = mu * conc
    beta = (1.0 - mu) * conc
    site_expectation = float(np.sum(1.0 - betabinom.pmf(0, depth, alpha, beta)))
    rng = np.random.default_rng(seed)
    sim_reads = np.empty(N_SIM)
    sim_sites = np.empty(N_SIM)
    for i in range(N_SIM):
        draw = betabinom.rvs(depth, alpha, beta, random_state=rng)
        sim_reads[i] = draw.sum()
        sim_sites[i] = np.count_nonzero(draw)
    return {
        "p_reads": float(np.mean(sim_reads <= obs_reads)),
        "p_sites": float(np.mean(sim_sites <= obs_sites)),
        "lo": float(np.percentile(sim_reads, 2.5)),
        "hi": float(np.percentile(sim_reads, 97.5)),
        "site_expectation": site_expectation,
    }


def detection_row(name, g, vaf_col, slope, error, conc, seed) -> dict:
    depth = g["N_depth"].to_numpy(float)
    vaf = g[vaf_col].fillna(0).to_numpy(float)
    obs_reads = int(g["N_alt"].fillna(0).sum())
    expected_reads = float(np.sum(depth * slope * vaf))
    positive = int((g["N_alt"].fillna(0) >= 1).sum())
    p_bin = site_probabilities(depth, vaf, slope, 0.0)
    p_err = site_probabilities(depth, vaf, slope, error)
    bin_tail = poisbinom_tails(p_bin, positive)
    err_tail = poisbinom_tails(p_err, positive)
    bb = betabinom_tests(depth, vaf, obs_reads, positive, slope, conc, seed)
    return {
        "group": name,
        "n_sites": len(g),
        "n_positive_sites": positive,
        "N_alt_observed": obs_reads,
        "N_alt_expected_at_slope": expected_reads,
        "expected_positive_binomial": bin_tail["expected_positive"],
        "expected_positive_binomial_plus_error": err_tail["expected_positive"],
        "expected_positive_betabinom": bb["site_expectation"],
        "poisson_binomial_p_lower": bin_tail["p_lower_leq_observed"],
        "poisson_binomial_p_lower_plus_error": err_tail["p_lower_leq_observed"],
        "beta_binomial_p_slope_fixed": bb["p_reads"],
        "beta_binomial_positive_site_p": bb["p_sites"],
        "betabinom_predictive_lo": bb["lo"],
        "betabinom_predictive_hi": bb["hi"],
        "source_vaf": vaf_col,
        "slope": slope,
        "error_rate": error,
        "note": "Slope held fixed. Predictive limits are the 2.5 and 97.5 percentiles of the beta-binomial read sum. Not the slope-bootstrap test.",
    }


def main() -> None:
    root = PROJECT_DIR
    var_path = root / "06_variants" / "somatic_filtered.tsv"
    ccf_path = root / "09_ccf_clonality" / "ccf.tsv"
    prior = pd.read_csv(root / "08_lookback_N" / "background_prior.tsv", sep="\t").iloc[0]
    slope_tab = pd.read_csv(root / "08_lookback_N" / "tumour_in_normal_estimates.tsv", sep="\t")
    slope = float(slope_tab.loc[slope_tab["patient"] == "DH00023", "TiN_shared"].iloc[0])
    error = float(prior["mean_error_rate"])
    conc = float(prior["alpha0"]) + float(prior["beta0"])

    V = pd.read_csv(var_path, sep="\t")
    if len(V) != 380:
        raise SystemExit(f"expected 380 primary rows, found {len(V)}")
    t_old, j_old = assign_states(V, joint=False)
    reproduced = [spatial_class(t, j) for t, j in zip(t_old, j_old)]
    # The reads-only labels are recomputed from depth and VAF. They are not read back
    # from spatial_class, which this script sets to the joint rule.

    t_new, j_new = assign_states(V, joint=True)
    new_class = [spatial_class(t, j) for t, j in zip(t_new, j_new)]
    both = V[["var_id", "patient", "chrom", "pos", "ref", "alt", "gene", "hgvsp_short", "consequence"]].copy()
    both["T_state_reads_only"] = t_old
    both["J_state_reads_only"] = j_old
    both["spatial_class_reads_only"] = reproduced
    both["T_state"] = t_new
    both["J_state"] = j_new
    both["spatial_class"] = new_class
    both["changed"] = both["spatial_class"] != both["spatial_class_reads_only"]

    log = root / "00_admin" / "decision_log.md"
    text = log.read_text()
    if "| D065 |" not in text:
        text += (
            "| D065 | 2026-10-16 | Primary absence uses the full presence rule: "
            "P(alt >= max(3, ceil(0.03*depth))) >= 0.80 at max(half the other VAF, 0.03). "
            "Presence itself is unchanged. The reads-only rule P(alt >= 3) is kept in "
            "spatial_class_both_rules.tsv and is not described as a test of the presence rule. "
            "DH00023 lineages and the normal-read tests are recomputed on the joint-rule classes. "
            "The previous 18/41/46 partition and its p = 0.011 stay in the reads-only files. | "
            "At depth 200 the presence rule needs six alternate reads, while P(alt >= 3) only asks about three. "
            "The follow-up audit reproduced 31 T-enriched and 18 J-enriched under the joint rule. | "
            "Keep P(alt >= 3) as the primary absence call |\n"
        )
    if "| D066 |" not in text:
        text += (
            "| D066 | 2026-10-16 | The nine variants with gnomAD allele frequency above 0.001 stay in the 380. "
            "gnomad_af_review.tsv records normal reads, whether the site recurs, and the class counts "
            "with those nine removed. Removal is a sensitivity, not a new primary set. | "
            "A population frequency is not by itself germline evidence, and adjacent tissue is not blood. "
            "The chr9:64407964 call is in two patients. | Drop the nine from the primary table |\n"
        )
    log.write_text(text)

    out_both = root / "06_variants" / "spatial_class_both_rules.tsv"
    both.to_csv(out_both, sep="\t", index=False)
    counts = pd.concat([
        summarise(both, "spatial_class_reads_only"),
        summarise(both, "spatial_class"),
    ], ignore_index=True)
    out_counts = root / "12_tables" / "supplementary" / "spatial_class_both_rules_counts.tsv"
    counts.to_csv(out_counts, sep="\t", index=False)

    V["T_state"] = t_new
    V["J_state"] = j_new
    V["spatial_class"] = new_class
    V.to_csv(var_path, sep="\t", index=False)

    ccf = pd.read_csv(ccf_path, sep="\t")
    mapped = both.set_index("var_id")
    for col in ("T_state", "J_state", "spatial_class"):
        ccf[col] = ccf["var_id"].map(mapped[col])
    ccf.to_csv(ccf_path, sep="\t", index=False)

    for rel in (
        "12_tables/main/table2_driver_alterations.tsv",
        "12_tables/main/table2_driver_context.tsv",
        "12_tables/main/table_oncoprint_variants.tsv",
    ):
        path = root / rel
        if not path.exists():
            continue
        tab = pd.read_csv(path, sep="\t")
        if "var_id" in tab.columns:
            tab["spatial_class"] = tab["var_id"].map(mapped["spatial_class"]).fillna(tab["spatial_class"])
        elif {"patient", "chrom", "pos", "ref", "alt"}.issubset(tab.columns):
            key = both.set_index(["patient", "chrom", "pos", "ref", "alt"])["spatial_class"]
            idx = pd.MultiIndex.from_frame(tab[["patient", "chrom", "pos", "ref", "alt"]])
            mapped_class = pd.Series(idx.map(key), index=tab.index)
            tab["spatial_class"] = mapped_class.where(mapped_class.notna(), tab["spatial_class"])
        tab.to_csv(path, sep="\t", index=False)

    d = V[V["patient"] == "DH00023"].copy()
    old_d = both[both["patient"] == "DH00023"]
    d = d.merge(old_d[["var_id", "spatial_class_reads_only"]], on="var_id")
    d["lineage_joint"] = [ratio_lineage(s, t, j) for s, t, j in zip(d["spatial_class"], d["T_vaf"], d["J_vaf"])]
    d["lineage_reads_only"] = [
        ratio_lineage(s, t, j) for s, t, j in zip(d["spatial_class_reads_only"], d["T_vaf"], d["J_vaf"])
    ]
    d["lineage_vaf_pattern"] = [
        presence_lineage(tp, jp, t, j)
        for tp, jp, t, j in zip(d["T_present"], d["J_present"], d["T_vaf"], d["J_vaf"])
    ]
    lineage_path = root / "09_ccf_clonality" / "DH00023_lineages.tsv"
    reads_only_lineage = root / "09_ccf_clonality" / "DH00023_lineages_reads_only.tsv"
    if reads_only_lineage.exists():
        old_lineage = pd.read_csv(reads_only_lineage, sep="\t")
    else:
        old_lineage = pd.read_csv(lineage_path, sep="\t")
        old_lineage.to_csv(reads_only_lineage, sep="\t", index=False)
    keep_cols = [c for c in old_lineage.columns if c in d.columns or c in ("record", "region")]
    variants = d.copy()
    variants["record"] = "variant"
    variants["region"] = ""
    variants["spatial_class"] = variants["spatial_class"]
    variants["lineage"] = variants["lineage_joint"]
    variants["vaf_ratio"] = np.where((variants["T_vaf"] > 0) & (variants["J_vaf"] > 0), variants["T_vaf"] / variants["J_vaf"], np.nan)
    for col in old_lineage.columns:
        if col not in variants.columns:
            variants[col] = np.nan
    lineage_out = variants[old_lineage.columns]
    lineage_out.to_csv(root / "09_ccf_clonality" / "DH00023_lineages.tsv", sep="\t", index=False)
    assign_path = root / "09_ccf_clonality" / "DH00023_lineage_assignments.tsv"
    d[["var_id", "gene", "hgvsp_short", "spatial_class_reads_only", "spatial_class",
       "lineage_reads_only", "lineage_joint", "lineage_vaf_pattern",
       "T_vaf", "J_vaf", "N_vaf", "N_alt", "N_depth", "T_present", "J_present"]].to_csv(assign_path, sep="\t", index=False)

    ccf23 = ccf[ccf["patient"] == "DH00023"].merge(
        d[["var_id", "lineage_joint", "lineage_reads_only", "lineage_vaf_pattern"]], on="var_id"
    )
    neutral = (
        (ccf23["T_tcn"].round() == 2) & (ccf23["T_lcn"].round() == 1)
        & (ccf23["J_tcn"].round() == 2) & (ccf23["J_lcn"].round() == 1)
    )

    def fraction_rows(frame, lineage_col, label):
        rows = []
        trunk = frame[frame[lineage_col] == "trunk"]
        for region, vaf_col, ccf_col in (("T", "T_vaf", "T_ccf_raw"), ("J", "J_vaf", "J_ccf_raw")):
            trunk_v = trunk[vaf_col].replace(0, np.nan).median()
            for name, sub in frame.groupby(lineage_col):
                med = sub[vaf_col].median()
                rows.append({
                    "definition": label,
                    "copy_neutral": False,
                    "region": region,
                    "lineage": name,
                    "n": len(sub),
                    "median_vaf": med,
                    "vaf_ratio_to_trunk": (med / trunk_v) if pd.notna(trunk_v) and trunk_v else np.nan,
                    "median_uncapped_ccf": sub[ccf_col].median(),
                })
        neut = frame[neutral.reindex(frame.index).fillna(False)]
        trunk_n = neut[neut[lineage_col] == "trunk"]
        for region, vaf_col, ccf_col in (("T", "T_vaf", "T_ccf_raw"), ("J", "J_vaf", "J_ccf_raw")):
            trunk_v = trunk_n[vaf_col].replace(0, np.nan).median()
            for name, sub in neut.groupby(lineage_col):
                med = sub[vaf_col].median()
                rows.append({
                    "definition": label,
                    "copy_neutral": True,
                    "region": region,
                    "lineage": name,
                    "n": len(sub),
                    "median_vaf": med,
                    "vaf_ratio_to_trunk": (med / trunk_v) if pd.notna(trunk_v) and trunk_v else np.nan,
                    "median_uncapped_ccf": sub[ccf_col].median(),
                })
        return rows

    frac = pd.DataFrame(
        fraction_rows(ccf23, "lineage_joint", "joint_rule")
        + fraction_rows(ccf23, "lineage_reads_only", "reads_only")
        + fraction_rows(ccf23, "lineage_vaf_pattern", "vaf_pattern")
    )
    frac_path = root / "09_ccf_clonality" / "dh00023_lineage_fractions_d065.tsv"
    frac.to_csv(frac_path, sep="\t", index=False)

    groups = {
        "joint_T_enriched": (d["spatial_class"] == "T_private", "T_vaf", 1),
        "joint_J_enriched": (d["spatial_class"] == "J_private", "J_vaf", 2),
        "joint_T_lineage": (d["lineage_joint"] == "T-lineage", "T_vaf", 3),
        "joint_J_lineage": (d["lineage_joint"] == "J-lineage", "J_vaf", 4),
        "joint_trunk": (d["lineage_joint"] == "trunk", "T_vaf", 5),
        "reads_only_T_lineage": (d["lineage_reads_only"] == "T-lineage", "T_vaf", 6),
        "reads_only_J_lineage": (d["lineage_reads_only"] == "J-lineage", "J_vaf", 7),
        "reads_only_J_enriched": (d["spatial_class_reads_only"] == "J_private", "J_vaf", 8),
        "vaf_pattern_T_lineage": (d["lineage_vaf_pattern"] == "T-lineage", "T_vaf", 9),
        "vaf_pattern_J_lineage": (d["lineage_vaf_pattern"] == "J-lineage", "J_vaf", 10),
    }
    det_rows = []
    site_rows = []
    for name, (mask, vaf_col, seed) in groups.items():
        g = d.loc[mask]
        det_rows.append(detection_row(name, g, vaf_col, slope, error, conc, 20261016 + seed))
        for kind, err in (("binomial", 0.0), ("binomial_plus_error", error)):
            ps = site_probabilities(g["N_depth"], g[vaf_col].fillna(0), slope, err)
            obs = int((g["N_alt"].fillna(0) >= 1).sum())
            tails = poisbinom_tails(ps, obs)
            site_rows.append({
                "test": "poisson_binomial",
                "lineage": name,
                "n_sites": len(g),
                "n_sites_with_N_alt": obs,
                "N_alt_reads": int(g["N_alt"].fillna(0).sum()),
                "N_depth_reads": int(g["N_depth"].fillna(0).sum()),
                "source_vaf_median": float(g[vaf_col].median()) if len(g) else np.nan,
                "odds_ratio": np.nan,
                "p_fisher": np.nan,
                "note": (
                    f"{kind}; expected positive {tails['expected_positive']:.4f}; "
                    f"lower-tail P(X<=observed) {tails['p_lower_leq_observed']:.6g}; "
                    f"upper-tail P(X>=observed) {tails['p_upper_geq_observed']:.6g}. "
                    "Exact Poisson-binomial. Not the beta-binomial read test."
                ),
            })
    det = pd.DataFrame(det_rows)
    det_path = root / "08_lookback_N" / "dh00023_detection_joint_rule.tsv"
    det.to_csv(det_path, sep="\t", index=False)

    bars_path = root / "09_ccf_clonality" / "dh00023_vaf_ratio_bars.tsv"
    nobs_path = root / "09_ccf_clonality" / "dh00023_n_observed_expected.tsv"
    if not (root / "09_ccf_clonality" / "dh00023_vaf_ratio_bars_reads_only.tsv").exists():
        pd.read_csv(bars_path, sep="\t").to_csv(
            root / "09_ccf_clonality" / "dh00023_vaf_ratio_bars_reads_only.tsv", sep="\t", index=False
        )
        pd.read_csv(nobs_path, sep="\t").to_csv(
            root / "09_ccf_clonality" / "dh00023_n_observed_expected_reads_only.tsv", sep="\t", index=False
        )
    bar_rows = frac[(frac["definition"] == "joint_rule") & (~frac["copy_neutral"]) & (frac["lineage"] != "unassigned")]
    bar_rows = bar_rows.rename(columns={"n": "n_sites", "vaf_ratio_to_trunk": "vaf_ratio_uncapped", "median_vaf": "median_vaf"})
    bar_rows["vaf_ratio_capped"] = bar_rows["vaf_ratio_uncapped"].clip(upper=1)
    bar_rows[["region", "lineage", "n_sites", "vaf_ratio_uncapped", "vaf_ratio_capped", "median_vaf"]].to_csv(
        bars_path, sep="\t", index=False
    )
    nobs_rows = []
    for lineage, group in (("trunk", "joint_trunk"), ("J-lineage", "joint_J_lineage"), ("T-lineage", "joint_T_lineage")):
        row = det.loc[det["group"] == group].iloc[0]
        nobs_rows.append({
            "lineage": lineage,
            "source_vaf": row["source_vaf"],
            "n_sites": int(row["n_sites"]),
            "observed_alt": row["N_alt_observed"],
            "expected_alt": row["N_alt_expected_at_slope"],
            "poisson_lo": np.nan,
            "poisson_hi": np.nan,
            "betabinom_lo": row["betabinom_predictive_lo"],
            "betabinom_hi": row["betabinom_predictive_hi"],
        })
    pd.DataFrame(nobs_rows).to_csv(nobs_path, sep="\t", index=False)
    site_path = root / "08_lookback_N" / "branch_test_sitelevel.tsv"
    old_site = pd.read_csv(site_path, sep="\t")
    old_site = old_site[old_site["test"] != "poisson_binomial"]
    pd.concat([old_site, pd.DataFrame(site_rows)], ignore_index=True).to_csv(site_path, sep="\t", index=False)

    # gnomAD review. The nine stay in the primary table.
    high = V[V["gnomAD_AF"].fillna(0) > 0.001].copy()
    keys = ["chrom", "pos", "ref", "alt"]
    counts_key = V.groupby(keys)["patient"].nunique()
    review = []
    for r in high.itertuples(index=False):
        chrom = str(r.chrom).replace("chr", "")
        review.append({
            "var_id": r.var_id,
            "patient": r.patient,
            "gene": r.gene if pd.notna(r.gene) else "",
            "consequence": r.consequence,
            "gnomAD_AF": r.gnomAD_AF,
            "N_alt": r.N_alt,
            "N_depth": r.N_depth,
            "N_vaf": r.N_vaf,
            "J_vaf": r.J_vaf,
            "T_vaf": r.T_vaf,
            "spatial_class_joint": r.spatial_class,
            "n_patients_with_this_allele": int(counts_key.loc[(r.chrom, r.pos, r.ref, r.alt)]),
            "recurrent_chr9_64407964": chrom == "9" and int(r.pos) == 64407964,
        })
    drop_ids = set(high["var_id"])
    kept = both[~both["var_id"].isin(drop_ids)]
    sens = summarise(kept, "spatial_class")
    sens["rule"] = "joint_excluding_gnomad_gt_0.001"
    d_ex = d[~d["var_id"].isin(drop_ids) & (d["spatial_class"] == "T_private")]
    sens_note = detection_row("joint_T_enriched_excluding_gnomad", d_ex, "T_vaf", slope, error, conc, 99)
    review_path = root / "12_tables" / "supplementary" / "gnomad_af_review.tsv"
    pd.DataFrame(review).to_csv(review_path, sep="\t", index=False)
    sens_path = root / "12_tables" / "supplementary" / "gnomad_exclusion_sensitivity.tsv"
    sens.to_csv(sens_path, sep="\t", index=False)
    pd.DataFrame([sens_note]).to_csv(root / "08_lookback_N" / "dh00023_t_enriched_excluding_gnomad.tsv", sep="\t", index=False)

    un_path = root / "09_ccf_clonality" / "DH00023_unassigned.tsv"
    un = d[d["lineage_joint"] == "unassigned"][[
        "var_id", "gene", "hgvsp_short", "consequence", "spatial_class", "T_vaf", "J_vaf", "N_vaf", "N_alt", "N_depth"
    ]].copy()
    un["vaf_ratio"] = np.where((un["T_vaf"] > 0) & (un["J_vaf"] > 0), un["T_vaf"] / un["J_vaf"], np.nan)
    un.to_csv(un_path, sep="\t", index=False)

    for path in (out_both, out_counts, assign_path, frac_path, det_path, review_path, sens_path):
        write_provenance(path, inputs=[var_path, ccf_path], params={"decision": "D065"})

    print("COUNTS")
    print(counts.to_string(index=False))
    print("DH00023 lineage counts joint")
    print(d["lineage_joint"].value_counts().to_string())
    print("DH00023 lineage counts reads_only")
    print(d["lineage_reads_only"].value_counts().to_string())
    print("DH00023 lineage counts vaf_pattern")
    print(d["lineage_vaf_pattern"].value_counts().to_string())
    print("of reads-only T-lineage still joint T-enriched",
          int(((d["lineage_reads_only"] == "T-lineage") & (d["spatial_class"] == "T_private")).sum()),
          "of", int((d["lineage_reads_only"] == "T-lineage").sum()))
    print("of reads-only J-enriched still joint J-enriched",
          int(((d["spatial_class_reads_only"] == "J_private") & (d["spatial_class"] == "J_private")).sum()),
          "of", int((d["spatial_class_reads_only"] == "J_private").sum()))
    print("DETECTION")
    print(det[["group", "n_sites", "n_positive_sites", "N_alt_observed", "N_alt_expected_at_slope",
               "expected_positive_binomial", "expected_positive_betabinom",
               "poisson_binomial_p_lower", "beta_binomial_p_slope_fixed"]].to_string(index=False))
    changed = both[both["changed"] & both["gene"].isin(["KMT2C", "TBX3", "NF1", "PIK3CA", "TP53", "BRCA2"])]
    print("DRIVER CLASS CHANGES")
    print(changed[["patient", "gene", "hgvsp_short", "spatial_class_reads_only", "spatial_class"]].to_string(index=False)
          if len(changed) else "none of the named drivers changed class")
    print("gnomAD n", len(review))
    print("wrote", out_counts)


if __name__ == "__main__":
    main()
