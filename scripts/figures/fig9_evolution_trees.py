"""Per-patient clonal-sharing schematic (NOT a phylogeny): germline root -> shared trunk -> T and J branches.
Branch length proportional to number of somatic variants. Labels = cancer-related genes. Right: tumour content
is the purity used for CCF. Under root: tumour-derived DNA in adjacent normal and which branch it matches."""
import sys, textwrap, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, "scripts/figures")
from fig_style import setup, save, WIDTH, CLASS_COLOR, SPECIMEN_COLOR, INK2, GRID, panel_label
setup()
V = pd.read_csv("06_variants/somatic_filtered.tsv", sep="\t")
H = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "12_tables/supplementary/pathway_hits.tsv", sep="\t")
TIN = pd.read_csv("08_lookback_N/tumour_in_normal_estimates.tsv", sep="\t").set_index("patient")
REC = pd.read_csv("07_cnv_purity/purity_reconciled.tsv", sep="\t").set_index("sample")
def lineage(p):
    if p not in TIN.index: return ""
    r = TIN.loc[p]
    if r.get("n_shared", 0) < 5: return "\ntoo few variants to test"
    parts, notes = ["trunk"], []
    for cls, name in (("T_private", "tumour"), ("J_private", "junction")):
        n = r.get(f"{cls}_n", 0)
        if not n: continue
        label = str(r.get(f"{cls}_p_label", ""))
        exp = r.get(f"{cls}_N_alt_expected", np.nan)
        if label == "too few reads to test" or (pd.notna(exp) and exp < 5):
            notes.append(f"{name} branch: too few reads to test")
        else:
            shown = label if label else ""
            notes.append(f"{name} branch {shown}".strip())
    if len(parts) == 1 and not notes:
        notes.append("no enriched branches to test")
    return "; ".join(notes) if notes else "trunk only"
EXTRA = {"HERC5","MRC1","NLRP10","RIPK2","RBBP4","SIN3B","TLE4"}
TRUNC = {"stop_gained","frameshift_variant","start_lost","splice_acceptor_variant","splice_donor_variant"}
V["csq"] = V.consequence.str.split("&").str[0]
key = set(H.gene) | EXTRA
def labels(df):
    d = df[df.gene.isin(key) & df.csq.isin(TRUNC | {"missense_variant"})].drop_duplicates("gene")
    d = d.assign(pri=d.is_driver_gene.astype(int)*3 + d.csq.isin(TRUNC).astype(int)*2).sort_values("pri", ascending=False)
    out = []
    for r in d.itertuples():
        if r.gene == "KMT2C" and pd.notna(getattr(r, "hgvsp_short", None)) and str(r.hgvsp_short):
            out.append(f"KMT2C {r.hgvsp_short}")
        else:
            out.append(r.gene + ("*" if r.csq in TRUNC else ""))
    return out[:5]

pts = sorted(V.patient.unique())
fig, axes = plt.subplots(len(pts), 1, figsize=(WIDTH["double"], 160/25.4), gridspec_kw={"hspace": 0.35})
SCALE = 2.2 / 95
for ax, p in zip(axes, pts):
    g = V[V.patient == p]
    cnt = g.spatial_class.value_counts()
    n_sh, n_t, n_j = cnt.get("shared", 0), cnt.get("T_private", 0), cnt.get("J_private", 0)
    n_other = 0
    if p == "DH00023":
        lin_path = "09_ccf_clonality/DH00023_lineages.tsv"
        lin = pd.read_csv(lin_path, sep="\t")
        variants = lin[lin["record"] == "variant"]
        n_sh = int((variants["lineage"] == "trunk").sum())
        n_j = int((variants["lineage"] == "J-lineage").sum())
        n_t = int((variants["lineage"] == "T-lineage").sum())
        n_other = int((~variants["lineage"].isin(["trunk", "J-lineage", "T-lineage"])).sum())
        g = g.merge(variants[["var_id", "lineage"]], on="var_id", how="left")
        g["spatial_class"] = g["lineage"].map(
            {"trunk": "shared", "J-lineage": "J_private", "T-lineage": "T_private"})
        g = g[g["spatial_class"].notna()]
    n_tu, n_ju = (0, 0) if p == "DH00023" else (cnt.get("T_unresolved", 0), cnt.get("J_unresolved", 0))
    sh = g[g.spatial_class == "shared"]
    def used_purity(specimen):
        sid = f"{p}_{specimen}"
        if sid not in REC.index or pd.isna(REC.loc[sid, "used_purity"]):
            return np.nan
        return float(REC.loc[sid, "used_purity"])
    tc, jc = used_purity("T"), used_purity("J")
    x0, y0 = 2.45, 0.0
    xt = x0 + max(n_sh * SCALE, 0.08)
    ax.plot([x0, xt], [y0, y0], color=CLASS_COLOR["shared"], lw=3, solid_capstyle="butt")
    ax.plot(x0, y0, "o", color=SPECIMEN_COLOR["N"], ms=6, zorder=5)
    for side, n, nu, col, lab, cont, cls in ((+1, n_t, n_tu, SPECIMEN_COLOR["T"], "T", tc, "T_private"),
                                             (-1, n_j, n_ju, SPECIMEN_COLOR["J"], "J", jc, "J_private")):
        xk = xt + 0.1
        xe = xk + max(n * SCALE, 0.05)
        ax.plot([xt, xk, xe], [y0, side*0.55, side*0.55], color=col, lw=2.2, solid_capstyle="butt")
        xu = xe + nu * SCALE
        if nu:
            ax.plot([xe, xu], [side*0.55]*2, color=col, lw=1.4, ls=(0, (2, 1.5)), alpha=0.6)
        ax.plot(xu, side*0.55, "o", color=col, ms=6, zorder=5)
        branch_name = ("T-lineage" if lab == "T" else "J-lineage") if p == "DH00023" else ("T-enriched" if lab == "T" else "J-enriched")
        txt = f"{n} {branch_name}"
        if p != "DH00023" and nu:
            txt += f" (+{nu} unresolved)"
        if np.isfinite(cont): txt += f"  ·  purity {cont*100:.0f}%"
        ax.text(xu + 0.1, side*0.55, txt, va="center", fontsize=6.2, color="#1a1a1a")
        genes = labels(g[g.spatial_class == cls])
        if genes:
            ax.text(xk, side*0.55 + side*0.33, ", ".join(genes), fontsize=5.6, fontstyle="italic", color=col, va="center")
    tg = labels(sh)
    tin = TIN.loc[p, "TiN_shared"] if p in TIN.index else np.nan
    ax.text(0.0, 0.78, p, fontsize=7, fontweight="bold", va="center")
    trunk_name = "Trunk" if p == "DH00023" else "Trunk (shared T+J)"
    trunk_txt = f"{trunk_name}: {n_sh}"
    if p == "DH00023":
        trunk_txt += f"  ·  {n_other} unassigned ({n_sh}+{n_j}+{n_t}+{n_other}=121)"
    ax.text(0.0, 0.48, trunk_txt, fontsize=6.0, va="center", color=CLASS_COLOR["shared"])
    if tg:
        half = (len(tg) + 1) // 2
        ax.text(0.0, 0.22, ", ".join(tg[:half]) + ("," if len(tg) > half else "") + ("\n" + ", ".join(tg[half:]) if len(tg) > half else ""),
                fontsize=5.6, fontstyle="italic", color=CLASS_COLOR["shared"], va="center")
    n_note = textwrap.fill(
        f"N: relative VAF slope {tin:.3f}, not a percent of tumour DNA and not a cell fraction. {lineage(p)}".replace("\n", " "), width=62)
    ax.text(0.0, -0.72, n_note, fontsize=5.2, color=SPECIMEN_COLOR["N"], va="center")
    ax.set_xlim(-0.02, 7.4); ax.set_ylim(-1.0, 1.0); ax.axis("off")
fig.text(0.5, 0.035,
         "Inferred illustration, not a phylogeny and not observed cells. Branch length is a variant count "
         "(DH00023: 18 trunk, 26 J-lineage, 19 T-lineage, 58 unassigned), not time and not a cell fraction. "
         "It does not show which junction cells carry FGFR1, whether alterations share a cell, or which non-cancer cells are present. "
         "Dashed segments are unresolved calls. * truncating. Purity is the value used for CCF.",
         ha="center", fontsize=5.2, color=INK2)
save(fig, "Fig9_clonal_sharing_per_patient_PRELIMINARY",
     note="Inferred illustration only. Not a phylogeny, not observed cells, and not a map of FGFR1 or immune-cell types.")
