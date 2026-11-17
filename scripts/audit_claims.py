#!/usr/bin/env python3
"""Check every number the README quotes against the committed results/ artefacts.

Prose is where numbers rot, and an audit that only searches for a substring
certifies claims it never checked (a "0.98" can be satisfied by any 0.98 in the
file). So every claim here is ANCHORED: it names the phrase that must
immediately precede the value, and the value is recomputed from a CSV in
results/. A claim whose anchor no longer exists is a failure, not a silent pass,
and a missing artefact fails the whole run rather than shrinking it.

What is NOT audited, by design: citations and external facts (Squair et al.
2021, GEO GSE96583, the 38.4 MB download size), commit hashes, and the design
constants the author chose (500 genes, 8 donors, 400 cells, SD 0.3, the 12-gene
ISG panel). Those are inputs, not results.

Exit 0 when every claim holds, 1 otherwise.
    python scripts/audit_claims.py [--strict] [--readme PATH]
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"

REQUIRED = (
    "sim_null.csv", "sim_power.csv", "sim_sample_sweep.csv", "sim_normalisation.csv",
    "sim_nbglm_summary.csv", "false_positives.csv", "real_analysis_summary.csv",
    "de_CD14_Monocytes.csv", "celltype_composition_paired.csv",
    "celltype_composition_proportions.csv", "celltype_stratified_counts.csv",
    "celltype_analysis_summary.csv", "megakaryocyte_summary.csv",
    "megakaryocyte_threshold_sweep.csv", "megakaryocyte_donor_counts.csv",
    "mixed_vs_others_summary.csv", "mixed_null_draws.csv",
    "nbglm_real_summary.csv", "nbglm_real_null_balanced.csv",
)

# the value must start within this many characters after the anchor ends
# (markdown emphasis, a sign, a space); anything further is "moved", a failure
SLACK = 3


def pct(x: float, nd: int = 0) -> str:
    """Percent, rounded half-up as a reader would (307/2,000 = 15.35% -> 15.4%);
    float formatting would round the binary 15.3499... down."""
    d = Decimal(repr(round(100 * float(x), 9))).quantize(Decimal(1).scaleb(-nd), ROUND_HALF_UP)
    return f"{d}%"


def thousands(n: float) -> str:
    return f"{int(round(n)):,}"


def ratio_range(lo: float, hi: float, nd: int) -> list[str]:
    a, b = f"{lo:.{nd}f}", f"{hi:.{nd}f}"
    return [f"{a}–{b}", f"{a}-{b}"]


def build_claims() -> list[tuple[str, list[str], str, str]]:
    """(label, accepted spellings, source, anchor regex)."""
    rd = lambda n: pd.read_csv(RES / n)  # noqa: E731
    C: list[tuple[str, list[str], str, str]] = []

    def add(label, spellings, source, anchor):
        if isinstance(spellings, str):
            spellings = [spellings]
        C.append((label, spellings, source, anchor))

    # ---- §1 simulated null (headline)
    sn = rd("sim_null.csv")
    pc_mean, pb_mean = sn.per_cell.mean(), sn.pseudobulk.mean()
    share = pc_mean / sn.n_genes.iloc[0]
    add("headline per-cell null share", pct(share), "sim_null.csv per_cell/n_genes mean",
        r"a per-cell Wilcoxon test called \*\*")
    add("null caption per-cell count", f"{pc_mean:.1f}", "sim_null.csv per_cell mean",
        r"per-cell Wilcoxon calls ")
    add("null caption per-cell share", pct(share, 1), "sim_null.csv per_cell/n_genes mean",
        r"calls \d+\.\d of 500 genes \(")
    add("null caption pseudobulk count", f"{pb_mean:.1f}", "sim_null.csv pseudobulk mean",
        r"pseudobulk paired t calls ")
    add("null caption replicates", str(sn.rep.nunique()), "sim_null.csv reps",
        r"BH q < 0\.05, (?=\d+ replicates\))")
    add("§1 per-cell count", f"{pc_mean:.1f}", "sim_null.csv per_cell mean",
        r"the figure at the top: per-cell Wilcoxon ")
    add("§1 pseudobulk count", f"{pb_mean:.1f}", "sim_null.csv pseudobulk mean",
        r"the figure at the top: per-cell Wilcoxon \d+\.\d of 500, pseudobulk ")
    add("limits 73%", pct(share), "sim_null.csv", r"\*\*The ")

    # ---- §5 real null
    fp = rd("false_positives.csv")
    fp_share = fp.per_cell / fp.n_genes_sc
    add("headline real per-cell median share", pct(fp_share.median()),
        "false_positives.csv median per_cell/n_genes_sc",
        r"the per-cell test called a median of ")
    add("headline real pseudobulk median", f"{fp.pseudobulk.median():.0f}",
        "false_positives.csv pseudobulk median", r"and pseudobulk a median of ")
    add("headline real splits", str(len(fp)), "false_positives.csv rows",
        r"8 donors split at random ")
    add("limits real median share", pct(fp_share.median()), "false_positives.csv",
        r"the real-data figure is ")
    add("limits real mean share", pct(fp_share.mean()), "false_positives.csv",
        r"the real-data figure is \d+% \(median\) to ")
    add("§5 median count", f"{fp.per_cell.median():.0f}", "false_positives.csv",
        r"Per-cell Wilcoxon: median ")
    add("§5 max count", f"{fp.per_cell.max():.0f}", "false_positives.csv",
        r"false positives per split \(max ")
    add("§5 median share", pct(fp_share.median(), 1), "false_positives.csv",
        r"false positives per split \(max \d+; ")
    add("§5 genes tested", thousands(fp.n_genes_sc.iloc[0]), "false_positives.csv",
        r"per split \(max \d+; [\d.]+% of the ")
    add("§5 mean count", f"{fp.per_cell.mean():.1f}", "false_positives.csv",
        r"genes tested\), mean ")
    add("§5 mean share", pct(fp_share.mean(), 1), "false_positives.csv",
        r"genes tested\), mean \d+\.\d \(")
    add("§5 pseudobulk median", f"{fp.pseudobulk.median():.0f}", "false_positives.csv",
        r"Pseudobulk: median ")
    add("§5 pseudobulk max", f"{fp.pseudobulk.max():.0f}", "false_positives.csv",
        r"Pseudobulk: median \d+, max ")
    add("§5 note median share", pct(fp_share.median(), 1), "false_positives.csv",
        r"over the \d+ splits in `results/false_positives.csv`\. ")
    add("§5 note median ratio", f"{fp.per_cell.median():.0f}/{thousands(fp.n_genes_sc.iloc[0])}",
        "false_positives.csv", r"is the \*median\* share, ")
    add("§5 note mean share", pct(fp_share.mean(), 1), "false_positives.csv",
        r"/2,000; ")
    add("§5 distinct counts", str(fp.per_cell.nunique()), "false_positives.csv",
        r"and only ")
    add("§5 repeated counts", {6: "six"}.get(int((fp.per_cell.value_counts() == 2).sum()),
                                         str(int((fp.per_cell.value_counts() == 2).sum()))),
        "false_positives.csv", r"occur among the 20, ")
    add("§5 distinct splits", str(math.comb(8, 4) // 2), "C(8,4)/2",
        r"eight donors have only ")

    # ---- §2 power at one point
    sp = rd("sim_power.csv")
    g = sp.groupby("method")[["n_called", "tp", "fp", "power", "fdp"]].mean()
    pcw, pbw = g.loc["per_cell"], g.loc["pseudobulk"]
    nrep = str(sp.rep.nunique())
    add("headline pseudobulk TP", f"{pbw.tp:.1f}", "sim_power.csv pseudobulk tp",
        r"with 8 donors it found ")
    add("headline per-cell TP", f"{pcw.tp:.1f}", "sim_power.csv per_cell tp",
        r"true genes against ")
    add("headline per-cell FDP", pct(pcw.fdp), "sim_power.csv per_cell fdp",
        r"whose hit list was ")
    add("limits pseudobulk TP", f"{pbw.tp:.1f}", "sim_power.csv", r"pseudobulk recovers ")
    add("limits per-cell TP", f"{pcw.tp:.1f}", "sim_power.csv",
        r"pseudobulk recovers \d+\.\d of 50 true genes against ")
    add("limits replicates", nrep, "sim_power.csv reps",
        r"400 cells per sample, BH q < 0\.05, (?=\d+ replicates: pseudobulk)")
    add("limits per-cell FDP", pct(pcw.fdp), "sim_power.csv",
        r"not comparable anyway, because ")
    add("glance caption per-cell TP", f"{pcw.tp:.1f}", "sim_power.csv",
        r"finds no more of the true genes \(")
    add("glance caption pseudobulk TP", f"{pbw.tp:.1f}", "sim_power.csv",
        r"finds no more of the true genes \(\d+\.\d vs ")
    add("§2 replicates", nrep, "sim_power.csv", r"mean of (?=\d+ replicates: per-cell)")
    add("§2 per-cell called", f"{pcw.n_called:.1f}", "sim_power.csv",
        r"per-cell Wilcoxon calls (?=\d+\.\d genes \()")
    add("§2 per-cell TP", f"{pcw.tp:.1f}", "sim_power.csv",
        r"per-cell Wilcoxon calls \d+\.\d genes \(")
    add("§2 per-cell FP", f"{pcw.fp:.1f}", "sim_power.csv",
        r"per-cell Wilcoxon calls \d+\.\d genes \(\d+\.\d true, ")
    add("§2 per-cell power", f"{pcw.power:.2f}", "sim_power.csv",
        r"false; power (?=1\.)")
    add("§2 per-cell FDP", f"{pcw.fdp:.2f}", "sim_power.csv",
        r"power \d\.\d\d, observed FDP \*\*")
    add("§2 pseudobulk called", f"{pbw.n_called:.1f}", "sim_power.csv",
        r"and pseudobulk calls ")
    add("§2 pseudobulk TP", f"{pbw.tp:.1f}", "sim_power.csv",
        r"and pseudobulk calls \d+\.\d \(")
    add("§2 pseudobulk FP", f"{pbw.fp:.1f}", "sim_power.csv",
        r"and pseudobulk calls \d+\.\d \(\d+\.\d true, ")
    add("§2 pseudobulk power", f"{pbw.power:.2f}", "sim_power.csv",
        r"\d+\.\d false; power (?=0\.)")
    add("§2 pseudobulk FDP", f"{pbw.fdp:.2f}", "sim_power.csv",
        r"power 0\.\d\d, FDP \*\*")
    add("§2 caption FP", f"{pcw.fp:.0f}", "sim_power.csv", r"per-cell also returns ")
    add("§2 caption FDP", pct(pcw.fdp), "sim_power.csv", r"not real, so ")

    # ---- donor grid (paired t and GLM at 6 / 3 donors)
    ng = rd("sim_nbglm_summary.csv").set_index(["donors", "method"])
    t6 = pct(ng.loc[(6, "paired_t_logcpm"), "power_mean"])
    t3 = pct(ng.loc[(3, "paired_t_logcpm"), "power_mean"])
    nb6 = pct(ng.loc[(6, "nb_glm"), "power_mean"])
    add("headline t at 6 donors", t6, "sim_nbglm_summary.csv 6 paired_t power_mean",
        r"With 6 donors the same paired t found ")
    add("limits t at 6 donors", t6, "sim_nbglm_summary.csv", r"in this repo's own grid \(")
    add("limits t at 3 donors", t3, "sim_nbglm_summary.csv",
        r"in this repo's own grid \(\d+% at 6 donors, ")
    add("glance t at 6 donors", t6, "sim_nbglm_summary.csv",
        r"with 6 donors the paired t finds ")
    add("§2 t at 6 donors", t6, "sim_nbglm_summary.csv",
        r"Pseudobulk's power falls with fewer donors \(")
    add("§2 t at 3 donors", t3, "sim_nbglm_summary.csv",
        r"Pseudobulk's power falls with fewer donors \(\d+% at 6 and ")
    add("§9 GLM at 6 donors", nb6, "sim_nbglm_summary.csv 6 nb_glm power_mean",
        r"it finds ")
    add("§9 t at 6 donors", t6, "sim_nbglm_summary.csv",
        r"of changed genes at 6 donors against ")

    # ---- §3 sample-variation sweep
    sw = rd("sim_sample_sweep.csv").set_index("sample_sd")
    sw_share = sw.per_cell / sn.n_genes.iloc[0]
    add("headline sweep low", pct(sw_share.loc[0.0]), "sim_sample_sweep.csv sd 0",
        r"§3 sweeps it from ")
    add("headline sweep high", pct(sw_share.loc[0.5]), "sim_sample_sweep.csv sd 0.5",
        r"§3 sweeps it from \d+% to ")
    add("limits sweep low", pct(sw_share.loc[0.0]), "sim_sample_sweep.csv",
        r"§3 shows ")
    add("limits sweep high", pct(sw_share.loc[0.5]), "sim_sample_sweep.csv",
        r"§3 shows \d+% at SD 0 and ")
    add("§3 caption sweep row", " / ".join(f"{v:.1f}" for v in sw.per_cell),
        "sim_sample_sweep.csv per_cell", r"0\.00 / 0\.10 / 0\.20 / 0\.30 / 0\.50: ")
    add("§3 caption pseudobulk", f"{sw.pseudobulk.max():.1f}", "sim_sample_sweep.csv",
        r"; pseudobulk ")

    # ---- §4 normalisation
    nm = rd("sim_normalisation.csv").pivot(index="n_donors", columns="normalisation",
                                           values="fdp")
    add("§4 median-of-ratios FDP row", " / ".join(f"{v:.3f}" for v in nm.median_ratio),
        "sim_normalisation.csv median_ratio fdp", r"median-of-ratios ")
    add("§4 CPM FDP row",
        [" / ".join(f"{v:.3f}" for v in nm.cpm[:-1]) + f" / **{nm.cpm.iloc[-1]:.3f}**"],
        "sim_normalisation.csv cpm fdp", r"total-count CPM ")

    # ---- §6 positive control and ranking
    ra = rd("real_analysis_summary.csv").set_index("cell_type")
    mono = ra.loc["CD14+ Monocytes"]
    de14 = rd("de_CD14_Monocytes.csv").set_index("gene")
    for gene, anchor in (("CXCL10", r"\(CXCL10 \+"), ("IFIT1", r", IFIT1 \+"),
                         ("RSAD2", r", RSAD2 \+"), ("ISG15", r", ISG15 \+")):
        add(f"§6 {gene} log2FC", f"{de14.loc[gene, 'log2fc']:.1f}",
            f"de_CD14_Monocytes.csv {gene} log2fc", anchor)
    add("§6 CD14 median ISG rank by q", str(int(mono.median_isg_rank_by_q)),
        "real_analysis_summary.csv CD14 median_isg_rank_by_q",
        r"sits at \*\*rank ")
    add("§6 CD14 ISGs in top 12", "none" if mono.isgs_top12_by_q == 0 else
        str(int(mono.isgs_top12_by_q)), "real_analysis_summary.csv CD14 isgs_top12_by_q",
        r"in that type, and ")
    add("§6 median rank by q across types",
        str(int(ra.median_isg_rank_by_q.median())), "real_analysis_summary.csv",
        r"moves the median interferon-gene rank from ")
    add("§6 median rank by effect across types",
        str(int(ra.median_isg_rank_by_effect.median())), "real_analysis_summary.csv",
        r"interferon-gene rank from \d+ to ")
    n_isg = 12 * len(ra)
    add("§6 ISGs in top 12 by q", f"{int(ra.isgs_top12_by_q.sum())} of {n_isg}",
        "real_analysis_summary.csv", r"in the top 12 from ")
    add("§6 ISGs in top 12 by effect", f"{int(ra.isgs_top12_by_effect.sum())} of {n_isg}",
        "real_analysis_summary.csv", r"in the top 12 from \d+ of \d+ to ")

    def sci(x: float) -> str:
        e = math.floor(math.log10(x))
        sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴"
                            "⁵⁶⁷⁸⁹")
        return f"{x / 10 ** e:.1f} × 10" + str(e).translate(sup)
    add("§6 min q CD14", sci(mono.min_qval), "real_analysis_summary.csv CD14 min_qval",
        r"the smallest q in this cell type is ")
    add("§6 min q all types", sci(ra.min_qval.min()), "real_analysis_summary.csv min",
        r"across all seven types it is ")

    # ---- §7 composition
    cp = rd("celltype_composition_paired.csv").set_index(rd(
        "celltype_composition_paired.csv").columns[0])
    l2 = lambda ct: math.log2(cp.loc[ct, "mean_stim"] / cp.loc[ct, "mean_ctrl"])  # noqa
    sgn = lambda v: f"{v:+.2f}".replace("-", "−")  # noqa: E731
    for ct, anchor in (("CD14+ Monocytes", r"CD14\+ monocytes "),
                       ("Megakaryocytes", r"\), megakaryocytes "),
                       ("CD4 T cells", r"\), CD4 T "), ("NK cells", r"\), NK "),
                       ("B cells", r"\), B ")):
        a, b = cp.loc[ct, "mean_ctrl"], cp.loc[ct, "mean_stim"]
        add(f"§7 composition {ct}", f"{a:.3f} → {b:.3f} ({sgn(l2(ct))})",
            "celltype_composition_paired.csv", r"with the log₂ ratio stim / ctrl: [^.]*?"
            + anchor if ct == "CD14+ Monocytes" else anchor)
    cs = rd("celltype_analysis_summary.csv").iloc[0]
    add("§7 largest shift", f"{cs.max_composition_shift_log2:.2f}",
        "celltype_analysis_summary.csv max_composition_shift_log2",
        r"The largest shift is ")
    add("§7 largest shift fold", f"{2 ** cs.max_composition_shift_log2:.1f}×",
        "celltype_analysis_summary.csv", r"The largest shift is \d\.\d\d log₂ \(")
    mk = cp.loc["Megakaryocytes"]
    add("§7 Mk share ctrl -> stim",
        f"{100 * mk.mean_ctrl:.2f}% → {100 * mk.mean_stim:.2f}%",
        "celltype_composition_paired.csv Megakaryocytes", r"\(megakaryocytes, ")
    prop = rd("celltype_composition_proportions.csv")
    pv = prop.pivot(index="donor", columns="condition", values="CD14+ Monocytes")
    add("§7 CD14 falls in n donors", f"{int((pv.stim < pv.ctrl).sum())} of {len(pv)}",
        "celltype_composition_proportions.csv", r"CD14\+ monocytes also fall, in ")
    add("§7 composition paired t Mk", f"{mk.t:.2f}".replace("-", "−"),
        "celltype_composition_paired.csv Megakaryocytes t", r"paired t = ")
    add("§7 composition p Mk", f"{mk.pval:.2f}", "celltype_composition_paired.csv pval",
        r"paired t = \S+, p = ")
    add("§7 largest shift again", f"{cs.max_composition_shift_log2:.2f}",
        "celltype_analysis_summary.csv", r"The composition shift is small \(largest ")

    sc = rd("celltype_stratified_counts.csv").set_index("cell_type")
    for ct in sc.index:
        r = sc.loc[ct]
        add(f"§7 table {ct}",
            f"{int(r.n_donors)} | {thousands(r.genes_tested)} | "
            f"{thousands(r.sig_within_type_bh)} ({thousands(r.sig_global_bh)})",
            "celltype_stratified_counts.csv", r"\| " + re.escape(ct) + r" \| ")
    for ct, anchor in (("CD14+ Monocytes", r"respond differently \(CD14\+ monocytes "),
                       ("CD8 T cells", r"significant genes, CD8 T cells "),
                       ("CD14+ Monocytes", r"the programme's \*size\* is cell-type "
                                           r"specific: "),
                       ("CD8 T cells", r"genes in CD14\+ monocytes against ")):
        add(f"§7 sig count {ct}", thousands(sc.loc[ct, "sig_within_type_bh"]),
            "celltype_stratified_counts.csv", anchor)
    add("§7 pooled", thousands(cs.pooled_sig), "celltype_analysis_summary.csv pooled_sig",
        r"Pooled over all cells, \*\*")
    add("§7 stratified union", thousands(cs.stratified_union_sig),
        "celltype_analysis_summary.csv stratified_union_sig", r"stratified, \*\*")
    add("§7 pooled only", thousands(cs.pooled_only_sig),
        "celltype_analysis_summary.csv pooled_only_sig", r"significant somewhere\. \*\*")
    add("§7 stratified only", thousands(cs.stratified_only_sig),
        "celltype_analysis_summary.csv stratified_only_sig",
        r"pooled and in no cell type and \*\*")
    add("§7 marker share", pct(cs.pooled_sig_marker_share, 2),
        "celltype_analysis_summary.csv pooled_sig_marker_share", r"Only ")
    add("§7 marker share pooled-only n", thousands(cs.pooled_only_sig),
        "celltype_analysis_summary.csv", r"Only [\d.]+% of the ")
    add("§7 marker genes", str(round(cs.pooled_sig_marker_share * cs.pooled_only_sig)),
        "celltype_analysis_summary.csv share x pooled_only", r"pooled-only genes \(")
    add("§7 pooled describes neither", thousands(cs.pooled_sig),
        "celltype_analysis_summary.csv", r"A single pooled number \(")
    add("§5 Kang cells", thousands(cs.n_cells), "celltype_analysis_summary.csv n_cells",
        r"stimulated aliquot \(")
    add("§7 types tested", str(int(cs.n_cell_types_tested)),
        "celltype_analysis_summary.csv n_cell_types_tested", r"8 cell types are covered: ")

    # ---- megakaryocytes
    ms = rd("megakaryocyte_summary.csv").iloc[0]
    dc = rd("megakaryocyte_donor_counts.csv").set_index("donor")
    add("Mk caption cells", str(int(ms.n_cells)), "megakaryocyte_summary.csv n_cells",
        r"megakaryocyte_donors\.png\)\n\n\*")
    add("Mk caption donors", str(int(ms.donors_with_any_cell)),
        "megakaryocyte_summary.csv donors_with_any_cell", r"\d+ cells across ")
    add("Mk largest donor cells", str(int(dc.loc[ms.largest_donor, "total"])),
        "megakaryocyte_donor_counts.csv", r"`" + re.escape(ms.largest_donor) + r"` supplies ")
    add("Mk largest donor share", pct(ms.largest_donor_share),
        "megakaryocyte_summary.csv largest_donor_share", r"supplies \d+ \(")
    tsw = rd("megakaryocyte_threshold_sweep.csv").set_index("min_cells_per_condition")
    one_donor = sorted(k for k, v in tsw.donors_qualifying.items() if v == 1)
    add("Mk sweep >=25..>=3 leave 1", f"≥{max(one_donor)} through ≥{min(one_donor)}"
        f" cells leave 1 donor", "megakaryocyte_threshold_sweep.csv", r"\n(?=≥)")
    add("Mk sweep >=2", str(int(tsw.loc[2, "donors_qualifying"])),
        "megakaryocyte_threshold_sweep.csv", r"≥2 leaves ")
    add("Mk sweep >=1", str(int(tsw.loc[1, "donors_qualifying"])),
        "megakaryocyte_threshold_sweep.csv", r"≥1 leaves ")
    add("Mk profile range", f"{int(ms.cells_per_profile_min)} to {int(ms.cells_per_profile_max)}",
        "megakaryocyte_summary.csv cells_per_profile_min/max", r"pairs profiles of ")
    add("Mk share again", pct(ms.largest_donor_share), "megakaryocyte_summary.csv",
        r"replicates, with ")
    add("Mk genes called", f"{int(ms.genes_significant)} of {int(ms.genes_tested)}",
        "megakaryocyte_summary.csv", r"and it calls ")
    add("Mk genes surviving", str(int(ms.genes_tested)), "megakaryocyte_summary.csv",
        r"At the loosest threshold ")
    add("Mk paired SD", f"{ms.median_paired_sd_log2:.2f}", "megakaryocyte_summary.csv",
        r"median paired log₂ SD is ")
    add("Mk MDE exact", f"{ms.mde_log2_exact_noncentral_t:.2f}", "megakaryocyte_summary.csv",
        r"at 80% power is \*\*")
    add("Mk MDE central", f"{ms.mde_log2_at_80pct_power:.2f}", "megakaryocyte_summary.csv",
        r"the figure to use; ")
    add("Mk MDE fold", f"{ms.mde_fold:.1f}×", "megakaryocyte_summary.csv mde_fold",
        r"the figure to use; \d\.\d\d log₂, ")
    add("Mk MDE Bonferroni", f"{ms.mde_log2_exact_bonferroni:.2f}",
        "megakaryocyte_summary.csv", r"central-t approximation; \*\*")
    add("Mk Bonferroni genes", str(int(ms.genes_tested)), "megakaryocyte_summary.csv",
        r"with Bonferroni over ")
    add("Mk largest observed", f"{ms.largest_observed_abs_log2fc:.2f}",
        "megakaryocyte_summary.csv largest_observed_abs_log2fc",
        r"The largest observed \|log₂FC\| is ")
    add("Mk earlier sentence", f"{ms.mde_log2_at_80pct_power:.2f}",
        "megakaryocyte_summary.csv", r"An earlier sentence called ")

    # ---- §8 mixed model
    mx = rd("mixed_vs_others_summary.csv")
    cost = mx.mixed_s / mx.pseudobulk_s
    add("§8 heading cost", f"{cost.min():.0f}–{cost.max():.0f}×",
        "mixed_vs_others_summary.csv mixed_s/pseudobulk_s", r"the same test at ")
    calls = mx.mixed_sig / mx.pseudobulk_sig
    add("§8 call ratio", [s + "×" for s in ratio_range(calls.min(), calls.max(), 2)],
        "mixed_vs_others_summary.csv mixed_sig/pseudobulk_sig",
        r"wherever the donor variance is estimated, ")
    bshare = mx.mixed_only_calls_boundary.sum() / mx.mixed_only_calls.sum()
    add("§8 boundary share", pct(bshare), "mixed_vs_others_summary.csv boundary/mixed_only",
        r"the calls \(")
    add("§8 caption boundary share", pct(bshare), "mixed_vs_others_summary.csv",
        r"and carry ")
    add("§8 boundary stat max", f"{mx.mixed_stat_over_t_max_boundary.max():.2f}×",
        "mixed_vs_others_summary.csv mixed_stat_over_t_max_boundary", r"run up to ")
    nd = rd("mixed_null_draws.csv")
    bal = nd[nd.balance == "balanced"]
    add("§8 null draws equal", f"{int((bal.fp_pseudobulk == bal.fp_mixed).sum())} of {len(bal)}",
        "mixed_null_draws.csv balanced fp_pseudobulk == fp_mixed",
        r"the same false-positive count as pseudobulk in ")
    words = {5: "five", 6: "six", 7: "seven", 8: "eight"}
    add("§8 types", words.get(len(mx), str(len(mx))), "mixed_vs_others_summary.csv rows",
        r"On all ")

    # ---- §9 NB GLM
    nb = rd("nbglm_real_summary.csv")
    more = nb.n_nb / nb.n_pb - 1
    add("§9 more calls", f"{100 * more.min():.0f}–{100 * more.max():.0f}%",
        "nbglm_real_summary.csv n_nb/n_pb - 1", r"the negative-binomial GLM calls ")
    rt = nb.t_nb_s / nb.t_pb_s
    add("§9 runtime ratio", f"{rt.min():.1f}–{rt.max():.0f}×",
        "nbglm_real_summary.csv t_nb_s/t_pb_s", r"ranks the interferon genes higher, at ")
    nbn = rd("nbglm_real_null_balanced.csv").set_index(["cell_type", "method"])
    bglm, mglm = nbn.loc[("B cells", "nb_glm")], nbn.loc[("CD14+ Monocytes", "nb_glm")]
    add("§9 B draws with FP", f"{int(bglm.n_draws_with_any_fp)} of {int(bglm.draws)}",
        "nbglm_real_null_balanced.csv B nb_glm", r"false positives on B cells \(")
    add("§9 CD14 draws with FP", f"{int(mglm.n_draws_with_any_fp)} of {int(mglm.draws)}",
        "nbglm_real_null_balanced.csv CD14 nb_glm", r"and CD14\+ monocytes \(")
    add("§9 caption B draws", str(int(bglm.draws)), "nbglm_real_null_balanced.csv",
        r"every one of the ")
    add("§9 caption B mean FP", f"{bglm.mean_fp_small_effect:.2f}",
        "nbglm_real_null_balanced.csv B mean_fp_small_effect",
        r"at least one false positive \(")
    add("§9 caption CD14 draws", f"{int(mglm.n_draws_with_any_fp)} of {int(mglm.draws)}",
        "nbglm_real_null_balanced.csv", r"on CD14\+ monocytes ")
    add("§9 caption CD14 mean FP", f"{mglm.mean_fp_small_effect:.2f}",
        "nbglm_real_null_balanced.csv CD14 mean_fp_small_effect", r"do \(")
    pbn = nbn.xs("pseudobulk", level="method").loc[["B cells", "CD14+ Monocytes"]]
    add("§9 pseudobulk none on either",
        "none" if pbn.n_draws_with_any_fp.sum() == 0 else "SOME",
        "nbglm_real_null_balanced.csv pseudobulk B/CD14", r"pseudobulk gives ")
    return C


def build_doc_claims() -> dict[str, list[tuple[str, list[str], str, str]]]:
    """The two call tables moved verbatim from README §8/§9 into docs/.

    The rest of docs/ (prose, the null tables) is not audited; see the README's
    Reproducing note.
    """
    rd = lambda n: pd.read_csv(RES / n)  # noqa: E731
    out: dict[str, list] = {"docs/mixed_model.md": [], "docs/nb_glm.md": []}
    mx = rd("mixed_vs_others_summary.csv")
    for _, r in mx.iterrows():
        out["docs/mixed_model.md"].append((
            f"mixed table {r.cell_type}",
            [f"{int(r.n_paired_donors)} | {thousands(r.percell_sig)} | "
             f"{thousands(r.pseudobulk_sig)} | {thousands(r.mixed_sig)} | "
             f"{r.mixed_sig / r.pseudobulk_sig:.2f}\u00d7 | {r.mixed_pct_boundary:.1f}%"],
            "mixed_vs_others_summary.csv", r"\| " + re.escape(r.cell_type) + r" \| "))
    out["docs/mixed_model.md"].append((
        "mixed-only calls", [f"{int(mx.mixed_only_calls.sum())} mixed-only calls, "
                             f"{int(mx.mixed_only_calls_boundary.sum())} on boundary"],
        "mixed_vs_others_summary.csv sums", r"80% of the extras on boundary fits:\*\* "))
    nb = rd("nbglm_real_summary.csv")
    for _, r in nb.iterrows():
        out["docs/nb_glm.md"].append((
            f"nb table {r.cell_type}",
            [f"{int(r.donors)} | {thousands(r.n_tested_pb)} | {thousands(r.n_pb)} | "
             f"{thousands(r.n_nb)} | {thousands(r.n_both)} | {r.jaccard:.3f}"],
            "nbglm_real_summary.csv", r"(?i)\| " + re.escape(r.cell_type) + r" \| "))
    return out


def check(readme: str, claims) -> list[tuple]:
    bad = []
    for label, spellings, source, anchor in claims:
        if not anchor:
            if not any(s in readme for s in spellings):
                bad.append((label, spellings, source, "no anchor and not found"))
            continue
        hits = list(re.finditer(anchor, readme))
        if not hits:
            bad.append((label, spellings, source, f"anchor {anchor!r} matched nowhere"))
            continue
        ok = False
        for m in hits:
            after = readme[m.end():m.end() + SLACK + max(map(len, spellings))]
            lead = len(after) - len(after.lstrip("* "))
            if lead <= SLACK and any(after[lead:].startswith(s) for s in spellings):
                ok = True
                break
        if not ok:
            got = readme[hits[0].end():hits[0].end() + 24].replace("\n", " ")
            bad.append((label, spellings, source, f"anchor matched; text there is {got!r}"))
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--readme", default=str(ROOT / "README.md"))
    ap.add_argument("--strict", action="store_true",
                    help="also fail on any claim without an anchor (none exist today; "
                         "this guards a future addition that lacks one)")
    args = ap.parse_args()

    missing = [n for n in REQUIRED if not (RES / n).exists()]
    if missing:
        print(f"MISSING ARTEFACT ({len(missing)}): {', '.join(missing)}")
        print("The claims that depend on it cannot be checked, which is a failure.")
        return 1

    readme = Path(args.readme).read_text(encoding="utf-8")
    claims = build_claims()
    bad = check(readme, claims)
    for doc, doc_claims in build_doc_claims().items():
        text = (ROOT / doc).read_text(encoding="utf-8")
        bad += [(f"{doc}: {b[0]}",) + b[1:] for b in check(text, doc_claims)]
        claims += doc_claims
    unanchored = sum(1 for c in claims if not c[3])
    total = len(claims)
    print(f"{total - len(bad)}/{total} claims hold in place"
          + (f"; {unanchored} unanchored" if unanchored else "; all anchored to a phrase"))
    if bad:
        print("\nWRONG, MOVED OR MISSING:")
        for label, spellings, source, why in bad:
            print(f"  {label:38s} expected {spellings} ({source}) :: {why}")
        print("\nEither the README is stale or the artefacts changed. Neither ships.")
        return 1
    if args.strict and unanchored:
        print(f"\n{unanchored} claim(s) are not bound to a phrase.")
        return 1
    print("every audited README number matches its committed artefact, in place")
    return 0


if __name__ == "__main__":
    sys.exit(main())
