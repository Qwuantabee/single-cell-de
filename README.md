# Are cells replicates? Differential expression tested at the donor level

[![tests](https://github.com/Qwuantabee/single-cell-de/actions/workflows/tests.yml/badge.svg)](https://github.com/Qwuantabee/single-cell-de/actions/workflows/tests.yml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Qwuantabee/single-cell-de/blob/main/LICENSE)
[![python](https://img.shields.io/badge/python-3.10--3.12-blue.svg)](https://github.com/Qwuantabee/single-cell-de/blob/main/pyproject.toml)

**No.** Counting every cell as a replicate, a per-cell Wilcoxon test called **73% of all genes significant** (BH q < 0.05) on simulated data where nothing changed; pseudobulk (one profile per donor and condition, paired t-test) called none. On real control cells (Kang 2018, CD14+ monocytes, 8 donors split at random 20 times) the per-cell test called a median of 15% of genes and pseudobulk a median of 0.

What pseudobulk costs in power was measured at one point: at a 2.72-fold effect with 8 donors it found 48.8 of 50 true genes against 50.0 for the per-cell test, whose hit list was 87% false. With 6 donors the same paired t found 68% of true genes in the simulated donor grid ([docs/nb_glm.md](docs/nb_glm.md)). The 73% depends on the simulator's assumed sample-to-sample variation (§3 sweeps it from 0% to 81%).

![null false positives](figures/null_false_positives.png)

*Simulated null, 500 genes, 8 donors: per-cell Wilcoxon calls 364.7 of 500 genes (72.9%); pseudobulk paired t calls 0.0 (BH q < 0.05, 10 replicates).*

## The argument

Cells from one donor are not independent observations. Treating 6,400 cells from 8 donors as 6,400 samples inflates the effective sample size by three orders of magnitude, and the p-values stop meaning what they say. The replicate is the donor, because the donor is what was randomised. This is established (Squair et al. 2021, Nat Commun 12:5692) but still the default in most tutorials, so it is worth demonstrating, not citing. The sections below measure it on simulated and real data, then hold two alternatives (a mixed model, a negative-binomial GLM) to the same standard.

## Limits, stated plainly

- **The 73% is a simulator setting, not a constant of nature.** It is measured at one assumed sample-level SD (0.3 on the natural-log scale, a simulator default chosen by the author; no fit to data is in this repo). §3 shows 0% at SD 0 and 81% at SD 0.5; the real-data figure is 15% (median) to 20% (mean) on one cell type of one dataset.
- **"No power lost" is measured at one point.** One effect size (natural-log effect 1.0: 2.72-fold, 1.44 log₂), 8 donors, 400 cells per sample, BH q < 0.05, 10 replicates: pseudobulk recovers 48.8 of 50 true genes against 50.0. Smaller effects were not simulated against the per-cell test. With fewer donors the pseudobulk paired t does lose power in this repo's own grid (68% at 6 donors, 0% at 3; [docs/nb_glm.md](docs/nb_glm.md)), and per-cell "power" is not comparable anyway, because 87% of its calls are false.
- **One real dataset** (Kang 2018, 8 donors, one stimulus). The real-data nulls can show anti-conservatism, not calibration (§8, §9), and several have 15–70 draws.
- **Real data has no ground truth**: power and false discovery proportion are simulated only; the simulator draws negative-binomial counts, which is the NB GLM's own model.
- The donor counts, gene filters, ISG panel (12 genes) and the |log₂FC| < 0.5 "small-effect" null set are the author's choices.

## Corrections to this repo's own earlier claims

Each is described in the section where it applies.

| what was wrong | where | commit |
|---|---|---|
| Composition table had ctrl and stim swapped, so every sign was reversed | §7 | `6d5fb63` |
| Null reference was 25 of 500 (uncorrected p) where BH expects about 0.05 of 500 | §1, §8 | `b72cde6`, `6f3e0fd` |
| Mixed model reported 12–272× null inflation: a normal reference at 4–8 pairs, and unpaired profiles in the mixed arm | §8 | `7f1e46d`, `96f8026` |
| GLM and pseudobulk log₂FC "agreed to three decimals" while on different normalisations | §9 | `8f09882` |
| Nulls resampled plans that do not exist as distinct draws; GLM "calibrated" on one type and ten draws | §8, §9 | `eb13b69`, `4b90d2a` |
| Simulated effect labelled 2-fold / log₂ 1.0; it is 1.0 on the natural-log scale (2.72-fold) | §2, docs/nb_glm.md | `5a6d99e`, `2ae8f56` |
| Three result files no script writes or reads (`mixed_cd4.csv`, the CD4-only mixed run on the withdrawn normal reference, with p = 0 rows; `mixed_percell_cd4.csv`; `sim_donor_sweep.csv`) sat beside the corrected ones | `results/` | withdrawn from the tree in this revision; still in history at `8431d1c`, `d8fcc0d` |

The simulator fix described under "How the simulation got fixed" and the first-version normalisation numbers in §4 predate the first commit, so they are narrated here and cannot be checked from the history.

## At a glance

![at a glance](figures/at_a_glance.png)

*A: per-cell testing manufactures false positives on null data, simulated and real; pseudobulk does not. B: at the one effect size and donor count simulated (2.72-fold, 8 donors, BH q < 0.05) per-cell testing finds no more of the true genes (50.0 vs 48.8 of 50); with 6 donors the paired t finds 68%, so this is not a general power result (§2, Limits). C: the donor-intercept mixed model is the pseudobulk test at far higher cost (§8). D: the NB GLM calls more genes but is anti-conservative on B cells and CD14+ monocytes (§9). C and D correct earlier claims of this repo.*

## 1. False positives when nothing is happening

500 genes, 8 donors, 400 cells per donor per condition, zero genes truly changed: every call is a false positive by construction (no gene has a condition effect; each sample gets an independent shift with SD 0.3 on the natural-log scale, the simulator default, swept in §3). The result is the figure at the top: per-cell Wilcoxon 364.7 of 500, pseudobulk 0.0.

The dashed line is the *uncorrected* p < 0.05 expectation (5% × 500 = 25). Under the Benjamini-Hochberg correction both methods use at q = 0.05, a calibrated test on all-null data calls a mean of about **0.05 of 500**, not 25: BH bounds the false-discovery *proportion* among rejections, so a correct procedure returns essentially nothing. That distinction matters for how the mixed-model null in §8 is read.

## 2. It is not buying power (at the effect size tested)

The obvious defence of per-cell testing is sensitivity. With 50 of 500 genes given a real effect of 1.0 on the natural-log scale (2.72-fold, 1.44 log₂; an earlier version of this sentence called it 2-fold, but `simulate.py` applies the effect inside `exp`), 8 donors, 400 cells per sample, BH q < 0.05, mean of 10 replicates: per-cell Wilcoxon calls 376.5 genes (50.0 true, 326.5 false; power 1.00, observed FDP **0.87**) and pseudobulk calls 51.3 (48.8 true, 2.5 false; power 0.98, FDP **0.05**).

![power vs fdp](figures/power_vs_fdp.png)

*Both recover the real genes; per-cell also returns 326 that are not real, so 87% of its hit list is wrong. It is not finding more, it is calling more.*

This is one effect size at one donor count, not a power curve. Pseudobulk's power falls with fewer donors (68% at 6 and 0% at 3 for the same paired t in the simulated donor grid, [docs/nb_glm.md](docs/nb_glm.md)), and no smaller effect was simulated against the per-cell test, so "loses no power" holds here and is not shown elsewhere. A per-cell test that calls three quarters of all genes will always contain the true ones; that is not power in any useful sense.

## 3. Why — and the test that proves it

The per-cell test ignores *sample-level* variation: everything that differs between two aliquots of one donor's cells (processing day, capture efficiency, that individual's response). It sees thousands of cells agree and reads one sample's idiosyncrasy as strong evidence. If that is the cause, the test should behave when sample-level variation is zero and degrade as it grows.

![sample variation sweep](figures/sample_variation_sweep.png)

*Per-cell false positives (of 500 null genes) at sample-level SD 0.00 / 0.10 / 0.20 / 0.30 / 0.50: 0.0 / 138.2 / 306.0 / 365.4 / 404.2; pseudobulk 0.0 throughout (mean of 5 replicates per point, natural-log SD).*

It is exactly calibrated at zero and collapses at once. The per-cell test is not broken in general: it validly tests "did these cells come from the same distribution", which nobody asks. The question is whether the *condition* did something, and that has 8 replicates, not 6,400. Per-cell testing is fine for comparing cell types within one sample, where there is no sample-level confound.

## 4. Normalisation, and a bias that grows with study size

Total-count normalisation assumes every sample holds the same amount of RNA. When 50 genes go up strongly they take up more of the library, so every other gene's share falls though nothing happened to it. In the first version of the simulation this shifted every null gene by −0.07 log₂ CPM, consistently enough across donors that the paired test called 163 of them significant (a 57% false discovery rate from an otherwise correct method). Median-of-ratios size factors (as in DESeq2) fixed it, and with realistic sample-level noise the two normalisations looked identical at 8 donors.

![normalisation](figures/normalisation.png)

*False discovery proportion at 8 / 20 / 40 donors: median-of-ratios 0.043 / 0.059 / 0.059, total-count CPM 0.049 / 0.098 / **0.192** (mean of 6 replicates per cell).*

Composition bias is a *bias*: the same in every sample, so more donors do not average it out while sampling noise shrinks. A choice invisible at 8 donors dominates the error at 40, against the intuition that more data makes methodological choices matter less.

## 5. The same experiment on real data

Kang et al. 2018: PBMCs from 8 lupus patients, each split into a control and an interferon-β-stimulated aliquot (24,673 cells, GEO GSE96583). Taking **control cells only**, where no donor differs from any other by treatment, and randomly splitting the 8 donors into two fake "conditions" 20 times, in CD14+ monocytes:

![real null](figures/real_null.png)

*Per-cell Wilcoxon: median 307 false positives per split (max 988; 15.4% of the 2,000 genes tested), mean 395.9 (19.8%). Pseudobulk: median 0, max 1. The simulation was not exaggerating.*

(Share = false positives / genes tested per split, over the 20 splits in `results/false_positives.csv`. 15.4% is the *median* share, 307/2,000; 19.8% is the *mean* share. An earlier version of this row mixed the two statistics. The 20 splits are random draws, not 20 distinct splits: eight donors have only 35 distinct 4-vs-4 splits, and only 14 distinct per-cell counts occur among the 20, six of them twice. This is also an unpaired 4-vs-4 comparison of whole donors, a different design from the paired one used everywhere else.)

## 6. The positive control, and a ranking problem

A method that calls nothing is trivially free of false positives, so check that pseudobulk still finds the real response. Interferon-β is a well-characterised programme, so the answer is known in advance. In CD14+ monocytes the effects are enormous (CXCL10 +9.8 log₂, IFIT1 +8.6, RSAD2 +8.3, ISG15 +7.7) and all are significant. But ranked by q-value, the median known interferon gene sits at **rank 441** in that type, and none is in the top 12.

![isg ranking](figures/isg_ranking.png)

*Ranking the same significant genes by effect size moves the median interferon-gene rank from 50 to 10 across the seven cell types, and the interferon genes in the top 12 from 9 of 84 to 51 of 84.*

The cause is the degrees of freedom. A paired t-test on 8 donors has 7 df, which floors the p-value (the smallest q in this cell type is 1.5 × 10⁻⁶; across all seven types it is 6.6 × 10⁻⁷, on CD4 T cells). Hundreds of genes pile against the floor and their order is decided by noise, so a modest but very consistent change outranks a 400-fold one that wobbles between donors. Use the p-value as a filter and the effect size as the ranking.

Three "rank" definitions appear in this README and are not interchangeable; compare only inside one.

<details><summary>The three rank definitions</summary>

(a) The figure (`results/real_analysis_summary.csv`, `real_analysis.py`): rank *within the genes already significant* at BH q < 0.05, by q or by |log₂FC|, median over the 12 ISGs truncated to an integer. (b) §7 (`results/celltype_isg_recovery.csv`, `celltype_analysis.py`): the same definition from a different script; its medians differ by 1–2 for three types (CD4 73 vs 74, CD8 31 vs 32, FCGR3A+ monocytes 172 vs 170), and the cause was not isolated (re-sorting with different algorithms did not reproduce it, so it is not tie-breaking alone). (c) §9 (`results/nbglm_real_summary.csv`): rank among *all tested genes* by q, median not truncated (hence 73.5).

</details>

## 7. Cell type first, or the composition shift is the result

If the *mix* of cell types differs between conditions, a pooled comparison compares two cell mixtures, not two expression states: every marker of whichever type became more common comes out "up", with valid p-values, for a question nobody asked. On Kang the composition moves only slightly. Mean proportion ctrl → stim, with the log₂ ratio stim / ctrl: CD14+ monocytes 0.254 → 0.238 (−0.10), megakaryocytes 0.015 → 0.012 (−0.30), CD4 T 0.440 → 0.448 (+0.03), NK 0.067 → 0.068 (+0.01), B 0.097 → 0.102 (+0.07).

The largest shift is 0.30 log₂ (1.2×) and a *decrease* (megakaryocytes, 1.51% → 1.22% of cells); CD14+ monocytes also fall, in 7 of 8 donors. Every sign is the reverse of an earlier version of this table, whose `mean_ctrl` and `mean_stim` columns had been read in the wrong order (commit `6d5fb63`). Interferon-β does not rearrange PBMC composition much, so this dataset is *not* the case where pooling fails for that reason; a disease contrast, where composition effects usually bite, would be.

<details><summary>Sign convention</summary>

The log₂ ratio above is stim / ctrl, the opposite of `results/celltype_composition_shift.csv` (`log2_ratio` = ctrl / stim; `mean_prop_a` is ctrl, `mean_prop_b` is stim; `tests/test_celltype_pseudobulk.py` pins this). Proportions come from `celltype_composition_paired.csv`, the source of record (no log₂ column).

</details>

Stratifying changes the answer anyway, because cell types respond differently (CD14+ monocytes 2,952 significant genes, CD8 T cells 88):

<details><summary>Per-type counts (table)</summary>

| cell type | donors | genes tested | significant, BH within type (global) |
|---|---|---|---|
| CD14+ Monocytes | 8 | 5,364 | 2,952 (2,578) |
| CD4 T cells | 8 | 7,048 | 986 (1,191) |
| FCGR3A+ Monocytes | 5 | 2,673 | 746 (717) |
| Dendritic cells | 5 | 2,078 | 565 (551) |
| B cells | 8 | 3,804 | 563 (661) |
| NK cells | 7 | 2,436 | 401 (459) |
| CD8 T cells | 5 | 2,081 | 88 (187) |
| Megakaryocytes | 1 | — | not testable |

</details>

**All 8 cell types are covered: 7 by an expression test, 8 by the composition test.** Megakaryocytes get a composition result (the largest shift above; paired t = −0.91, p = 0.39 on all 8 donors) and no expression result, and that exclusion is not a threshold that could be relaxed (an earlier version implied it was, by quoting only the ≥25-cell rule).

![megakaryocyte donors](figures/megakaryocyte_donors.png)

*132 cells across 8 donors, but `patient_1039` supplies 111 (84%); no cell-count threshold leaves a usable design.*

≥25 through ≥3 cells leave 1 donor, ≥2 leaves 2 (a paired t on 1 df) and ≥1 leaves 4; the loosest pairs profiles of 1 to 58 cells as equally informative replicates, with 84% of the type from one donor, and it calls 0 of 85 genes. That is a power statement, not a biological one: megakaryocytes may respond, this dataset cannot say.

<details><summary>Megakaryocyte power numbers</summary>

At the loosest threshold 85 genes survive filtering, the median paired log₂ SD is 0.80, and the minimum detectable effect at 80% power is **1.71 log₂** (exact noncentral-t, the figure to use; 1.67 log₂, 3.2×, by the central-t approximation; **7.75 log₂** with Bonferroni over 85 genes). The largest observed |log₂FC| is 3.49 (ISG15): above the per-gene minimum and well below the Bonferroni one, so one gene responding that much is not ruled out but nothing could survive correction. (An earlier sentence called 1.67 "larger than anything in the interferon programme here"; the observed ISG15 value contradicts that.) Files: `results/megakaryocyte_threshold_sweep.csv`, `megakaryocyte_summary.csv` (`expression_test_reliable = False`), `megakaryocyte_donor_counts.csv`.

</details>

Pooled over all cells, **3,211** genes are significant; stratified, **3,512** distinct genes are significant somewhere. **763** are significant pooled and in no cell type and **1,064** in at least one type and not pooled, so the pooled answer is a *different* answer, not a superset or subset. This analysis does **not** establish why. The composition shift is small (largest 0.30 log₂), so composition is unlikely to be the whole explanation; different responses being averaged, and different gene filters, are other candidates that were not separated. Only 0.52% of the 763 pooled-only genes (4 genes) are among the top-30 markers of any cell type (`pooled_sig_marker_share` in `results/celltype_analysis_summary.csv`), where a composition artefact would put marker genes.

All 12 known ISGs are significant in all 7 testable types, but the programme's *size* is cell-type specific: 2,952 genes in CD14+ monocytes against 88 in CD8 T cells. A single pooled number (3,211) describes neither.

## 8. Mixed-effects models: the same test at 138–721× the cost

The mixed model `expression ~ condition + (1 | donor)` on one pseudobulk profile per donor and condition is the paired t-test plus a donor variance component. On all seven cell types it gives identical effects and statistics wherever the donor variance is estimated, 1.02–1.09× the calls (80% of the extras on genes whose variance collapsed to zero), and the same false-positive count as pseudobulk in 220 of 220 balanced donor-swap draws. An earlier version of this section reported 12–272× null inflation; that was an implementation artefact and is withdrawn (commits `7f1e46d`, `96f8026`). A *per-cell* mixed model is not evaluated.

![boundary genes](figures/mixed_boundary_genes.png)

*Interior fits (donor variance estimated) track the paired t; boundary fits (variance collapsed to 0) run up to 1.41× higher and carry 80% of the mixed-only calls.*

Full section, tables, the donor-swap null and the withdrawn claims: [docs/mixed_model.md](docs/mixed_model.md).

## 9. The NB GLM on real data: more calls, and false positives on two cell types

On the same pseudobulk matrices the negative-binomial GLM calls 7–91% more genes than the paired t-test and ranks the interferon genes higher, at 5.5–23× the runtime. Under a balanced donor-swap null, where pseudobulk makes none, the GLM makes false positives on B cells (35 of 35 draws) and CD14+ monocytes (14 of 35). So the extra calls are an unvalidated gain in calls, not demonstrated power. In simulation (NB counts, the GLM's own model) it finds 97% of changed genes at 6 donors against 68% for the paired t.

![NB GLM null](figures/nbglm_null_anticonservative.png)

*On B cells every one of the 35 balanced draws gives the NB GLM at least one false positive (2.11 per draw) and on CD14+ monocytes 14 of 35 do (1.14 per draw); pseudobulk gives none on either.*

Full section (call tables, log₂FC agreement, the 17 separated fits, the ISG ranks, the coin-flip families) and the simulated donor grid: [docs/nb_glm.md](docs/nb_glm.md).

## How the simulation got fixed

The first version of this study produced the opposite of the published result: the per-cell test looked perfectly calibrated and pseudobulk looked broken. Both were artefacts. The simulator applied each donor's random effect identically to that donor's control and stimulated samples and gave every sample the same number of cells, so donor effects cancelled exactly when cells were pooled and the per-cell test had nothing to trip over. Adding an independent per-(donor, condition) effect (`sample_sd`) reproduced the known behaviour at once. The failure was quiet: the code ran, produced plausible tables, and said something false. A known answer to check against caught it.

## Layout

```
src/scpipe/
    pseudobulk.py   aggregation, median-of-ratios size factors, filtering
    de.py           per-cell Wilcoxon, pseudobulk t-test, BH correction
    simulate.py     negative-binomial counts with donor and sample effects
    preprocess.py   QC thresholds, standard scanpy pipeline
    data.py         real-dataset loading
scripts/
    simulation_study.py     the four simulation experiments
    false_positive_test.py  null experiment on real data
    real_analysis.py        ctrl vs stim per cell type, positive control
    nbglm_study.py          GLM vs t-test, simulated (power, FDP)
    nbglm_real.py           GLM vs t-test, real Kang counts (§9)
    mixed_vs_others.py      mixed-model comparison, 7 cell types (§8)
    bh_null_expectation.py  BH all-null reference simulation (§8)
    celltype_analysis.py    stratified analysis, composition shift (§7)
    megakaryocyte_check.py  why megakaryocytes get no expression test (§7)
    make_figures.py         every figure, from results/*.csv
    audit_claims.py         every README number against its results/ artefact
docs/               §8 and §9 in full, and the simulated GLM donor grid
results/ figures/   the CSVs behind every number; every plot
tests/              124 tests, including calibration tests that assert the
                    per-cell test fails and stops failing when sample-level
                    variance is removed, and driver tests that pin the
                    seven-type roster and the balanced null
```

## Reproducing

```
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps   # exact versions; or: pip install -e ".[dev]"
pytest -q                                   # 124 tests, no dataset needed
python scripts/audit_claims.py --strict     # every README number against results/
# simulation (no download)
python scripts/simulation_study.py          # §1-§4
python scripts/nbglm_study.py               # simulated GLM vs paired t donor grid
python scripts/bh_null_expectation.py       # BH all-null reference
# real data (downloads Kang 2018, 38.4 MB, into data/ on first use)
python scripts/false_positive_test.py --cell-type "CD14+ Monocytes"   # §5
python scripts/real_analysis.py             # §6
python scripts/celltype_analysis.py         # §7
python scripts/megakaryocyte_check.py       # §7
python scripts/mixed_vs_others.py           # §8; --types "CD4 T cells" for one
python scripts/nbglm_real.py                # §9
python scripts/make_figures.py              # every figure from results/*.csv
```

Expected runtime on one desktop CPU (12 threads, Python 3.12, measured for this revision in a fresh clone): `pytest` about 50 s and the auditor under 1 s; the simulation scripts about 2 min (`simulation_study.py` 71 s, `nbglm_study.py` 34 s); `false_positive_test.py`, `real_analysis.py`, `celltype_analysis.py` and `megakaryocyte_check.py` under 1 min each, including the download; `make_figures.py` 9 s. Two scripts are slow: `mixed_vs_others.py` took 14 min for CD4 T cells, 9 min for CD14+ monocytes and 9 min for B cells (plan on most of an hour for all seven types; `--types` runs one), and `nbglm_real.py` runs for over 15 min, mostly in its donor-swap nulls (neither was timed to completion for this revision).

Tested in CI on Python 3.10 (lower bounds from `pyproject.toml`) and 3.12 (pinned `requirements.txt`); other versions are untested. `make_figures.py` skips a figure whose CSV is missing. The dataset URL `figshare.com/ndownloader/files/34464122` looks like a download link and is not: it is an async dispatch endpoint that answers `202` with `Content-Length: 0`, so a naive download writes an empty file and fails later with an unhelpful h5py error. `src/scpipe/data.py` uses the scverse host and checks the HDF5 magic bytes and md5.

## Note on how this was built

Written with Claude (Anthropic) as a coding assistant. I set the question, chose the validation strategy, and interpreted the results. The three findings I would defend in detail (the per-cell failure is caused specifically by unmodelled sample-level variance, composition bias grows in relative importance as a study gets larger, and p-value ranking is close to arbitrary at these replicate counts) all came from following up results that did not match what was expected, instead of accepting output that looked reasonable.

## Part of a set

Six repos, one method: take the number a field reports, measure what a practitioner uses, publish
the gap.

- [admet-split-truth](https://github.com/Qwuantabee/admet-split-truth): does a validation split predict error on the next compounds?
- [mlip-fidelity](https://github.com/Qwuantabee/mlip-fidelity): do ML potentials rank near-degenerate conformers correctly?
- [epistasis-ceiling](https://github.com/Qwuantabee/epistasis-ceiling): do single-mutant effects add up on a double mutant?
- [variant-effect-prediction](https://github.com/Qwuantabee/variant-effect-prediction): does a cleaner assay target help zero-shot ESM-2?
- single-cell-de: are cells replicates?
- [enzyme-kinetics](https://github.com/Qwuantabee/enzyme-kinetics): how wrong is a Lineweaver-Burk Km?

Overview: https://github.com/Qwuantabee
