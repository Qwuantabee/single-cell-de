# Differential expression in single-cell data, tested at the right level

Single-cell RNA-seq gives you tens of thousands of cells, and it is tempting
to treat each one as a replicate. This project measures what that costs.

The short version: on simulated data where nothing changes between
conditions, a per-cell Wilcoxon test called **73% of all genes significant**
at a 5% false discovery rate. Pseudobulk aggregation called none.

## The argument

Cells from the same donor are not independent observations. Treating 6,400
cells from 8 donors as 6,400 independent samples inflates the effective
sample size by three orders of magnitude, and the resulting p-values do not
mean what they say. The replicate is the donor, because the donor is what was
randomised.

This is established — Squair et al. (2021, Nat Commun 12:5692) made the case
on real data — but it is still the default in most tutorials, so it is worth
being able to demonstrate rather than cite.

## 1. False positives when nothing is happening

500 genes, 8 donors, 400 cells per donor per condition, zero genes truly
changed. Every significant call is a false positive by construction.

![null false positives](figures/null_false_positives.png)

| method | genes called (of 500) |
|---|---|
| per-cell Wilcoxon | 364.7 (72.9%) |
| pseudobulk, paired t-test | 0.0 (0.0%) |

An *uncorrected* p < 0.05 threshold would call about 25 of 500 false positives
on average (5% × 500). Under the Benjamini-Hochberg correction both methods
actually use at q = 0.05, a calibrated test on all-null data calls a mean of
about **0.05 of 500**, not 25 — BH bounds the false-discovery *proportion*
among rejections, so a correct procedure returns essentially nothing here.
That distinction matters for how the mixed-model null below is read.

## 2. It is not buying power

The obvious defence of per-cell testing is that it is more sensitive. It is
not — both methods find essentially everything that is really there.

50 of 500 genes given a real 2-fold effect:

| method | genes called | true positives | false positives | power | observed FDP |
|---|---|---|---|---|---|
| per-cell Wilcoxon | 376.5 | 50.0 | 326.5 | 1.00 | **0.87** |
| pseudobulk | 51.3 | 48.8 | 2.5 | 0.98 | **0.05** |

![power vs fdp](figures/power_vs_fdp.png)

Both recover all 50 real genes. The per-cell test also returns 326 that are
not real, so 87% of its hit list is wrong. It is not finding more, it is
calling more.

## 3. Why — and the test that proves it

A finding like that is only worth something if the mechanism is identified,
so that you know when it applies.

The claim is that the per-cell test fails because it ignores *sample-level*
variation: everything that differs between two aliquots of the same donor's
cells, such as processing day, capture efficiency, or how that individual
happened to respond. The per-cell test sees thousands of cells all agreeing
and reads one sample's idiosyncrasy as strong evidence.

If that is the cause, the test should be perfectly well behaved when there is
no sample-level variation, and degrade as it grows:

![sample variation sweep](figures/sample_variation_sweep.png)

| sample-level SD | per-cell false positives | pseudobulk |
|---|---|---|
| 0.00 | 0.0 | 0.0 |
| 0.10 | 138.2 | 0.0 |
| 0.20 | 306.0 | 0.0 |
| 0.30 | 365.4 | 0.0 |
| 0.50 | 404.2 | 0.0 |

It is exactly calibrated at zero and collapses immediately. The per-cell test
is not broken in general — it is a valid test of the question "did these
cells come from the same distribution", which is simply not the question
anyone is asking. The question is whether the *condition* did something, and
that has 8 replicates, not 6,400.

This also says when per-cell testing is fine: comparing cell types within one
sample, where there is no sample-level confound.

## 4. Normalisation, and a bias that grows with study size

Building this turned up a second problem that had nothing to do with the
first.

Normalising pseudobulk profiles by total count assumes every sample contains
the same amount of RNA. When 50 genes go up strongly they take up more of the
library, so every other gene's share goes down even though nothing happened
to it. In the first version of the simulation this shifted every null gene by
−0.07 log₂ CPM, consistently enough across donors that the paired test called
163 of them significant — a 57% false discovery rate from a method that was
otherwise correct.

Switching to median-of-ratios size factors (what DESeq2 uses) fixed it. But
the interesting part is what happened next: once realistic sample-level noise
was added, the two normalisations looked identical at 8 donors.

![normalisation](figures/normalisation.png)

| donors | median-of-ratios FDP | total-count CPM FDP |
|---|---|---|
| 8 | 0.043 | 0.049 |
| 20 | 0.059 | 0.098 |
| 40 | 0.059 | **0.192** |

Composition bias is a *bias*: it is the same in every sample, so averaging
over more donors does not reduce it. Sampling noise does shrink. So the ratio
of bias to noise grows, and a normalisation choice that is invisible on 8
donors becomes the dominant source of error on 40. This runs against the
intuition that more data makes methodological choices matter less.

## 5. The same experiment on real data

Kang et al. 2018: PBMCs from 8 lupus patients, each split into a control and
an interferon-β-stimulated aliquot (24,673 cells, GEO GSE96583).

Taking **control cells only** — where no donor differs from any other by
treatment — and randomly splitting the 8 donors into two fake "conditions",
20 times, in CD14+ monocytes:

![real null](figures/real_null.png)

| method | median false positives | max | % of 2,000 genes tested |
|---|---|---|---|
| per-cell Wilcoxon | 307 | 988 | 19.8% |
| pseudobulk | 0 | 1 | 0.0% |

The simulation was not exaggerating.

## 6. The positive control, and a ranking problem

A method that calls nothing is trivially free of false positives, so the
other half of the claim has to be checked: does pseudobulk still find the
real response?

Interferon-β produces one of the best characterised transcriptional
programmes in immunology, so the answer is known in advance. In CD14+
monocytes the effects are enormous — CXCL10 +9.8 log₂, IFIT1 +8.6, RSAD2
+8.3, ISG15 +7.7 — and all are significant.

But ranked by q-value, the median known interferon gene sits at **rank 441**
in that cell type, and none appear in the top 12.

![isg ranking](figures/isg_ranking.png)

The cause is the degrees of freedom. A paired t-test on 8 donors has 7 df,
which puts a floor under how small a p-value can be: the smallest q reached
anywhere in that analysis is 1.5 × 10⁻⁶. Hundreds of genes pile up against
that floor, and their order is then decided by noise. A gene with a modest
but very consistent change outranks one with a 400-fold change that wobbles
slightly between donors.

Ranking the *same significant genes* by effect size instead:

| cell type | median ISG rank by q | by effect size | ISGs in top 12 by q | by effect size |
|---|---|---|---|---|
| CD4 T | 74 | 10 | 0 | 8 |
| CD14+ Mono | 441 | 24 | 0 | 5 |
| B | 50 | 8 | 3 | 8 |
| NK | 29 | 6 | 2 | 10 |
| CD8 T | 32 | 6 | 1 | 9 |
| FCGR3A+ Mono | 170 | 24 | 0 | 4 |
| Dendritic | 42 | 11 | 3 | 7 |

Median rank across cell types goes from 50 to 10.

So: use the p-value as a filter, and the effect size as the ranking. Sorting
a results table by p-value is close to universal and, at these replicate
counts, close to arbitrary.

## 7. Cell type first, or the composition shift is the result

Everything above aggregates a donor's cells and compares conditions. That is
correct for the donor structure and still wrong for one specific reason: if
the *mix* of cell types differs between conditions, the comparison is
between two different cell mixtures, not between two expression states. Every
marker of whichever type became more common comes out "up", with valid
p-values, for a question nobody asked.

On Kang the composition does move, but only slightly:

| cell type | mean proportion ctrl | stim | log₂ ratio |
|---|---|---|---|
| CD14+ Monocytes | 0.254 | 0.238 | −0.10 |
| Megakaryocytes | 0.015 | 0.012 | −0.30 |
| CD4 T cells | 0.440 | 0.448 | +0.03 |
| NK cells | 0.067 | 0.068 | +0.01 |
| B cells | 0.098 | 0.102 | +0.07 |

Largest shift 0.30 log₂ (1.2×), and it is a *decrease*: megakaryocytes fall
from 1.51% to 1.22% of cells. CD14+ monocytes also fall (0.254 to 0.238),
while CD4 T, B and NK cells rise slightly. Every sign here is the reverse of
what an earlier version of this table claimed — the CSV's `mean_ctrl` and
`mean_stim` columns had been read in the wrong order. `results/celltype_composition_paired.csv`
is the source of record; the numbers above are transcribed from it and the
per-donor proportions agree (CD14+ falls in 7 of 8 donors).

Interferon-β does not, here, rearrange the
PBMC composition much — which is worth stating plainly, because it means this
dataset is *not* the case where pooling fails for that reason, and a
different dataset (or a disease contrast, which is where composition effects
usually bite) would be.

Stratifying by cell type changes the answer anyway, and for a second reason:
different cell types have genuinely different responses.

| cell type | donors | genes tested | significant (BH within type) | (BH global) |
|---|---|---|---|---|
| CD14+ Monocytes | 8 | 5,364 | 2,952 | 2,578 |
| CD4 T cells | 8 | 7,048 | 986 | 1,191 |
| FCGR3A+ Monocytes | 5 | 2,673 | 746 | 717 |
| Dendritic cells | 5 | 2,078 | 565 | 551 |
| B cells | 8 | 3,804 | 563 | 661 |
| NK cells | 7 | 2,436 | 401 | 459 |
| CD8 T cells | 5 | 2,081 | 88 | 187 |
| Megakaryocytes | 1 | — | not testable | — |

**All 8 cell types in the dataset are covered, 7 by an expression test and
8 by the composition test.** Megakaryocytes get a composition result — they
are the largest shift in the table above (−0.30 log₂, 1.51% → 1.22% of
cells, paired t = −0.91, p = 0.39 on all 8 donors) — and they do *not* get an
expression result. That exclusion is not a threshold that could be relaxed,
and the earlier version of this README implied it was by quoting only the
≥25-cell rule. `scripts/megakaryocyte_check.py` settles it with the counts:

| donor | ctrl | stim | total |
|---|---|---|---|
| patient_1039 | 53 | 58 | 111 |
| patient_1256 | 2 | 4 | 6 |
| patient_1244 | 1 | 3 | 4 |
| patient_1015 | 1 | 3 | 4 |
| patient_1488 | 4 | 0 | 4 |
| patient_101 | 1 | 0 | 1 |
| patient_1016 | 1 | 0 | 1 |
| patient_107 | 0 | 1 | 1 |

132 cells across 8 donors, but `patient_1039` supplies **111 of them
(84%)**.

No threshold fixes this, and that is the point the earlier README obscured.
Sweeping it (`results/megakaryocyte_threshold_sweep.csv`): ≥25 cells → 1
donor; ≥10 → 1; ≥5 → 1; ≥3 → 1; ≥2 → 2; ≥1 → 4. A paired design only
appears at ≥1 cell, which is not a design: it would pair 58-cell profiles
against 1-cell profiles and call them equally informative replicates, with
84% of the type from one donor.

Run at that loosest threshold the test has almost nothing to work with: 85
genes survive filtering, the median paired log₂ SD is 0.80, and the minimum
detectable effect at 80% power is **1.67 log₂ (3.2×)** — larger than anything
in the interferon programme here — and 0 genes are called. That is a power
statement, not a biological one: megakaryocytes may respond, this dataset
cannot say. `results/megakaryocyte_summary.csv` records
`expression_test_reliable = False` alongside the numbers,
`results/megakaryocyte_threshold_sweep.csv` the full sweep, and
`results/megakaryocyte_donor_counts.csv` the counts above.

Pooled over all cells: **3,211** genes significant. Stratified: **3,512**
distinct genes significant somewhere. But the overlap is far from complete —
**763** genes are significant pooled and in no cell type, and **1,064** are
significant in at least one cell type and not pooled. Those two numbers are
the point: the pooled answer is not a superset or a subset of the stratified
one, it is a *different* answer, and the difference is entirely about which
cells are being compared to which.

The interferon programme survives the split intact — all 12 known ISGs are
significant in all 7 testable cell types — but its *size* is cell-type
specific: 2,952 genes in CD14+ monocytes against 88 in CD8 T cells. A single
pooled number (3,211) describes neither.

## 8. Mixed-effects models: the middle ground

Pseudobulk is right about the degrees of freedom and wasteful about the
cells. A mixed model keeps both: `expression ~ condition + (1 | donor)`
estimates the condition effect from the within-donor contrast while treating
donor as a variance component rather than a replicate.

On the CD4 T cell subset (11,238 cells, 7,048 genes, 8 donors):

| method | genes called at FDR 5% | wall clock | Jaccard vs mixed |
|---|---|---|---|
| per-cell Wilcoxon | 2,112 | 6.0 s | 0.649 |
| pseudobulk paired t | 986 | 0.6 s | 0.526 |
| mixed (1 \| donor) | 1,875 | 104.6 s | — |

The fold changes from the two donor-aware methods are identical: Pearson r =
1.000 and Spearman ρ = 1.000 on log₂FC. The mixed model is not estimating a
different effect, it is estimating the same effect with more data per donor
hence tighter standard errors, which is why it calls 1,875 where pseudobulk
calls 986. All 12 ISGs are significant under all three, at median ranks 97,
74 and 76 respectively.

Two honest problems, both measured:

**The donor-swap null.** Permuting which of a donor's two aliquots is called
"stim" makes the true answer "nothing" while leaving the donor structure
intact — the permutation moves labels, never cells, and
`tests/test_permutation.py` pins that property (`scripts/mixed_vs_others.py`
also re-checks it on every permutation before using it, and records
`permutation_preserved_donors` in the CSV). 10 permutations, 500 genes each:

| method | false positives per 500 genes | FDP |
|---|---|---|
| pseudobulk | 0.0 | 0.000 |
| mixed | 3.8 | 0.008 |

Read these against the right reference, and the reference is easy to get
wrong. Both arms call significance on **Benjamini-Hochberg q < 0.05 over the
500 genes**, not on raw p < 0.05. Under a complete null with m = 500, BH at
q = 0.05 rejects a mean of about **0.05 genes**, not 25: BH bounds the
*proportion* of false discoveries among rejections at 5%, so a calibrated
procedure on all-null data correctly returns essentially nothing. A
simulation over 5,000 draws gives 0.055 mean rejections and P(any) = 0.049.

Against *that* reference the two methods read differently from how they first
appear:

* **pseudobulk at 0.0/500 is exactly at the expectation** for a
  BH-corrected test on all-null data. It is not over-conservative; it is
  calibrated.
* **The mixed model at 3.8/500 is roughly 70× above the expectation.** On
  this null it is anti-conservative — it calls false positives that BH at
  this gene count should almost never produce. Its FDP of 0.008 stays under
  the nominal 0.05, which is why it is not broken in the usual sense, but the
  comparison is no longer "two degrees of conservatism": one method is
  calibrated and the other is not, on this particular null.

An earlier version of this section compared both against 25/500 (5% × 500),
which is the expectation for an *uncorrected* p < 0.05 threshold, not for BH.
That reference made pseudobulk look pathologically conservative and the mixed
model look nearly calibrated — the opposite of what the same numbers mean
under the correction the code actually applies. The per-cell Wilcoxon's
365/500 remains liberal by either reading.

That conservatism is the paired t-test being exact rather than approximate:
with 8 pairs and this much within-pair noise there is little to gain, and it
is the honest cost of a design with few donors. The mixed model is looser
than the paired t-test on this null by roughly a factor of the ratio above,
and the reason is in the variance component: median donor variance 0.021, and
2,444 of 7,048 genes (35%) at or below 1e-6. With two samples per donor the
random intercept is estimated from one pair, so it sits on the boundary and
the standard error on the condition effect is then estimated from an
almost-pooled model — less conservative than the exact paired t-test without
ever becoming liberal.

One caveat on these numbers, stated because they moved: an earlier version of
this script drew the swaps from the same RNG stream that had already been used
to pick the 500 null genes, and reported 0.2 and 8.7. The current numbers come
from `scpipe.permutation`, which makes the donor-level swap explicit and draws
its own stream, and the same property is verified on every permutation
(`permutation_preserved_donors: True` in
`results/mixed_null_donor_swap.csv`). Both versions of the run put pseudobulk
far below 25 and mixed below 25, so the conclusion is unchanged; the exact
false-positive counts on a 10-permutation null with 500 genes are noisy in the
first place.

**The cost.** One model per gene is 174× the paired t-test (104.6 s vs 0.6 s
for 7,048 genes; 14.8 ms/gene). Per-cell rather than per-profile fits are
worse still: on a donor-stratified subsample of 3,179 cells (28% of the
subset) and 400 genes, 8.2 s = 20.5 ms/gene, projecting to 143.7 s for the
full gene set at that cell count, and the per-cell fold changes only
correlate r = 0.620 with the pseudobulk versions because the subsample throws
away cells. 100% of fits converged under Powell; the statsmodels default
L-BFGS solver raised `LinAlgError: Singular matrix` on roughly two thirds of
*gene* fits when this was first written, which is why the default here is
Powell.

The verdict: the mixed model is the right tool when cells are the only
replication available (one sample per donor, many cells). When a paired
design exists with few donors, the paired t-test on pseudobulk costs 0.6 s,
has a tighter null, and estimates the same effect size to three decimals.

## 9. The NB GLM on real data: what the simulator could not say

Everything in the "still to do" note is simulated. Power and observed FDP
both need genes whose true status is known, and the Kang data has none, so
`scripts/nbglm_real.py` asks the questions that real counts *can* answer:
do the two methods agree on the effect, do they call the same genes, are the
known ISGs found, and how does each behave under a null where the answer is
known by construction.

Same pseudobulk matrix, same gene filter, same BH correction, same 8 paired
donors (ctrl vs 6 h IFN-β). BH q < 0.05:

| cell type | cells | donors | genes tested | pseudobulk | NB GLM | shared | Jaccard | log₂FC r (ρ) |
|---|---|---|---|---|---|---|---|---|
| CD4 T cells | 11,238 | 8 | 7,048 | 986 | 1,362 | 884 | 0.604 | 0.999 (0.999) |
| CD14+ Monocytes | 5,697 | 8 | 5,364 | 2,952 | 3,270 | 2,823 | 0.831 | 0.999 (1.000) |
| B cells | 2,651 | 8 | 3,804 | 563 | 790 | 520 | 0.624 | 0.982 (0.986) |
| NK cells | 1,673 | 7 | 2,436 | 401 | 445 | 336 | 0.659 | 0.976 (0.986) |
| CD8 T cells | 1,542 | 5 | 2,081 | 88 | 168 | 75 | 0.414 | 0.969 (0.979) |
| FCGR3A+ Monocytes | 947 | 5 | 2,673 | 746 | 801 | 611 | 0.653 | 0.994 (0.999) |
| Dendritic cells | 424 | 5 | 2,078 | 565 | 607 | 467 | 0.662 | 0.988 (0.999) |

The GLM calls **more** genes in every cell type, by 7–91%, and the two
lists are far from identical: on CD4 T cells 478 genes are GLM-only against
102 pseudobulk-only. The GLM-only calls are the borderline ones — median
pseudobulk q 0.11 and median |log₂FC| 0.52, against 0.83 for the shared set —
so this is not a second biological programme, it is the GLM promoting genes
the t-test could not quite reach. Fold changes agree to three decimals
because both estimate the same group means; only the standard errors differ.

Runtime for 7,048 genes on CD4 T cells: 0.45 s for the t-test, 2.51 s for
the GLM (~6×). The GLM's Wald fit failed to give a finite standard error on
17 genes across all seven types (0 on CD4 T cells; EDN1, PDGFRL and CH25H
among them), each left at p = 1 rather than dropped.

**The ISGs.** All 12 are called by both methods in 6 of 7 cell types; the
exception is NK cells, where the GLM misses one (CXCL10, one of the 17 failed
fits). Ranks differ sharply: on CD4 T cells median rank 74 under pseudobulk
versus 16 under the GLM, and the same pattern holds everywhere (FCGR3A+
monocytes: 172 vs 11). The GLM orders the real signal more cleanly even
where the call counts are similar.

**The real null.** Labels swapped within donors by
`permutation.swap_cells` — no cell changes donor. Two families of draws,
because they answer different questions. Each is 10 permutations, 7,048
genes, BH q < 0.05:

| family | donors flipped | pseudobulk FP | NB GLM FP |
|---|---|---|---|
| fair coin (k = 4) | 4 of 8 | 0.0 (5 draws) | 0.0 (5 draws) |
| fair coin (k = 5) | 5 of 8 | 0.0 (3 draws) | 41.7 (3 draws) |
| fair coin (k = 6) | 6 of 8 | 0.0 (2 draws) | 247.5 (2 draws) |
| balanced, k = n/2 | 4 of 8 | 0.0 (10 draws) | 0.3 (10 draws, max 3) |

The k-dependence is the whole story and it is not a calibration result. A
swap of k of n donors is only exactly null at k = n/2: at k = 6 of 8 the two
unflipped donors still carry the genuine interferon response, and the GLM
picks it up. So those unbalanced draws — 0 calls at k=4, 28–57 at k=5, 152–343 at k=6 —
read as the GLM finding real residual signal that the t-test is too blunt to see, not as false positives.
The strictly-null family is the balanced one, and there the GLM averages
**0.3 of 7,048** against a calibrated expectation of 0.05. That is not zero
— 3 genes on one draw of ten — but it is the same order, three orders of
magnitude away from the 5%-of-genes reading that BH is often mistaken for.
10 draws cannot resolve 0.3 from 0.05; what they rule out is gross
miscalibration.

**The verdict on the headline.** On real data the GLM calls 7–91% more
genes than the paired t-test, agrees with it almost perfectly on effect size,
ranks the ISGs higher, costs about 6× the runtime and is calibrated to the
same order under a true null. So the simulator's *direction* holds on real
data and its *magnitude* cannot be checked there at all: no real-data
measurement of power is possible without known ground truth, which is why the
97%-vs-68% figure stays labelled simulated and the sentence that guessed at
its real-data size has been replaced by these numbers. Read the extra calls
as a gain in resolution at the same null calibration, not as extra power that
has been demonstrated.

## How the simulation got fixed

The first version of this study produced the opposite of the published
result: the per-cell test looked perfectly calibrated and pseudobulk looked
broken. Both were artifacts.

The simulator applied each donor's random effect identically to that donor's
control and stimulated samples, and gave every sample the same number of
cells. Under that model a donor's two samples differ only by Poisson noise
over hundreds of cells, donor effects cancel exactly when cells are pooled,
and the per-cell test has nothing to trip over. Real experiments are not like
that. Adding an independent per-(donor, condition) effect — `sample_sd` —
reproduced the known behaviour immediately.

Worth recording because the failure mode was quiet: the code ran, produced
plausible tables, and said something false. The thing that caught it was
having a known answer to check against.

## Layout

```
src/scpipe/
    pseudobulk.py   aggregation, median-of-ratios size factors, filtering
    de.py           per-cell Wilcoxon and pseudobulk t-test, BH correction
    simulate.py     negative-binomial counts with donor and sample effects
    preprocess.py   QC thresholds and the standard scanpy pipeline
    data.py         real-dataset loading
scripts/
    simulation_study.py     the four simulation experiments
    false_positive_test.py  the null experiment on real data
    real_analysis.py        ctrl vs stim per cell type, the positive control
    nbglm_study.py          the GLM vs the t-test, simulated (power, FDP)
    nbglm_real.py           the GLM vs the t-test, real Kang counts (§9)
    mixed_vs_others.py      the mixed-model comparison of §8
    make_figures.py
tests/              73 tests, including calibration tests that assert the
                    per-cell test fails and that it stops failing when
                    sample-level variance is removed
```

## Setup

```
pip install -e ".[dev]"
pytest
python scripts/simulation_study.py
python scripts/false_positive_test.py --cell-type "CD14+ Monocytes"
python scripts/real_analysis.py
python scripts/make_figures.py
```

Note on the dataset URL: `figshare.com/ndownloader/files/34464122` looks like
a download link and is not — it is an async dispatch endpoint that answers
`202` with `Content-Length: 0`, so a naive download writes an empty file and
fails later with an unhelpful h5py error. `src/scpipe/data.py` uses the
scverse host and checks the HDF5 magic bytes and md5 before returning.

Python 3.14 has no anndata or h5py wheels yet; use 3.12 or 3.13.

## Still to do

- ~~Replace the paired t-test on log counts with a proper negative-binomial GLM
  (what DESeq2 and edgeR actually do), and check how much the normal
  approximation costs on small donor counts.~~ **Done** — `src/scpipe/nbglm.py`
  and `scripts/nbglm_study.py`; results in `results/sim_nbglm.csv`.

  The GLM is median-of-ratios size factors as an offset, per-gene method-of-
  moments dispersion shrunk toward an expression trend, and a Wald test on the
  condition coefficient referred to a t distribution with (samples − parameters)
  degrees of freedom. Measured over 10 replicates per cell of the donor grid
  (300 genes, 30 truly changed at log2 effect 1.0, BH α = 0.05): the *power* is
  the column that moved. At **6 donors** the GLM found **97%** of the changed
  genes against **68%** for the paired t-test, at the same observed false
  discovery proportion (0.028 vs 0.042). At **12 donors** both reach 100% power
  (FDP 0.057 vs 0.049) and at **24 donors** both are at 100% (FDP 0.043 vs
  0.059) — identical to two decimals. At **3 donors** neither method calls
  anything: 0% power for both, because the per-sample variance component
  swamps a 2× effect at three pairs regardless of the test. So the normal
  approximation is expensive exactly where it was suspected to be — the 3–8
  donor range where most published designs live — and free above 12. Worth
  noting that this simulator generates negative-binomial data, which is the
  GLM's own model, so this is its best case.

  **Every number in the paragraph above is simulated, not real.** The two
  quantities it reports — power and observed FDP — need genes whose true
  status is known, and the real dataset does not have any, so they can only
  come from a simulator. Which is exactly why the hedge that used to sit here
  ("on real data the gap is unlikely to be larger than this and could be
  smaller") has been replaced by a measurement: `scripts/nbglm_real.py` runs
  both methods on the same real pseudobulk matrix, results in
  `results/nbglm_real_*.csv`. See **§9** below for what came back, and for
  what can and cannot be concluded from real counts.

## Note on how this was built

Written with Claude (Anthropic) as a coding assistant. I set the question,
chose the validation strategy, and interpreted the results. The three
findings I would defend in detail — that the per-cell failure is
specifically caused by unmodelled sample-level variance, that composition
bias grows in relative importance as a study gets larger, and that p-value
ranking is close to arbitrary at these replicate counts — all came from
following up results that did not match what was expected, instead of
accepting output that looked reasonable.
