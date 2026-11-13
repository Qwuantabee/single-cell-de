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

Correct FDR control would call about 25 or fewer.

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
    make_figures.py
tests/              19 tests, including calibration tests that assert the
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
  GLM's own model, so this is its best case; on real data the gap is unlikely
  to be larger than this and could be smaller.
- Per-cell-type analysis: aggregate within cell type, since a composition
  shift between conditions can masquerade as expression change.
- Mixed-effects models as a middle ground between per-cell and pseudobulk.

## Note on how this was built

Written with Claude (Anthropic) as a coding assistant. I set the question,
chose the validation strategy, and interpreted the results. The three
findings I would defend in detail — that the per-cell failure is
specifically caused by unmodelled sample-level variance, that composition
bias grows in relative importance as a study gets larger, and that p-value
ranking is close to arbitrary at these replicate counts — all came from
following up results that did not match what was expected, instead of
accepting output that looked reasonable.
