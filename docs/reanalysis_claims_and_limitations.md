# Revised Claims and Publication Gaps

## Claims Supported by This Reanalysis

The deposited spatial data contain tissue-state associations after reducing
multiple ROIs to patient-profile/state means and blocking on inferred profile.
Specific Normal-to-ADM and ADM-to-PDAC contrasts describe cross-sectional
tissue differences. They do not prove a temporal sequence, causality, or
cell-type-specific regulation.

The deposited diagnostic and PLCO cohorts allow a direct assessment of
transportability on the same assay. Under the documented preprocessing and
refitting strategy, diagnostic discrimination is much stronger than
prediagnostic discrimination. None of the tested panels establishes a useful
prediagnostic classifier in the full deposited PLCO sample.

Age/sex adjustment materially weakens the bulk tissue findings. The data
therefore do not establish a robust CP-specific discovery set. This may reflect
confounding, loss of precision in a small sample, model assumptions, or a
combination. Failure to pass FDR is not proof of no biological effect.

## Claims to Withdraw or Narrow

| Historical claim | Revised treatment |
| --- | --- |
| A validated malignant-switch gene set | Report the direct adjusted PDAC-versus-CP contrast; do not infer specificity from a non-significant CP-versus-normal test |
| Independent spatial n = 48 | Report 48 ROIs and eight inferred patient profiles |
| Six ROIs and all three states in every patient | Deposited ROI counts are uneven; only six profiles contain all three states |
| STAT1 is a significant transient ADM peak | Its Normal-to-ADM increase is supported; its subsequent decrease is not significant |
| 82 healthy plasma controls | The source studies report 46 healthy and 36 pancreatitis controls; chronic pancreatitis is a subset of pancreatitis |
| 15 normal controls in GSE268771 | Report three healthy and 12 benign participants |
| The six-miRNA CV AUC is fully selection-adjusted | Historical feature selection occurred before this CV; retain it as a conditional benchmark |
| The published three-miRNA signature failed replication | This is a refit using our preprocessing and all 96 deposited PLCO samples, not the original published score or lead-time subset |
| Tissue genes and plasma miRNAs form a validated regulatory mechanism | The cohorts are unpaired and database links do not establish reciprocal regulation |
| PPI hubs are the principal discovery | Keep the old network analysis as historical supporting material |

## Remaining Constraints

1. GEO does not provide explicit spatial patient identifiers. Clinical-profile
   grouping is plausible and matches eight reported patients, but must be
   confirmed before a manuscript treats it as definitive.
2. GSE143754 has six CP specimens and notable age differences. A linear age term
   and binary sex adjustment cannot eliminate all confounding. The source
   matrix was normalized by its authors; raw CEL normalization was not repeated.
   Adjacent-normal/tumor pairing is not identifiable from unique patient IDs
   in the deposited metadata. The primary adjusted and unadjusted PDAC-versus-CP
   models therefore use only those 17 specimens. All-group and normal contrasts
   are exploratory; their independence assumption cannot be verified.
3. The GeoMx assay measures a cancer-focused panel. Results cannot be generalized
   to the complete transcriptome. Patient-profile blocking absorbs profile-level
   stable differences, including age/sex; it cannot identify cell composition.
4. GSE259327 sample-level control subtypes, cancer stage, age/sex, CA19-9,
   PLCO lead time, and matched-pair identifiers are absent from the downloaded
   public files. Their absence prevents trustworthy subgroup performance,
   clinical adjustment, paired-cohort modeling, and CA19-9 benchmarking.
   The [newer source-cohort paper](https://doi.org/10.1038/s41598-025-29548-4)
   states that clinical, metabolite, and CA19-9 data require a request to the authors.
5. PLCO is a balanced case-control sample. Calibration intercepts, Brier scores,
   and observed proportions describe that sample and are not population
   screening-risk calibration. Decision curves would require a justified
   target prevalence and intended clinical setting; none is claimed here.
6. Cross-validation uses one deterministic outer split. CIs calculated from
   out-of-fold scores are conditional intervals and omit uncertainty from
   repeated pipeline training and split selection.
7. PLCO has now been examined. Further model redesign must not treat it as an
   untouched confirmatory set. A future model requires a new independent
   public validation cohort or explicit exploratory status.
8. EV filtering intentionally excludes very sparse miRNAs. It cannot refute
   source-paper claims based on low-abundance features, radiomics, or different
   filtering and normalization.
9. GSE268771 is dominated by advanced cancer and includes treatment exposure.
   Its differential results do not validate early detection.

## Manuscript Direction

The [independent-cohort extension](external_validation_methods_results.md) adds
35 matched tumor/normal donor pairs, 14 paired human culture donors, exploratory
bulk contrasts, and two serum studies. STAT1 increases recur in paired tumor
and culture data. This strengthens tissue-state association evidence without
establishing the claimed transient peak or CP specificity. GSE91035 has only
two CP samples; GSE101462 has perfect disease/chip confounding and was not fitted
as validation. Fixed serum panels do not materially improve the binary-CA19-9,
age/sex baseline in the separate GSE85589 nested-CV analysis. These new data do
not change the locked PLCO results or establish prediagnostic performance.

A defensible working title is **Patient-aware tissue-state associations and
prediagnostic transportability of plasma miRNA models in public PDAC cohorts**.
The revised paper should be presented as an exploratory reproducibility study.
The methods and results in the companion report can replace the corresponding
sections, but a full journal manuscript still needs an author-reviewed
Introduction, Discussion, declarations, and a focused novelty assessment.

Before submission, confirm spatial patient identities, verify source-paper
preprocessing and score reconstruction, and seek independent public tissue
evidence with an inflammatory comparator. These are substantive evidence gaps,
not formatting tasks. Neither a journal quartile nor acceptance can be inferred
from the current reanalysis.
