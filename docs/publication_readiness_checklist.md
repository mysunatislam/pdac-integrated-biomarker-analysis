# Submission Readiness Checklist

Updated 2026-10-05. This is an author-facing work checklist, not a completed
reporting-guideline submission form or a prediction of journal acceptance.
The completed E-MTAB-1791 analysis and its limitations are recorded in the generated
[publication audit](publication_audit_methods_results.md).

Numerical/output checks and visual review of both new PDF figures are complete.
The full larger-cohort script was rerun from the committed matrices without
changing its result tables; see the
[same-machine reproduction record](../reanalysis/provenance/publication_reproduction_check.txt).
This does not replace a clean restoration on another machine or OS.

## Scientific Claims

| Item | Evidence now available | Required treatment in the manuscript |
| --- | --- | --- |
| Spatial repeated observations | 48 ROIs aggregated to 22 profile/state means with blocking and sensitivity checks | Report eight inferred clinical profiles; do not use 48 independent patients or assume six ROIs in every case |
| Spatial identity crosswalk | Eight unique source-paper clinical profiles match the deposited metadata | Describe corroborated inference, not directly deposited ROI/patient identifiers |
| STAT1 pattern | ADM/Normal association and uncertain PDAC/ADM contrast | No confirmed transient peak, time trajectory or causal progression claim |
| Cellular composition | Named T-cell/fibroblast proxy models with marker-coverage checks | Report conditional sensitivity, failed spatial fibroblast adjustment and inability to resolve cellular origin |
| CP-specific tissue expression | Small age/sex-adjusted cohort, exploratory two-CP cohort, and a larger cohort with clinical/donor limits | Keep cohort limitations distinct; significant PDAC/CP contrasts do not by themselves prove malignancy specificity |
| Historical malignant-switch list | Original significance-list subtraction is not current inference | Remove as the discovery basis; do not restore old claims merely because a gene is significant in another cohort |
| miRNA feature selection | Selection nested within training folds; fixed panels benchmarked separately | Distinguish coefficient CV for historical panels from full nested selection performance |
| PLCO evaluation | Diagnostic models evaluated on all deposited prediagnostic samples without retuning | Report negative/weak transfer; do not match the published lead-time-specific score by tuning this already examined cohort |
| Clinical comparator | Fixed serum panels assessed against binary CA19-9 plus age/sex in shared nested CV | No compelling increment; conditional CV comparisons are exploratory, not a prospective test of clinical benefit |
| Cross-analyte evidence | Serum and EV analyses are separate from diagnostic plasma coefficients | Do not call re-estimated serum associations validation of the locked plasma model |
| Tissue-to-blood regulation | No paired tissue/blood measurements or tested mechanistic bridge | Keep evidence streams separate; database target records do not establish regulation in patients |
| Novelty | Focused primary-literature overlap assessment | Center the reproducibility/transportability question and cite source studies; do not promise first discovery |

## Reproducibility And Reporting

Use [TRIPOD+AI's official scope](https://www.tripod-statement.org/scope/) for the
prediction-model components, including regression models. Consider applicable
[STARD items](https://www.equator-network.org/reporting-guidelines/stard/) for
diagnostic-accuracy reporting. Applicability and exact manuscript page/line
locations must be completed after the manuscript is finalized. Neither
checklist is a substitute for an unbiased study design.

- Link every headline number to a saved table and the script that generated it.
- Define measured analyte, platform, normalization, eligibility, excluded records,
  repeated sampling, sample size and missing clinical fields separately by cohort.
- Report pointwise confidence intervals and clearly named multiplicity families.
- Describe preprocessing learned outside CV and the limits of conditional CV.
- Supply source URLs/hashes, selected probe mappings, model specifications,
  normalization records, package lockfile and both PNG/PDF figure exports.
- Preserve all failed/negative tests that bear on the hypothesis; do not retain
  only statistically significant results under the label "effective".
- Run the numerical/output checks and inspect final PDF renders.
- Perform a clean restoration/reproduction on a second machine or OS before
  submission. That independent environment check has not yet been performed.

## Manuscript And Author Actions

- Replace the old title, abstract, methods, results and discussion with the
  corrected study question and actual results. The original Word manuscript
  should not be submitted unchanged.
- Reconcile main/supplementary tables and captions with the current outputs.
- Credit every original dataset publication and distinguish this reanalysis
  from published scores, preprocessing, cohorts and objectives.
- Consolidate limitations into a specific section without deleting material
  uncertainty from the methods or conclusions.
- Obtain author/statistical review of the final computational results.
- Confirm authorship order, affiliations, correspondence, contributions,
  funding, conflicts and data-reuse/ethics statements. Do not invent approval
  numbers or assume that all public-data studies have identical requirements.
- Disclose software/AI assistance as required by the selected journal.
- Choose a journal with an appropriate scope for a retrospective public-data
  reproducibility study, then complete its formatting, data/code and reporting
  requirements. No journal has been selected and no submission has been made.

## Decision

The repository can support a substantially revised computational manuscript.
It does not establish a clinically validated biomarker, cellular mechanism or
guaranteed publication. Submission readiness requires the manuscript and
author actions above as well as the final analysis checks. Unavailable patient
identifiers or clinical fields must be disclosed and claims narrowed rather
than silently treated as resolved.
