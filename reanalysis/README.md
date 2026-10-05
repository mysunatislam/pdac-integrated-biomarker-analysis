# Corrected Public-Data Reanalysis

This directory contains the current analysis. The historical project remains
under the original `scripts/` and `results/` paths.

Read [the methods/results report](../docs/reanalysis_methods_results.md) first.
It is generated from the final tables, so reported counts agree with the
deposited outputs. [Claims and limitations](../docs/reanalysis_claims_and_limitations.md)
explains how the manuscript must change.

The [additional-cohort report](../docs/external_validation_methods_results.md)
documents paired human tissue/culture checks and independent serum analyses.
Its frozen families and cohort exclusions are specified in
[the dated protocol](../docs/external_validation_protocol.md).

The [publication audit](../docs/publication_audit_methods_results.md) adds a
source-paper case crosswalk, immune/stromal composition sensitivities and the
completed larger inflammatory-cohort analysis (195 PDAC / 59 CP specimens). The
[focused prior-art assessment](../docs/publication_prior_art.md) and
[submission checklist](../docs/publication_readiness_checklist.md) explain the
remaining manuscript work and unavailable metadata. These are post hoc audits,
not a prospective registration or a declaration of journal readiness.

## Reproduce

Use the instructions in [scripts/reanalysis/README.md](../scripts/reanalysis/README.md).
Run from the repository root. Inputs are already in `data/source/`.
The primary scripts assert sample identity, group sizes, feature mappings, and
numerical validity. The final validation script checks the figure pairs.
The lockfile records the packages used; exact restoration has not been tested
on a second operating system.

## Evidence Status

| Layer | Current interpretation |
| --- | --- |
| Bulk PDAC versus CP | Age/sex-adjusted primary test; unadjusted results are sensitivity evidence |
| Spatial Normal/ADM/PDAC | Patient-profile-blocked tissue-state association; grouping is inferred |
| EV discovery | Exploratory ranking after explicit detection filtering; no significant hits in the stated primary analysis |
| Diagnostic plasma | Fixed historical panels plus two selection pipelines with nested CV |
| PLCO plasma | Diagnostic models applied without PLCO outcome-driven tuning |
| Secondary EV cohort | Separate descriptive analysis because assay, analyte, stage mix, and treatment differ |
| GO enrichment | Supplementary; background restricted to jointly measurable, GO-annotated genes |
| PPI and multiMiR | Historical supporting material, without a new causal or clinical claim |
| Additional tissue/culture | Paired-donor associations and directional compatibility, not proof of in-vivo progression |
| Larger inflammatory cohort | 303 of 501 measurable frozen genes pass chip-adjusted family FDR; STAT1 threshold significance is sensitive to multiplicity and CP omissions; clinical/donor limits remain |
| Additional serum | Fixed marker associations and exploratory nested-CV CA19-9 increment; no compelling panel increment |

The spatial labels are tissue states sampled at one time, not longitudinal
stages of an observed patient's cancer development. The bulk comparison uses
chronic pancreatitis as a comparator but does not demonstrate absence of
inflammation effects or malignancy-specific expression.

## Source Studies

1. Chhatriya et al. (2020). Transcriptome analysis identifies putative multi-gene signature distinguishing benign and malignant pancreatic head mass. [Journal of Translational Medicine](https://doi.org/10.1186/s12967-020-02597-1). Dataset: GSE143754.
2. Yang et al. (2022). An integrated model of acinar to ductal metaplasia-related N7-methyladenosine regulators predicts prognosis and immunotherapy in pancreatic carcinoma based on digital spatial profiling. [Frontiers in Immunology](https://doi.org/10.3389/fimmu.2022.961457). Dataset: GSE208536.
3. Xu et al. (2025). Human multiethnic radiogenomics reveals low-abundancy microRNA signature in plasma-derived extracellular vesicles for early diagnosis and molecular subtyping of pancreatic cancer. [eLife](https://doi.org/10.7554/eLife.103737). Dataset: GSE304572.
4. Treekitkarnmongkol et al. (2024). Blood-Based microRNA Biomarker Signature of Early-Stage Pancreatic Ductal Adenocarcinoma With Lead-Time Trajectory in Prediagnostic Samples. [Gastro Hep Advances](https://doi.org/10.1016/j.gastha.2024.08.002). Dataset: GSE259327.
5. Kim et al. (2025). Evaluation of Exosome-derived Small RNAs as Potential Biomarkers for Pancreatic Ductal Adenocarcinoma Using Next-generation Sequencing. [Annals of Laboratory Medicine](https://doi.org/10.3343/alm.2025.0121). Dataset: GSE268771.
6. Yao et al. (2026). Machine learning-based multimodal biomarkers enable accurate diagnosis and early detection of pancreatic ductal adenocarcinoma. [Scientific Reports](https://doi.org/10.1038/s41598-025-29548-4). Additional published analysis of GSE259327.

Public GEO records: [GSE143754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143754),
[GSE208536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208536),
[GSE304572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304572),
[GSE259327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259327),
[GSE268771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268771).
