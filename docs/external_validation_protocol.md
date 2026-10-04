# Independent Public-Cohort Extension

Protocol recorded on 2026-10-04 before testing expression outcomes in the
additional cohorts. This is a dated analysis specification, not a registered
prospective protocol.

## Questions

1. Do the existing spatial tissue-state associations recur in independent human
   pancreatic tissue, including a chronic-pancreatitis comparator?
2. Is the Normal-to-ADM expression direction compatible with a paired human
   acinar-cell culture model of ADM?
3. Which circulating-miRNA markers are measurable in additional public blood
   studies, and what comparisons and clinical covariates do those studies permit?

## Frozen Reference

Use the existing GSE208536 all-target statistics without selecting new genes
using the external outcomes. Freeze the 506 omnibus-associated targets,
contrast significance flags, and effect directions. STAT1 is a named example
already examined, not an externally selected winner. Preserve the old diagnostic
plasma coefficients and PLCO results; do not tune to PLCO or reverse ROC direction.

## Eligibility and Analysis

- Human data with public processed expression and identifiable sample groups.
- Prioritize PDAC/CP within the same study. Do not compare cases from one study
  against controls from another and call batch correction an independent test.
- Average technical replicates within donor/tissue before inference. Use paired
  donor blocks where pairing is explicit. Do not infer pairing from column order.
- Keep adjacent normal, healthy donor tissue, pancreatitis, culture-derived ADM,
  serum, plasma, and EV measurements conceptually separate.
- Map features using the deposited platform annotation or versioned annotation
  packages; exclude ambiguous symbols and collapse duplicate probes by highest
  overall mean, without outcome-dependent probe selection.
- Analyze deposited log-scale microarrays with limma. For public RNA counts,
  use an established count filter, library normalization, and paired voom/limma
  model. Audit transformations before fitting.
- Correct discovery tests across all eligible genes within each study/contrast.
  Also report multiplicity across the frozen measurable spatial family, with its
  denominator explicit. Confidence intervals are pointwise unless stated otherwise.
- Report all eligible results, discordant directions, and negative findings.
  Do not remove a study for failing to replicate the preferred result.
- Blood studies on different analytes/platforms may support marker association,
  not direct coefficient transfer on incompatible scales. Stage, CA19-9, age,
  sex, treatment and control subtype analyses require actual sample-level metadata.

Before serum outcome testing, freeze the union of the historical six-miRNA
panel, the published three-miRNA panel, and the existing 20-miRNA EV candidate
pool. Univariate ROC direction comes from the diagnostic plasma mean difference,
not from optimizing the serum ROC. Serum contrasts use age/sex-adjusted limma;
report all-feature BH and frozen-family BH within and across the prespecified
healthy, benign, and other-cancer comparisons. Do not call a mixed pancreatic/
biliary benign group chronic pancreatitis.

GSE85589 deposits CA19-9 only as <37/>37, not continuous concentrations. As an
exploratory within-cohort analysis, compare age/sex/binary-CA19-9 ridge models
against the same baseline plus each frozen 3- or 6-miRNA panel. Use 10 outer,
5 inner stratified folds, lambda selection by inner binomial deviance, alpha=0,
and identical outer folds. No serum feature selection or PLCO refitting. Report
both all-available and CA19-9-low subgroup fixed-direction marker associations.
CV CIs are conditional on the resulting out-of-fold scores, not estimates of
all training/split uncertainty. Do not interpret case-control calibration as
screening-population risk or binary-CA19-9 analysis as superiority to an optimized
continuous CA19-9 test. Stage-II subgroup analyses in GSE59856 use deposited
pStage IIA/IIB and exclude ypStage (post-neoadjuvant) specimens; the absence of
yp coding alone does not establish a complete treatment history. They are
diagnostic, not prediagnostic validation. CA19-9-low subgroup ROC comparisons
restrict both cases and controls to deposited <37 values.

## Candidate Audit

GSE91035: bulk normal/benign/PDAC tissue; determine exact benign diagnoses from
deposited labels before assigning a CP group.

GSE15471: matched tumor/normal tissues, with technical replicates and a source
quality-control exclusion to resolve from sample metadata/publication.
Before fitting, exclude donor 51721 (GSM388111/GSM388150), identified as the
failed-QC pair in the Bioconductor KEGGdzPathwaysGEO documentation. Average
technical duplicates for donors 30162, 40728 and 41027, leaving 35 donor pairs.
This uses a published exclusion, not an exclusion selected from our gene results.

GSE179248: 14 paired human donors, acinar cells at culture day 0 and day 6.
This tests a culture response accompanying ADM, not longitudinal human cancer.

GSE59856 and GSE85589: serum miRNA studies; audit clinical labels, benign controls,
stage information and feature annotations before using their expression outcomes.

GSE101462 added during metadata-only audit: pancreatic tissue with pancreatitis
controls; restrict the primary cancer/pancreatitis comparison to FFPE samples if
confirmed by the deposited preservation labels. Retain fresh-frozen samples only
for explicitly labeled sensitivity analyses, not as preservation-matched controls.
The metadata audit subsequently found complete disease/chip confounding: all
10 pancreatitis specimens use Sentrix 100965680058; all six PDAC and four normal
specimens use 100965680053. Even within FFPE, the disease-plus-chip design is
rank deficient. Therefore do not fit an inferential PDAC/pancreatitis validation
model or claim that ComBat can resolve it. Retain the audit and source only.

RNA counts are mapped by unambiguous ENSEMBL-to-symbol mappings, summed where
multiple ENSEMBL IDs map to one symbol, filtered with edgeR filterByExpr (day
group), TMM normalized and analyzed with a donor-blocked voom/limma model.

The GSE91035 supplement's GeneSymbol field contains transcript labels and
date-converted gene names. The corrected mapping uses GPL22763 Entrez IDs first,
then exact current symbols or unambiguous aliases from org.Hs.eg.db. Unresolved
labels are excluded from gene-level inference; the expression values are
checked against the series matrix independently of the annotation. This repair
was made after the first exploratory run revealed annotation defects; that
superseded run is not the delivered analysis. No manual date-to-gene guessing.

E-MTAB-1791 was inspected through the public BioStudies API. Its processed files
are approximately 33 MB each for 457 samples, and pancreatitis also occurs in
pancreases bearing other tumors. Deferred in this extension pending a compact
matrix and a donor/clinical-history audit; no gene-association results were tested.

GSE28735/GSE62452 and GSE179248/GSE295071 require donor-overlap checks before
being counted as independent studies. GSE123375 is cultured fibroblasts, not
whole tissue. Mouse datasets are not human clinical validation.

## Sources

- [GSE91035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91035)
- [GSE15471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471)
- [GSE179248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179248)
- [GSE59856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59856)
- [GSE85589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85589)
- [GSE101462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101462)
- [E-MTAB-1791](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1791)
- [GSE15471 QC documentation](https://master.bioconductor.org/packages/devel/data/experiment/manuals/KEGGdzPathwaysGEO/man/KEGGdzPathwaysGEO.pdf)
