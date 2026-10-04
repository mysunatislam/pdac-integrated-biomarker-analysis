# Additional Public Cohorts: Methods and Results

Generated from the saved tables on 2026-10-04. This extends the corrected analysis; it is not a submission-ready manuscript or evidence of clinical validation.

## Main Findings

Independent paired human tissue and culture data support parts of the frozen spatial tissue-state pattern, including a STAT1 increase. These are association and culture-response results. They do not demonstrate a cancer-specific malignant switch, a transient ADM peak, tissue-to-miRNA regulation, or prediagnostic biomarker performance.

In a separate serum study, the fixed miRNA panels did not materially improve nested-CV discrimination over binary CA19-9, age and sex. The prior diagnostic-to-PLCO results and coefficients were preserved without further PLCO tuning.

## Cohort Eligibility

| Dataset | Deposited samples | Analysis unit and role |
| --- | --- | --- |
| GSE15471 | 78 arrays; 36 tumor/normal donor pairs and technical replicates | Average technical replicates; exclude published failed-QC donor 51721; analyze 35 donor pairs |
| GSE179248 | 28 RNA-seq samples | 14 paired donors, acinar cells at culture day 0/day 6; culture response accompanying ADM |
| GSE91035 | 50 arrays: 25 PDAC, 2 CP, 15 adjacent benign, 8 commercial normal RNA | Separate exploratory two-group contrasts; source overview's 27 PDAC count does not match deposited disease labels |
| GSE59856 | 571 serum samples: 100 PC, 150 healthy, 21 mixed benign pancreatic/biliary, 300 other cancers | Marker associations; 17 pStage IIA/IIB cases for a diagnostic subgroup |
| GSE85589 | 232 serum samples: 88 PC, 19 healthy, 10 cholelithiasis, 115 other cancers | Marker associations and exploratory binary-CA19-9 increment analysis |
| GSE101462 | 20 tissue arrays: 6 PDAC, 4 normal, 10 pancreatitis | Eligibility audit only: disease/chip perfectly confounded, including the 3-PDAC/10-pancreatitis FFPE subset |

GSE59856's benign group is not a CP-only cohort. GSE85589's 29 non-cancer controls are 19 healthy plus 10 cholelithiasis, not 29 healthy. PC is the deposited pancreatic-cancer label; sample-level histologic subtype is not used to relabel those serum participants as proven PDAC. Distinct datasets are treated as separate studies; no documented participant overlap was identified among the analyzed cohorts, but explicit cross-study donor identifiers are unavailable.

E-MTAB-1791 was inspected at metadata/file-list level and deferred pending a compact matrix and donor/clinical-history audit. GSE179248 and GSE295071 share donors and must not be counted as two independent replications. Known GSE28735/GSE62452 overlap also requires care in later extensions.

## Methods

### Frozen Families and Annotation

The dated specification is [external_validation_protocol.md](external_validation_protocol.md). We froze the existing 506 GSE208536 omnibus-associated genes and their spatial contrast directions before testing external outcomes. The current protocol is a transparent dated specification, not a registered prospective protocol. The named STAT1 example was already examined in the original reanalysis.

Affymetrix probes use hgu133plus2.db. GSE91035's annotated supplement contains transcript IDs, zero labels, and spreadsheet-style date conversions; treating all of these as gene symbols would inflate the gene universe. The delivered analysis instead maps GPL22763 Entrez IDs to current symbols, then resolves remaining exact symbols/unambiguous aliases through org.Hs.eg.db. No date-to-gene guesses are made. Probe selection uses the highest overall mean across the study, not the disease effect. GSE91035 STAT1 is unresolved/unavailable after this annotation procedure, so no STAT1 effect is reported for that cohort. This annotation repair followed a superseded exploratory run and is disclosed in the protocol.

Serum mature MIMAT accessions are mapped using the compact GPL18941 miRBase-v20 annotation. GSE85589's _st probe suffix is removed before the accession match. Ambiguous/multiple accessions and features with any missing values are excluded; duplicate names use the highest overall mean. Thus all-feature BH is across the eligible annotated/complete miRNAs, not every array feature. Coverage denominators and probe-level exclusions are deposited. Complete mapping does not establish signal above the assay detection limit; constant/background-level markers can remain uninformative with AUC 0.5.

### Tissue Models

Deposited log2 microarrays use robust, trend-moderated limma. GSE15471 blocks explicitly on donor after technical-replicate averaging. GSE91035 uses separate PDAC-versus-CP/normal/benign two-group models. Age/sex and benign/tumor pairing are unavailable; commercially sourced normal RNA differs in procurement. None of these GSE91035 contrasts is a definitive clinical-specificity test.

GSE179248 contains integer featureCounts read counts. Unambiguous ENSEMBL-to-symbol mappings are summed by symbol, filtered with edgeR filterByExpr by day group, TMM normalized and modeled with voom/limma, donor + day. The 14 paired donors contribute 28 specimens, not 28 independent patients. Culture day, matrix exposure, isolation and other culture changes are inseparable from the observed ADM-associated response; this is not longitudinal human tumor development.

All eligible genes receive BH correction separately within each study/contrast. A second, explicitly labeled BH family covers the measurable subset of the frozen 506 genes in each external contrast. Same-direction counts and Spearman correlations are descriptive; genes are not independent biological replicates. No pooled cross-platform effect, pooled CP conclusion, or binomial significance test of gene concordance is claimed. Pointwise 95% CIs are not simultaneous confidence bands.

### Serum Models

The fixed family is the 27 unique miRNAs in the historical six, published three, and existing exploratory EV candidate20 pool. The candidate pool was not FDR-significant in EV discovery and is not relabeled as a validated panel. Age/sex-adjusted limma compares PC against healthy, mixed benign, and pooled other cancers within each study. All eligible annotated features receive contrast-level BH; the frozen family also receives BH within each contrast and across all three contrasts per dataset.

GSE59856 has unknown age for one PC case and all 52 liver-cancer controls: adjusted contrasts use 99 PC and 150 healthy/21 benign/248 other-cancer controls. Unadjusted single-marker ROC curves retain all otherwise eligible samples. ROC score direction is locked to the previously analyzed diagnostic-plasma mean difference, never selected to maximize serum AUC. Different analytes and scales prevent direct plasma-coefficient transfer. AUC CIs are pointwise DeLong intervals, with multiple subgroup/marker comparisons retained as exploratory rather than confirmatory wins.

The GSE59856 subgroup uses 17 deposited pStage IIA/IIB cases and excludes ypStage codes; this does not guarantee a complete untreated history or establish screening performance. GSE85589's CA19-9-low subgroup uses 23 PC cases and restricts controls to the same <37 category (19 healthy and 7 cholelithiasis). Stage is unavailable in GSE85589.

GSE85589 deposits CA19-9 only as <37/>37/NA. Two benign participants lack CA19-9 and are excluded from incremental modeling, leaving 88 PC and 27 non-cancer controls. Three ridge logistic models (glmnet alpha=0) use baseline age, sex and binary CA19-9, then add the fixed published three or historical six. Shared stratified 10 outer folds estimate predictions; five inner folds choose lambda by binomial deviance from a fixed 60-value grid. Predictor standardization is fit by glmnet within training folds. No marker selection, serum-ROC sign tuning, or PLCO retuning is performed.

The single deterministic outer split yields conditional out-of-fold-score CIs and DeLong comparisons; they omit pipeline/split uncertainty. The deposited arrays were normalized by their authors before this analysis; raw-array normalization was not repeated within CV folds. These CV estimates are therefore conditional on the supplied processed assay and are not end-to-end validation of a deployable preprocessing pipeline. Brier scores describe a case-control sample, not screening-risk calibration. These are newly fitted serum models with internal validation, not external validation of the original plasma coefficients, and cannot establish superiority over optimally modeled continuous CA19-9.

## Tissue Results

First, genome-wide/eligible-gene counts (not restricted to the spatial family):

| Dataset | Comparison | Cases | Controls | Genes_Tested | All_Gene_BH_Significant |
| --- | --- | --- | --- | --- | --- |
| GSE15471 | PDAC_vs_paired_Normal | 35 | 35 | 20823 | 15541 |
| GSE91035 | PDAC_vs_CP | 25 |  2 | 19100 |    46 |
| GSE91035 | PDAC_vs_normal | 25 |  8 | 19100 |  8877 |
| GSE91035 | PDAC_vs_benign | 25 | 15 | 19100 |  4597 |
| GSE179248 | Culture_D6_vs_D0 | 14 | 14 | 16466 | 10519 |

The next table counts external all-gene BH significance only among the frozen measurable spatial genes. The dual column additionally requires the matched spatial contrast to pass its original all-target BH threshold and the external effect to have the same sign. CP and benign comparisons are listed for transparency but do not reproduce a Normal-to-PDAC comparison directly.

| Dataset | Contrast | Frozen_Measurable | Same_Direction | External_AllGene_BH_Significant | Both_Significant_Same_Direction | Spearman |
| --- | --- | --- | --- | --- | --- | --- |
| GSE15471 | PDAC_vs_paired_Normal | 499 | 345 | 355 | 272 | 0.812 |
| GSE179248 | Culture_D6_vs_D0 | 484 | 357 | 351 |  35 | 0.685 |
| GSE91035 | PDAC_vs_CP | 469 | 280 |   0 |   0 | 0.336 |
| GSE91035 | PDAC_vs_benign | 469 | 379 | 170 | 160 | 0.822 |
| GSE91035 | PDAC_vs_normal | 469 | 382 | 294 | 268 | 0.835 |

Of the 53 measurable spatial Normal-to-ADM-significant genes in the culture dataset, 43 have the same direction and 35 additionally pass external all-gene BH <0.05. In the paired tumor cohort, 272 frozen genes are significant in both matched contrasts with concordant direction. The two-CP exploratory comparison has 46 all-gene BH-significant associations, but none are in the frozen spatial family; zero frozen-family BH-significant genes are observed. Its full results are retained. Two controls and missing clinical adjustment prevent a robust specificity claim; failure of the frozen family to reach significance is not evidence of equivalence.

### STAT1

| Dataset | Contrast | log2FC_95CI | All_Gene_BH |
| --- | --- | --- | --- |
| GSE15471 | PDAC_vs_paired_Normal | 0.940 (0.553 to 1.326) | 5.77e-05 |
| GSE179248 | Culture_D6_vs_D0 | 0.585 (0.400 to 0.771) | 1.46e-05 |

The independent culture increase is compatible with the spatial ADM-versus-Normal increase. The paired tumor increase supports a tumor/adjacent-normal association. Neither establishes the historical claim of a significant transient ADM peak: the original spatial PDAC-versus-ADM decrease remains non-significant, and the original age/sex-adjusted PDAC-versus-CP contrast remains non-significant.

## Serum Results

GSE59856 has 20/27 complete, unambiguously mapped frozen markers; the excluded markers are miR-130a-3p, miR-134-3p, miR-139-3p, miR-20a-5p, miR-210-3p, miR-2110, miR-23c. GSE85589 has 27/27. Missing or ambiguous measurement is not biological non-replication.

The following significance count uses the frozen-family BH correction across all three comparisons per dataset. Direction is compared with diagnostic plasma. Adjusted association and unadjusted single-marker discrimination are separate analyses.

| Dataset | Contrast | Frozen_Measurable | Adjusted_PC | Adjusted_Control | Family_BH_Significant | Significant_Same_Direction |
| --- | --- | --- | --- | --- | --- | --- |
| GSE59856 | PC_vs_Benign | 20 | 99 |  21 | 3 | 2 |
| GSE59856 | PC_vs_Healthy | 20 | 99 | 150 | 3 | 1 |
| GSE59856 | PC_vs_Other_cancer | 20 | 99 | 248 | 6 | 4 |
| GSE85589 | PC_vs_Benign | 27 | 88 |  10 | 0 | 0 |
| GSE85589 | PC_vs_Healthy | 27 | 88 |  19 | 0 | 0 |
| GSE85589 | PC_vs_Other_cancer | 27 | 88 | 115 | 5 | 3 |

### Binary CA19-9 Increment

| Model | CV_AUC_95CI | Brier |
| --- | --- | --- |
| CA19_9_age_sex | 0.893 (0.828 to 0.959) | 0.106 |
| CA19_9_plus_published3 | 0.900 (0.829 to 0.971) | 0.102 |
| CA19_9_plus_historical6 | 0.887 (0.829 to 0.946) | 0.121 |

Neither panel improves AUC convincingly: paired conditional DeLong BH q-values are 0.739 and 0.739 for the published-three and historical-six additions. The small control sample and lack of continuous CA19-9 limit this analysis.

## Figures

Each figure has PNG and PDF versions in reanalysis/figures/external/.

1. tissue_direction_concordance: all frozen genes measurable in paired tumor, commercial-normal, and paired culture comparisons; signs are descriptive.
2. STAT1_external_forest: pointwise CIs across distinct spatial, paired-tissue, culture, and adjusted-CP contrasts; no pooled estimate.
3. serum_fixed_marker_AUC: all nine historical/published markers, fixed plasma directions, healthy/benign/other-cancer comparisons; missing measurements remain blank.
4. serum_CA19_9_nested_CV: shared-fold out-of-fold ROC curves for baseline and two fixed-panel additions.

## Publication Interpretation

The added data strengthen a focused cross-cohort tissue-state reproducibility analysis and provide a negative clinical-increment benchmark in serum. They do not rescue the original early-detection or malignant-switch claims. A journal-facing paper should distinguish paired culture response, cross-sectional tissue association, diagnostic case-control discrimination, and prediagnostic transportability throughout.

Substantive gaps remain: confirm spatial patient IDs, obtain an adequately sized independent and batch-balanced cancer-versus-inflammatory tissue cohort, audit cell composition, and resolve sample-level PLCO lead time/clinical information. Public-only analysis does not require inventing lab validation, but neither dataset accumulation nor attractive ROC curves guarantees novelty or journal acceptance. The full manuscript still requires author review, declarations and a focused prior-art assessment.

## Reproduction and Sources

Run scripts/reanalysis/run_external.R after the original reanalysis. Dependencies are pinned in renv.lock, source hashes are saved, and scripts/reanalysis/17_validate_external.R checks sample pairing, BH recomputation, mapping, AUC/Brier recomputation, subgroup sizes and source hashes. Four figure pairs are additionally rendered for visual QA.

- [GSE15471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471) and [published-QC curation](https://master.bioconductor.org/packages/devel/data/experiment/manuals/KEGGdzPathwaysGEO/man/KEGGdzPathwaysGEO.pdf).
- [GSE179248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179248); [human ADM culture source article](https://doi.org/10.1016/j.gastha.2023.02.003).
- [GSE91035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91035) and [GPL22763](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL22763).
- [GSE59856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59856); [serum source article](https://doi.org/10.1371/journal.pone.0118220).
- [GSE85589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85589) and [GPL18941 mature-miRNA annotation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18941).
- [GSE101462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101462), excluded from inferential validation for perfect chip/disease confounding.
- [glmnet methods and implementation](https://glmnet.stanford.edu/articles/glmnet.html).
