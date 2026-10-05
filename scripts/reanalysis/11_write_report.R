#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
read_result <- function(path) read.csv(file.path("reanalysis", "results", path))
bulk <- read_result("bulk/GSE143754_gene_level_contrasts.csv")
spatial <- read_result("spatial/GSE208536_patient_blocked_all_targets.csv")
historical80 <- read_result("spatial/GSE208536_historical80_patient_aware.csv")
complete <- read_result("spatial/GSE208536_complete_profiles_all_targets.csv")
loo <- read_result("spatial/GSE208536_leave_one_profile_out.csv")
integrated <- read_result("integrated_tissue/analysis_flow_counts.csv")
ev <- read_result("ev/GSE304572_EV_differential.csv")
voom <- read_result("ev/GSE304572_voom_sensitivity.csv")
plasma <- read_result("plasma/diagnostic_to_PLCO_performance.csv")
comparison <- read_result("plasma/PLCO_paired_model_comparison.csv")
go <- read_result("enrichment/GO_BP_jointly_measurable_background.csv")
ages <- read_result("qc/bulk_age_summary.csv")
stat1 <- spatial[spatial$Gene == "STAT1", ]
stat1_complete <- complete[complete$Gene == "STAT1", ]
stat1_loo <- loo[loo$Gene == "STAT1", ]
md_table <- function(data) {
  header <- paste0("| ", paste(names(data), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |")
  rows <- apply(data, 1, function(row) {
    paste0("| ", paste(gsub("|", "/", row, fixed = TRUE), collapse = " | "), " |")
  })
  c(header, separator, rows)
}
performance <- data.frame(
  Model = gsub("_", " ", plasma$Model),
  Apparent_AUC = sprintf("%.3f", plasma$Apparent_AUC),
  Diagnostic_CV_AUC = sprintf("%.3f", plasma$Diagnostic_CV_AUC),
  PLCO_AUC_95CI = sprintf("%.3f (%.3f-%.3f)", plasma$PLCO_AUC,
                          plasma$PLCO_AUC_CI_low, plasma$PLCO_AUC_CI_high),
  PLCO_Brier = sprintf("%.3f", plasma$PLCO_Brier)
)
calibration <- data.frame(
  Model = gsub("_", " ", plasma$Model),
  Training_Threshold = sprintf("%.3f", plasma$Training_Youden_Threshold),
  PLCO_Sensitivity_95CI = sprintf("%.3f (%.3f-%.3f)", plasma$PLCO_Sensitivity,
                                  plasma$PLCO_Sensitivity_CI_low, plasma$PLCO_Sensitivity_CI_high),
  PLCO_Specificity_95CI = sprintf("%.3f (%.3f-%.3f)", plasma$PLCO_Specificity,
                                  plasma$PLCO_Specificity_CI_low, plasma$PLCO_Specificity_CI_high),
  Calibration_Intercept = sprintf("%.3f", plasma$PLCO_Calibration_Intercept),
  Calibration_Slope = sprintf("%.3f", plasma$PLCO_Calibration_Slope)
)
report <- c(
  "# Corrected Analysis: Methods and Results",
  "",
  "Public-data analysis reviewed on 2026-10-04. This report is generated from the final tables by scripts/reanalysis/11_write_report.R. It supplies replacement methods/results material and an evidence assessment; it is not an author-approved submission manuscript.",
  "",
  "## Working Title",
  "",
  "Patient-aware tissue-state associations and prediagnostic transportability of plasma miRNA models in public PDAC cohorts",
  "",
  "## Research Questions",
  "",
  "We assessed whether tissue-state associations persist when spatial ROIs are treated as repeated observations, whether PDAC-versus-chronic-pancreatitis tissue contrasts remain after clinical adjustment, and whether EV-ranked miRNA selection improves diagnostic-to-prediagnostic transportability relative to fixed published/historical panels and unrestricted plasma selection. These are exploratory secondary analyses of previously published data.",
  "",
  "## Cohort Audit",
  "",
  "| Dataset | Deposited analysis sample | Role |",
  "| --- | --- | --- |",
  "| GSE143754 | 9 adjacent normal, 6 CP, 11 PDAC tissues | Bulk tissue comparison |",
  "| GSE208536 | 48 ROIs: 16 Normal, 16 ADM, 16 PDAC; eight inferred profiles | Spatial tissue-state analysis |",
  "| GSE304572 | 65 PDAC, 10 CP, 10 IPMN plasma EV samples | Independent exploratory miRNA ranking |",
  "| GSE259327 diagnostic | 121 PDAC, 82 pooled controls | Plasma model training |",
  "| GSE259327 PLCO | 48 prediagnostic PDAC, 48 controls | Same-assay external evaluation |",
  "| GSE268771 | 51 PDAC, 12 benign, 3 healthy plasma EV samples | Separate descriptive sensitivity cohort |",
  "",
  "The GSE259327 source studies report 46 healthy and 36 pancreatitis controls; chronic pancreatitis is a subset rather than a synonym for all 36 pancreatitis participants. The deposited sample records label controls only as Control, so subtype-specific ROC estimates cannot be reconstructed. The secondary EV PDAC stages are I/II/III/IV = 2/4/12/33. Spatial TMA identifiers were matched exactly to all 48 workbook columns. Age, sex, grade, and stage combinations identify eight distinct clinical profiles, consistent with the source design, but GEO does not supply explicit patient IDs. ROI counts are uneven; six profiles contain all three states, one lacks ADM, and one lacks Normal.",
  "",
  "## Methods",
  "",
  "### Bulk Tissue",
  "",
  "We used the deposited invariant-set-normalized, log-scale GSE143754 matrix rather than repeating raw CEL preprocessing. Transcript-cluster identifiers were mapped with hta20transcriptcluster.db 8.8.0. Ambiguous symbol mappings and unannotated identifiers were excluded. For multiple probes mapped to one symbol, the probe with the greatest mean expression across all samples was retained, with probe ID as the tie breaker; this rule does not use outcome labels.",
  "",
  "The primary limma model was restricted to 11 PDAC and six CP specimens, with a linear age covariate and sex. This excludes possibly paired adjacent-normal specimens because explicit patient pairing is unavailable. Contrast standard errors were moderated with empirical Bayes, trend = TRUE and robust = TRUE. Unadjusted PDAC-versus-CP estimates used the same 17 specimens. A separate three-group adjusted model supplied exploratory PDAC-versus-normal and CP-versus-normal contrasts, plus an all-group PDAC-versus-CP sensitivity result; its independence assumption cannot be verified without pairing metadata. BH correction was applied across all analyzed genes separately for each contrast. FDR < 0.05 and absolute log2 difference >= 0.5 defined the operational CP-comparator candidate rule. No set subtraction based on non-significance was used, and inflammation-free or malignancy-specific expression was not inferred.",
  "",
  "### Spatial Tissue States",
  "",
  "Deposited normalized expression was transformed as log2(expression + 1). ROIs were averaged within each inferred patient profile and tissue state, producing 22 observed profile/state means. We fitted expression ~ profile + tissue state with a fixed profile block. The state omnibus test and Normal-to-ADM, ADM-to-PDAC, and Normal-to-PDAC contrasts were corrected separately across all 1,825 measured targets. Contrast CIs are pointwise 95% intervals, not multiplicity-adjusted intervals.",
  "",
  "The analysis was repeated in the six profiles containing all states and after leaving each profile out. A ROI-level random-intercept nlme model was fitted for the historical 80-gene family as a sensitivity analysis, with BH correction across that 80-gene family. Its results are not numerically comparable to all-panel FDR counts because both weighting and the tested family differ. Tissue-state labels describe cross-sectional regions, not observed longitudinal progression.",
  "",
  "### EV miRNA Ranking",
  "",
  "GSE304572 count columns were mapped exactly to the deposited participant labels. Counts were normalized to miRNA library totals and transformed as log2(CPM + 1). miRNAs required at least 1 CPM in at least ten samples. PDAC was compared with pooled CP/IPMN using a two-sided Wilcoxon test with exact = FALSE, followed by BH correction across detectable miRNAs. A count-aware limma-voom model using the same features and library sizes was included as sensitivity evidence. No TMM normalization or raw-read reprocessing was performed.",
  "",
  "The top 20 eligible measured plasma miRNAs were ranked using EV BH FDR and then absolute median difference. This creates an exploratory candidate pool despite the absence of significant EV hits. It is not a validated EV signature and does not directly link a tissue gene to a plasma miRNA.",
  "",
  "### Diagnostic Training and PLCO Evaluation",
  "",
  "Both GSE259327 workbooks were transformed as log2(raw count / deposited Total Counts * 1,000,000 + 1). Sample titles were matched to GEO labels; all 299 participants and the 2,102 shared assay features were reconciled. Four unpenalized logistic models were assessed: the historical six-miRNA panel, a refit of the published three-marker set (let-7i-5p, miR-130a-3p, miR-221-3p), EV-ranked selection, and unrestricted measured plasma-miRNA selection.",
  "",
  "The fixed panels received stratified ten-fold CV with seed 20261001. The two selection pipelines used stratified ten-fold outer CV with seed 20261002 and five-fold inner CV to choose 3, 6, or 9 features. Feature ranking by Wilcoxon P was repeated within every training fold. The unrestricted pipeline also repeated detection filtering within training folds: at least 1 CPM in max(5, ceiling(10% of training samples)). The external EV pool was formed from a separate cohort. Final feature counts were chosen by inner CV on the full diagnostic cohort with seed 20261025. All fits were required to converge with finite coefficients.",
  "",
  "Final models were fitted on all diagnostic samples and applied to PLCO without using PLCO outcomes for feature selection, coefficient fitting, threshold tuning, direction reversal, or outcome-driven batch correction. ROC direction was fixed so larger scores indicate PDAC. PLCO AUCs received DeLong 95% CIs, and model comparisons used paired DeLong tests with BH correction across six comparisons. Ten-fold out-of-fold AUC intervals are conditional DeLong intervals; they do not capture full pipeline/split uncertainty. The historical six-panel CV also excludes its earlier selection process.",
  "",
  "Youden thresholds were chosen from apparent diagnostic ROC curves and retained in PLCO. Sensitivity/specificity CIs used exact binomial intervals. Brier score and logistic calibration intercept/slope were estimated in the balanced PLCO sample; probabilities were clipped to 1e-6 to 1 - 1e-6 only for calibration logits. These measures are sample-specific and do not establish population screening-risk calibration.",
  "",
  "### Secondary Cohort and Enrichment",
  "",
  "GSE268771 deposited RPM values were transformed as log2(RPM + 1). PDAC was compared with all controls and with benign controls using separate Wilcoxon tests with BH correction over 1,175 miRNAs. Constant features received P = 1. Marker availability and fixed-direction univariate AUCs were descriptive; diagnostic plasma model coefficients were not transferred to this different EV assay.",
  "",
  "GO Biological Process enrichment used the spatial-associated genes measurable in both tissue datasets. The background was restricted to jointly measurable genes with GO BP annotation, rather than the whole genome. GOALL propagated annotations were used; terms with 10-500 background genes were tested by an upper-tail hypergeometric test and BH correction. The exact universe is deposited. This is supplementary context for a cancer-focused panel, without a claim of CP-specific pathway enrichment.",
  "",
  "## Results",
  "",
  "### Bulk Adjustment and Tissue Comparison",
  "",
  sprintf("The updated platform mapping retained %s unique genes. The unadjusted PDAC-versus-CP analysis yielded %s genes at FDR < 0.05 and %s after the absolute log2 difference >= 0.5 rule. After age/sex adjustment in the 17-specimen primary model, %s genes passed FDR < 0.05; the exploratory all-group adjusted model yielded %s. The adjusted CP-aware spatial overlap contains %s genes. Unadjusted overlap counts are explicitly exploratory and cannot support the previous validated malignant-switch claim.",
          nrow(bulk), sum(bulk$Unadjusted_PDAC_vs_CP_BH_FDR < 0.05),
          sum(bulk$Unadjusted_PDAC_vs_CP_BH_FDR < 0.05 &
                abs(bulk$Unadjusted_PDAC_vs_CP_log2FC) >= 0.5),
          sum(bulk$PDAC_vs_CP_BH_FDR < 0.05),
          sum(bulk$AllGroup_Adjusted_PDAC_vs_CP_BH_FDR < 0.05),
          integrated$Count[integrated$Metric == "Age-sex-adjusted CP-aware spatial overlap"]),
  "",
  "Bulk age distributions:",
  "",
  md_table(ages),
  "",
  "Tissue analysis flow:",
  "",
  md_table(integrated),
  "",
  "### Patient-Profile-Aware Spatial Associations",
  "",
  sprintf("The primary all-panel test identified %s of %s targets at omnibus FDR < 0.05. Restricting to the six complete profiles yielded %s targets. In the historical 80-gene set, %s passed the primary all-target FDR and %s passed the ROI random-intercept sensitivity FDR over 80 tests. These analyses support state associations while showing dependence on analysis scale and multiplicity.",
          sum(spatial$Patient_Blocked_F_BH_FDR_all_targets < 0.05), nrow(spatial),
          sum(complete$Patient_Blocked_F_BH_FDR_all_targets < 0.05),
          sum(historical80$Patient_Blocked_F_BH_FDR_all_targets < 0.05),
          sum(historical80$ROI_Mixed_Model_BH_FDR_80 < 0.05, na.rm = TRUE)),
  "",
  sprintf("STAT1 increased from Normal to ADM by %.3f log2 units (95%% CI %.3f to %.3f; all-target contrast FDR %.4f). Its ADM-to-PDAC difference was %.3f (95%% CI %.3f to %.3f; FDR %.4f), which does not establish a subsequent decline. The omnibus FDR was %.4f; the complete-profile omnibus FDR was %.4f. Leaving one profile out preserved the positive early-effect direction in all eight fits (effect range %.3f-%.3f), but the omnibus FDR ranged from %.4f to %.4f. STAT1 is therefore an early-associated example with finite-sample uncertainty, not a proven transient peak.",
          stat1$Normal_to_ADM_log2, stat1$Normal_to_ADM_CI_low, stat1$Normal_to_ADM_CI_high,
          stat1$Normal_to_ADM_BH_FDR_all_targets, stat1$ADM_to_PDAC_log2,
          stat1$ADM_to_PDAC_CI_low, stat1$ADM_to_PDAC_CI_high,
          stat1$ADM_to_PDAC_BH_FDR_all_targets, stat1$Patient_Blocked_F_BH_FDR_all_targets,
          stat1_complete$Patient_Blocked_F_BH_FDR_all_targets,
          min(stat1_loo$Normal_to_ADM_log2), max(stat1_loo$Normal_to_ADM_log2),
          min(stat1_loo$Patient_Blocked_F_BH_FDR_all_targets),
          max(stat1_loo$Patient_Blocked_F_BH_FDR_all_targets)),
  "",
  "### EV Discovery",
  "",
  sprintf("Of 2,369 measured miRNAs, %s met the detection filter; %s passed Wilcoxon BH FDR < 0.05 and %s passed voom BH FDR < 0.05. There were 363 detectable features measured in both plasma workbooks, from which the exploratory 20-feature EV ranking was formed. Negative results under this filter do not invalidate source-paper analyses of rare miRNAs or radiomics-derived signatures.",
          nrow(ev), sum(ev$BH_FDR < 0.05), sum(voom$adj.P.Val < 0.05)),
  "",
  "### Plasma Model Transportability",
  "",
  md_table(performance),
  "",
  "Diagnostic CV performance remained high, while every PLCO AUC CI included 0.5. All PLCO samples were below their model's diagnostic Youden threshold, producing zero sensitivity and unit specificity. This accompanies a large cohort-associated score shift; the files cannot separate lead-time biology from preanalytical or normalization differences.",
  "",
  md_table(calibration),
  "",
  sprintf("No pairwise PLCO AUC comparison survived BH correction (minimum adjusted P = %.3f). The EV-ranked strategy did not establish improvement over unrestricted selection or the refitted published marker set. The historical-six versus published-three-refit unadjusted paired DeLong P was %.4f. The full comparison table and model coefficients are deposited.",
          min(comparison$BH_FDR),
          comparison$Paired_DeLong_p[comparison$Model_1 == "Historical_six" &
                                      comparison$Model_2 == "Published_three_refit"]),
  "",
  "The source paper's published lead-time-specific AUC must not be compared as if it were obtained with these refitted coefficients, this normalization, and all 96 PLCO samples. Original score reconstruction and near-diagnosis subsets cannot be recovered from the current deposited clinical metadata. PLCO has now been examined, so future model redesign would require a fresh independent test set.",
  "",
  "### Secondary EV Cohort and Panel-Aware Enrichment",
  "",
  "GSE268771 contained eight of the nine historical/published miRNA markers, with 236 miRNAs passing PDAC-versus-benign FDR < 0.05. Forty-five of 51 PDAC participants were stage III/IV, and 32 were marked as receiving palliative treatment. These results are descriptive evidence in an advanced-disease cohort.",
  "",
  sprintf("GO BP enrichment tested %s terms with an annotated jointly measurable background of %s genes and %s selected spatial-associated genes. %s terms passed BH FDR < 0.05. The result concerns tissue-state association within the measurable panel and is not CP-specific enrichment.",
          nrow(go), unique(go$Background_Total), unique(go$Selected_Total),
          sum(go$BH_FDR < 0.05)),
  "",
  "## Interpretation",
  "",
  "The strongest retained observation is cross-sectional spatial tissue-state association after accounting for repeated ROIs. The current data do not establish a confounder-robust CP-specific gene set or a transferable prediagnostic miRNA panel. A credible manuscript should report the adjustment sensitivity and negative external evaluation directly, with historical network/database findings moved to supporting material. Publication readiness still depends on confirming spatial identities, resolving source-score/preprocessing differences, and establishing a focused contribution beyond the primary studies.",
  "",
  "The later [publication audit](publication_audit_methods_results.md) corroborates the eight inferred clinical profiles against the source article and reports composition sensitivities and the status of a larger inflammatory comparison. Read that follow-up alongside this original report; the direct patient identifiers and individual clinical-data limitations are not silently resolved.",
  "",
  "Detailed claim changes and remaining constraints are in [reanalysis_claims_and_limitations.md](reanalysis_claims_and_limitations.md). Original dataset publications and accession links are in [the reanalysis guide](../reanalysis/README.md). All source studies must be cited in the manuscript.",
  "",
  "## Figures and Captions",
  "",
  "Each listed figure has a same-stem PNG and PDF under reanalysis/figures.",
  "",
  "| Figure | Caption |",
  "| --- | --- |",
  "| spatial/STAT1_patient_stage | Mean log2 normalized STAT1 expression within inferred profile/state; lines join available states within profiles. |",
  "| bulk/PDAC_vs_CP_volcano | Age/sex-adjusted limma PDAC-versus-CP contrast; the dashed line is BH FDR = 0.05. |",
  "| integrated_tissue/spatial_tissue_state_heatmap | Top 20 jointly measurable spatial-associated genes, ranked by all-panel omnibus FDR; colors are within-gene z-scores of profile/state means. |",
  "| ev/EV_differential | Detectable EV miRNAs under the stated filter; Wilcoxon median differences and BH FDR. |",
  "| plasma/PLCO_ROC_models | Fixed-direction ROC curves for diagnostic-fitted models on all 96 deposited PLCO participants. |",
  "| plasma/score_distributions | Fitted model score distributions across diagnostic and PLCO cohorts, stratified by deposited outcome. |",
  "| secondary_ev/cohort_composition | Actual comparator and stage counts in GSE268771. |",
  "| qc/plasma_total_counts | Deposited plasma total-count distributions by cohort and outcome. |",
  "| qc/PLCO_calibration | Six equal-size score-rank bins per model in the balanced PLCO sample; exact binomial CIs for observed case fractions. |",
  "| qc/bulk_demographics | Deposited participant ages by tissue group and sex. |",
  "",
  "## Reproducibility",
  "",
  "The scripts, seed values, input MD5 hashes, package versions, sessionInfo, model coefficients, per-sample scores, fold selections, and PNG/PDF pairs are deposited. See scripts/reanalysis/README.md for execution. This workflow uses processed public data; it is not a reproduction of each primary study's raw-data processing. Output checks confirm structural/numerical validity, not biological truth or external validation of the workflow itself."
)
dir.create("docs", showWarnings = FALSE)
writeLines(report, "docs/reanalysis_methods_results.md")
cat("Methods/results report generated from final output tables.\n")
