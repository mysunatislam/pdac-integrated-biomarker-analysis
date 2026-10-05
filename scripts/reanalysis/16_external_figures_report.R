#!/usr/bin/env Rscript
source("scripts/reanalysis/external_helpers.R")
fig_dir <- "reanalysis/figures/external"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
export <- function(plot, name, width, height) {
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = width, height = height)
}
tissue <- read.csv(file.path(external_dir, "tissue_frozen_spatial_family.csv"))
summary <- read.csv(file.path(external_dir, "tissue_replication_summary.csv"))
subset <- tissue[tissue$Comparison %in% c("PDAC_vs_paired_Normal", "PDAC_vs_normal", "Culture_D6_vs_D0"), ]
subset$Cohort <- factor(paste(subset$Dataset, subset$Comparison),
  levels = c("GSE15471 PDAC_vs_paired_Normal", "GSE91035 PDAC_vs_normal", "GSE179248 Culture_D6_vs_D0"),
  labels = c("GSE15471: 35 donor pairs\nPDAC - adjacent normal",
             "GSE91035: 25 PDAC / 8 normal\nCommercial normal RNA (exploratory)",
             "GSE179248: 14 donor pairs\nCulture day 6 - day 0"))
scatter <- ggplot2::ggplot(subset, ggplot2::aes(Spatial_log2FC, log2FC)) +
  ggplot2::geom_hline(yintercept = 0, color = "#AAB1B8", linewidth = .3) +
  ggplot2::geom_vline(xintercept = 0, color = "#AAB1B8", linewidth = .3) +
  ggplot2::geom_point(ggplot2::aes(color = Same_direction), alpha = .65, size = 1.2) +
  ggplot2::facet_wrap(~Cohort, nrow = 1, scales = "free") +
  ggplot2::scale_color_manual(values = c("TRUE" = "#168273", "FALSE" = "#C04B58"),
                              labels = c("FALSE" = "Opposite direction", "TRUE" = "Same direction")) +
  ggplot2::labs(x = "Spatial contrast (log2 expression difference)",
                y = "External contrast (log2 expression difference)", color = NULL,
                title = "Frozen spatial associations in additional human datasets",
                subtitle = "Normal-to-PDAC reference for bulk tissue; Normal-to-ADM reference for culture") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(size = 9))
export(scatter, "tissue_direction_concordance", 10, 4)

spatial <- read.csv("reanalysis/results/spatial/GSE208536_patient_blocked_all_targets.csv")
spatial <- spatial[spatial$Gene == "STAT1", ]
bulk <- read.csv("reanalysis/results/bulk/GSE143754_gene_level_contrasts.csv")
bulk <- bulk[bulk$Gene == "STAT1", ]
external <- read.csv(file.path(external_dir, "tissue_all_gene_contrasts.csv"))
external <- external[external$Gene == "STAT1" & external$Dataset %in% c("GSE15471", "GSE179248"), ]
forest <- rbind(
  data.frame(Label = c("Spatial: ADM - Normal (8 inferred profiles)",
                       "Spatial: PDAC - ADM (8 inferred profiles)"),
             Effect = c(spatial$Normal_to_ADM_log2, spatial$ADM_to_PDAC_log2),
             Low = c(spatial$Normal_to_ADM_CI_low, spatial$ADM_to_PDAC_CI_low),
             High = c(spatial$Normal_to_ADM_CI_high, spatial$ADM_to_PDAC_CI_high), Role = "Spatial"),
  data.frame(Label = ifelse(external$Dataset == "GSE15471",
                            "GSE15471: PDAC - normal (35 donor pairs)",
                            "GSE179248: culture D6 - D0 (14 donor pairs)"),
             Effect = external$log2FC, Low = external$CI_low, High = external$CI_high, Role = "Additional cohort"),
  data.frame(Label = "GSE143754: PDAC - CP (age/sex adjusted)",
             Effect = bulk$PDAC_vs_CP_log2FC, Low = bulk$PDAC_vs_CP_CI_low,
             High = bulk$PDAC_vs_CP_CI_high, Role = "CP comparator")
)
forest$Label <- factor(forest$Label, levels = rev(forest$Label))
p <- ggplot2::ggplot(forest, ggplot2::aes(Effect, Label, color = Role)) +
  ggplot2::geom_vline(xintercept = 0, color = "#777777", linetype = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = Low, xmax = High), orientation = "y", width = .18) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_color_manual(values = c("Spatial" = "#168273", "Additional cohort" = "#C04B58", "CP comparator" = "#656A75")) +
  ggplot2::labs(x = "STAT1 log2 expression difference (pointwise 95% CI)", y = NULL,
                title = "STAT1 across distinct experimental comparisons",
                subtitle = "Different contrasts and assays; estimates are not pooled") +
  ggplot2::theme_classic(base_size = 10) + ggplot2::theme(legend.position = "none")
export(p, "STAT1_external_forest", 9.5, 4)

aucs <- read.csv(file.path(external_dir, "serum_fixed_direction_marker_AUC.csv"))
marker_reference <- read.csv(file.path(external_dir, "frozen_miRNA_reference.csv"))
nine <- marker_reference$miRNA[marker_reference$Historical6 | marker_reference$Published3]
aucs <- aucs[aucs$miRNA %in% nine & aucs$Subgroup == "All", ]
aucs$miRNA <- factor(aucs$miRNA, levels = rev(nine))
aucs$Comparison <- factor(aucs$Comparison, levels = c("PC_vs_Healthy", "PC_vs_Benign", "PC_vs_Other_cancer"),
                          labels = c("PC vs healthy", "PC vs benign disease", "PC vs other cancers"))
p <- ggplot2::ggplot(aucs, ggplot2::aes(AUC, miRNA)) +
  ggplot2::geom_vline(xintercept = .5, color = "#777777", linetype = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .15, color = "#7B818B") +
  ggplot2::geom_point(size = 1.8, color = "#168273") +
  ggplot2::facet_grid(Dataset ~ Comparison) +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::labs(x = "Fixed-direction single-marker AUC (pointwise 95% CI)", y = NULL,
                title = "Historical and published miRNA markers in serum",
                subtitle = "Direction fixed from diagnostic plasma; blank rows are not measurable with the annotation/filter") +
  ggplot2::theme_classic(base_size = 10) + ggplot2::theme(strip.text = ggplot2::element_text(size = 9))
export(p, "serum_fixed_marker_AUC", 9.5, 6.8)

scores <- read.csv(file.path(external_dir, "GSE85589_CA19_9_out_of_fold_scores.csv"))
cv <- read.csv(file.path(external_dir, "GSE85589_CA19_9_nested_CV_performance.csv"))
curves <- do.call(rbind, lapply(cv$Model, function(model) {
  r <- pROC::roc(scores$Case, scores[[model]], direction = "<", quiet = TRUE)
  row <- cv[cv$Model == model, ]
  data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities,
              Model = sprintf("%s (AUC %.3f)", gsub("_", " ", model), row$AUC))
}))
p <- ggplot2::ggplot(curves, ggplot2::aes(FPR, TPR, color = Model)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#777777") +
  ggplot2::geom_path(linewidth = .8) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::scale_color_manual(values = c("#656A75", "#C04B58", "#168273")) +
  ggplot2::labs(x = "False-positive rate", y = "True-positive rate", color = NULL,
                title = "Serum miRNAs added to binary CA19-9, age and sex",
                subtitle = "GSE85589: 88 PC / 27 controls; nested 10-by-5-fold CV") +
  ggplot2::theme_classic(base_size = 10) + ggplot2::theme(legend.position = "bottom") +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))
export(p, "serum_CA19_9_nested_CV", 6.5, 6)
cat("Four external PNG/PDF pairs exported.\n")

md_table <- function(d) {
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
    apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |")))
}
models <- read.csv(file.path(external_dir, "tissue_models.csv"))
all_tissue <- read.csv(file.path(external_dir, "tissue_all_gene_contrasts.csv"))
stat1 <- all_tissue[all_tissue$Gene == "STAT1", ]
cp_all <- all_tissue[all_tissue$Dataset == "GSE91035" & all_tissue$Comparison == "PDAC_vs_CP", ]
model_table <- models[, c("Dataset", "Comparison", "Cases", "Controls", "Genes_Tested", "All_Gene_BH_Significant")]
stat1_table <- data.frame(Dataset = stat1$Dataset, Contrast = stat1$Comparison,
                          log2FC_95CI = sprintf("%.3f (%.3f to %.3f)", stat1$log2FC, stat1$CI_low, stat1$CI_high),
                          All_Gene_BH = sprintf("%.3g", stat1$BH_all_genes))
summary_table <- data.frame(Dataset = summary$Dataset, Contrast = summary$Comparison,
                            Frozen_Measurable = summary$Measurable_Frozen,
                            Same_Direction = summary$Same_Direction,
                            External_AllGene_BH_Significant = summary$External_BH_All_Significant,
                            Both_Significant_Same_Direction = summary$Dual_AllGene_Significant_SameDirection,
                            Spearman = sprintf("%.3f", summary$Spearman_Effect_Correlation))
serum <- read.csv(file.path(external_dir, "serum_frozen_miRNA_family.csv"))
serum_summary <- do.call(rbind, lapply(split(serum, paste(serum$Dataset, serum$Comparison)), function(s) {
  data.frame(Dataset = s$Dataset[1], Contrast = s$Comparison[1], Frozen_Measurable = nrow(s),
             Adjusted_PC = s$Adjusted_Cases[1], Adjusted_Control = s$Adjusted_Controls[1],
             Family_BH_Significant = sum(s$BH_frozen_all_comparisons < .05),
             Significant_Same_Direction = sum(s$BH_frozen_all_comparisons < .05 & s$Same_Direction_As_Plasma))
}))
write.csv(serum_summary, file.path(external_dir, "serum_family_summary.csv"), row.names = FALSE)
cv_table <- data.frame(Model = cv$Model,
                       CV_AUC_95CI = sprintf("%.3f (%.3f to %.3f)", cv$AUC, cv$CI_low, cv$CI_high),
                       Brier = sprintf("%.3f", cv$Brier))
comparison <- read.csv(file.path(external_dir, "GSE85589_CA19_9_conditional_comparison.csv"))
coverage <- read.csv(file.path(external_dir, "serum_marker_coverage.csv"))
missing <- coverage$miRNA[coverage$Dataset == "GSE59856" & !coverage$Measurable_Complete]
report <- c(
  "# Additional Public Cohorts: Methods and Results", "",
  "Generated from the saved tables on 2026-10-04. This extends the corrected analysis; it is not a submission-ready manuscript or evidence of clinical validation.", "",
  "## Main Findings", "",
  "Independent paired human tissue and culture data support parts of the frozen spatial tissue-state pattern, including a STAT1 increase. These are association and culture-response results. They do not demonstrate a cancer-specific malignant switch, a transient ADM peak, tissue-to-miRNA regulation, or prediagnostic biomarker performance.", "",
  "In a separate serum study, the fixed miRNA panels did not materially improve nested-CV discrimination over binary CA19-9, age and sex. The prior diagnostic-to-PLCO results and coefficients were preserved without further PLCO tuning.", "",
  "## Cohort Eligibility", "",
  "| Dataset | Deposited samples | Analysis unit and role |",
  "| --- | --- | --- |",
  "| GSE15471 | 78 arrays; 36 tumor/normal donor pairs and technical replicates | Average technical replicates; exclude published failed-QC donor 51721; analyze 35 donor pairs |",
  "| GSE179248 | 28 RNA-seq samples | 14 paired donors, acinar cells at culture day 0/day 6; culture response accompanying ADM |",
  "| GSE91035 | 50 arrays: 25 PDAC, 2 CP, 15 adjacent benign, 8 commercial normal RNA | Separate exploratory two-group contrasts; source overview's 27 PDAC count does not match deposited disease labels |",
  "| GSE59856 | 571 serum samples: 100 PC, 150 healthy, 21 mixed benign pancreatic/biliary, 300 other cancers | Marker associations; 17 pStage IIA/IIB cases for a diagnostic subgroup |",
  "| GSE85589 | 232 serum samples: 88 PC, 19 healthy, 10 cholelithiasis, 115 other cancers | Marker associations and exploratory binary-CA19-9 increment analysis |",
  "| GSE101462 | 20 tissue arrays: 6 PDAC, 4 normal, 10 pancreatitis | Eligibility audit only: disease/chip perfectly confounded, including the 3-PDAC/10-pancreatitis FFPE subset |", "",
  "GSE59856's benign group is not a CP-only cohort. GSE85589's 29 non-cancer controls are 19 healthy plus 10 cholelithiasis, not 29 healthy. PC is the deposited pancreatic-cancer label; sample-level histologic subtype is not used to relabel those serum participants as proven PDAC. Distinct datasets are treated as separate studies; no documented participant overlap was identified among the analyzed cohorts, but explicit cross-study donor identifiers are unavailable.", "",
  "E-MTAB-1791 was initially deferred after metadata inspection. The later [publication audit](publication_audit_methods_results.md) records its current preparation/analysis status, clinical-history eligibility, donor limitations and the composition sensitivities for the existing tissue/culture cohorts. GSE179248 and GSE295071 share donors and must not be counted as two independent replications. Known GSE28735/GSE62452 overlap also requires care in later extensions.", "",
  "## Methods", "",
  "### Frozen Families and Annotation", "",
  "The dated specification is [external_validation_protocol.md](external_validation_protocol.md). We froze the existing 506 GSE208536 omnibus-associated genes and their spatial contrast directions before testing external outcomes. The current protocol is a transparent dated specification, not a registered prospective protocol. The named STAT1 example was already examined in the original reanalysis.", "",
  "Affymetrix probes use hgu133plus2.db. GSE91035's annotated supplement contains transcript IDs, zero labels, and spreadsheet-style date conversions; treating all of these as gene symbols would inflate the gene universe. The delivered analysis instead maps GPL22763 Entrez IDs to current symbols, then resolves remaining exact symbols/unambiguous aliases through org.Hs.eg.db. No date-to-gene guesses are made. Probe selection uses the highest overall mean across the study, not the disease effect. GSE91035 STAT1 is unresolved/unavailable after this annotation procedure, so no STAT1 effect is reported for that cohort. This annotation repair followed a superseded exploratory run and is disclosed in the protocol.", "",
  "Serum mature MIMAT accessions are mapped using the compact GPL18941 miRBase-v20 annotation. GSE85589's _st probe suffix is removed before the accession match. Ambiguous/multiple accessions and features with any missing values are excluded; duplicate names use the highest overall mean. Thus all-feature BH is across the eligible annotated/complete miRNAs, not every array feature. Coverage denominators and probe-level exclusions are deposited. Complete mapping does not establish signal above the assay detection limit; constant/background-level markers can remain uninformative with AUC 0.5.", "",
  "### Tissue Models", "",
  "Deposited log2 microarrays use robust, trend-moderated limma. GSE15471 blocks explicitly on donor after technical-replicate averaging. GSE91035 uses separate PDAC-versus-CP/normal/benign two-group models. Age/sex and benign/tumor pairing are unavailable; commercially sourced normal RNA differs in procurement. None of these GSE91035 contrasts is a definitive clinical-specificity test.", "",
  "GSE179248 contains integer featureCounts read counts. Unambiguous ENSEMBL-to-symbol mappings are summed by symbol, filtered with edgeR filterByExpr by day group, TMM normalized and modeled with voom/limma, donor + day. The 14 paired donors contribute 28 specimens, not 28 independent patients. Culture day, matrix exposure, isolation and other culture changes are inseparable from the observed ADM-associated response; this is not longitudinal human tumor development.", "",
  "All eligible genes receive BH correction separately within each study/contrast. A second, explicitly labeled BH family covers the measurable subset of the frozen 506 genes in each external contrast. Same-direction counts and Spearman correlations are descriptive; genes are not independent biological replicates. No pooled cross-platform effect, pooled CP conclusion, or binomial significance test of gene concordance is claimed. Pointwise 95% CIs are not simultaneous confidence bands.", "",
  "### Serum Models", "",
  "The fixed family is the 27 unique miRNAs in the historical six, published three, and existing exploratory EV candidate20 pool. The candidate pool was not FDR-significant in EV discovery and is not relabeled as a validated panel. Age/sex-adjusted limma compares PC against healthy, mixed benign, and pooled other cancers within each study. All eligible annotated features receive contrast-level BH; the frozen family also receives BH within each contrast and across all three contrasts per dataset.", "",
  "GSE59856 has unknown age for one PC case and all 52 liver-cancer controls: adjusted contrasts use 99 PC and 150 healthy/21 benign/248 other-cancer controls. Unadjusted single-marker ROC curves retain all otherwise eligible samples. ROC score direction is locked to the previously analyzed diagnostic-plasma mean difference, never selected to maximize serum AUC. Different analytes and scales prevent direct plasma-coefficient transfer. AUC CIs are pointwise DeLong intervals, with multiple subgroup/marker comparisons retained as exploratory rather than confirmatory wins.", "",
  "The GSE59856 subgroup uses 17 deposited pStage IIA/IIB cases and excludes ypStage codes; this does not guarantee a complete untreated history or establish screening performance. GSE85589's CA19-9-low subgroup uses 23 PC cases and restricts controls to the same <37 category (19 healthy and 7 cholelithiasis). Stage is unavailable in GSE85589. Its pooled other-cancer comparison includes 81 ICC specimens explicitly labeled validation as well as 34 initial other-cancer specimens. Diagnosis mixture, sampling context and unknown assay batches can affect that comparison; it is exploratory and is not a controlled pan-cancer specificity validation. Those other-cancer specimens are not used in the CA19-9 nested-CV models.", "",
  "GSE85589 deposits CA19-9 only as <37/>37/NA. Two benign participants lack CA19-9 and are excluded from incremental modeling, leaving 88 PC and 27 non-cancer controls. Three ridge logistic models (glmnet alpha=0) use baseline age, sex and binary CA19-9, then add the fixed published three or historical six. Shared stratified 10 outer folds estimate predictions; five inner folds choose lambda by binomial deviance from a fixed 60-value grid. Predictor standardization is fit by glmnet within training folds. No marker selection, serum-ROC sign tuning, or PLCO retuning is performed.", "",
  "The single deterministic outer split yields conditional out-of-fold-score CIs and DeLong comparisons; they omit pipeline/split uncertainty. The deposited arrays were normalized by their authors before this analysis; raw-array normalization was not repeated within CV folds. These CV estimates are therefore conditional on the supplied processed assay and are not end-to-end validation of a deployable preprocessing pipeline. Brier scores describe a case-control sample, not screening-risk calibration. These are newly fitted serum models with internal validation, not external validation of the original plasma coefficients, and cannot establish superiority over optimally modeled continuous CA19-9.", "",
  "## Tissue Results", "",
  "First, genome-wide/eligible-gene counts (not restricted to the spatial family):", "",
  md_table(model_table), "",
  "The next table counts external all-gene BH significance only among the frozen measurable spatial genes. The dual column additionally requires the matched spatial contrast to pass its original all-target BH threshold and the external effect to have the same sign. CP and benign comparisons are listed for transparency but do not reproduce a Normal-to-PDAC comparison directly.", "",
  md_table(summary_table), "",
  sprintf("Of the 53 measurable spatial Normal-to-ADM-significant genes in the culture dataset, 43 have the same direction and 35 additionally pass external all-gene BH <0.05. In the paired tumor cohort, 272 frozen genes are significant in both matched contrasts with concordant direction. The two-CP exploratory comparison has %d all-gene BH-significant associations, but none are in the frozen spatial family; zero frozen-family BH-significant genes are observed. Its full results are retained. Two controls and missing clinical adjustment prevent a robust specificity claim; failure of the frozen family to reach significance is not evidence of equivalence.", sum(cp_all$BH_all_genes < .05)), "",
  "### STAT1", "", md_table(stat1_table), "",
  "The independent culture increase is compatible with the spatial ADM-versus-Normal increase. The paired tumor increase supports a tumor/adjacent-normal association. Neither establishes the historical claim of a significant transient ADM peak: the original spatial PDAC-versus-ADM decrease remains non-significant, and the original age/sex-adjusted PDAC-versus-CP contrast remains non-significant.", "",
  "## Serum Results", "",
  sprintf("GSE59856 has 20/27 complete, unambiguously mapped frozen markers; the excluded markers are %s. GSE85589 has 27/27. Missing or ambiguous measurement is not biological non-replication.", paste(missing, collapse = ", ")), "",
  "The following significance count uses the frozen-family BH correction across all three comparisons per dataset. Direction is compared with diagnostic plasma. Adjusted association and unadjusted single-marker discrimination are separate analyses.", "",
  md_table(serum_summary), "",
  "### Binary CA19-9 Increment", "", md_table(cv_table), "",
  sprintf("Neither panel improves AUC convincingly: paired conditional DeLong BH q-values are %.3f and %.3f for the published-three and historical-six additions. The small control sample and lack of continuous CA19-9 limit this analysis.",
          comparison$BH_two_comparisons[1], comparison$BH_two_comparisons[2]), "",
  "## Figures", "",
  "Each figure has PNG and PDF versions in reanalysis/figures/external/.", "",
  "1. tissue_direction_concordance: all frozen genes measurable in paired tumor, commercial-normal, and paired culture comparisons; signs are descriptive.",
  "2. STAT1_external_forest: pointwise CIs across distinct spatial, paired-tissue, culture, and adjusted-CP contrasts; no pooled estimate.",
  "3. serum_fixed_marker_AUC: all nine historical/published markers, fixed plasma directions, healthy/benign/other-cancer comparisons; missing measurements remain blank.",
  "4. serum_CA19_9_nested_CV: shared-fold out-of-fold ROC curves for baseline and two fixed-panel additions.", "",
  "## Publication Interpretation", "",
  "The added data strengthen a focused cross-cohort tissue-state reproducibility analysis and provide a negative clinical-increment benchmark in serum. They do not rescue the original early-detection or malignant-switch claims. A journal-facing paper should distinguish paired culture response, cross-sectional tissue association, diagnostic case-control discrimination, and prediagnostic transportability throughout.", "",
  "The gaps identified after these analyses motivated the [publication follow-up](publication_audit_methods_results.md): clinical-profile corroboration, composition sensitivity and a larger inflammatory comparison. Direct spatial patient IDs, individual clinical fields for the larger cohort and sample-level PLCO lead time remain unavailable. The follow-up must be read before making cellular-origin or CP-specificity claims. Neither dataset accumulation nor attractive ROC curves guarantees novelty or journal acceptance; the full manuscript still requires author review and declarations.", "",
  "## Reproduction and Sources", "",
  "Run scripts/reanalysis/run_external.R after the original reanalysis. Dependencies are pinned in renv.lock, source hashes are saved, and scripts/reanalysis/17_validate_external.R checks sample pairing, BH recomputation, mapping, AUC/Brier recomputation, subgroup sizes and source hashes. Four figure pairs are additionally rendered for visual QA.", "",
  "- [GSE15471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471) and [published-QC curation](https://master.bioconductor.org/packages/devel/data/experiment/manuals/KEGGdzPathwaysGEO/man/KEGGdzPathwaysGEO.pdf).",
  "- [GSE179248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179248); [human ADM culture source article](https://doi.org/10.1016/j.gastha.2023.02.003).",
  "- [GSE91035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91035) and [GPL22763](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL22763).",
  "- [GSE59856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59856); [serum source article](https://doi.org/10.1371/journal.pone.0118220).",
  "- [GSE85589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85589) and [GPL18941 mature-miRNA annotation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18941).",
  "- [GSE101462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101462), excluded from inferential validation for perfect chip/disease confounding.",
  "- [glmnet methods and implementation](https://glmnet.stanford.edu/articles/glmnet.html)."
)
writeLines(report, "docs/external_validation_methods_results.md")
cat("External methods/results report generated.\n")
