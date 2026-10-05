#!/usr/bin/env Rscript

source("scripts/reanalysis/external_helpers.R")
out <- "reanalysis/results/publication"
fig <- "reanalysis/figures/publication"
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
export <- function(p, name, width, height) {
  ggplot2::ggsave(file.path(fig, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggplot2::ggsave(file.path(fig, paste0(name, ".pdf")), p, width = width, height = height)
}
models <- read.csv(file.path(out, "STAT1_composition_sensitivity.csv"))
labels <- c(Unadjusted_for_composition = "Base model", Plus_T_cells = "+ T-cell score",
            Plus_Fibroblasts = "+ Fibroblast score", Plus_both = "+ Both scores")
models$Model_Label <- factor(unname(labels[models$Variant]), levels = rev(unname(labels)))
models$Cohort <- paste(models$Dataset, models$Comparison)
facet_labels <- c("GSE208536 ADM_vs_Normal" = "Spatial: ADM - Normal\n22 profile/state means, 8 inferred cases",
                  "GSE208536 PDAC_vs_ADM" = "Spatial: PDAC - ADM\n22 profile/state means, 8 inferred cases",
                  "GSE208536 PDAC_vs_Normal" = "Spatial: PDAC - Normal\n22 profile/state means, 8 inferred cases",
                  "GSE15471 PDAC_vs_Normal" = "Paired tissue: PDAC - Normal\n35 donor pairs",
                  "GSE143754 PDAC_vs_CP" = "Bulk tissue: PDAC - CP\n11 PDAC / 6 CP; age and sex retained",
                  "GSE179248 Culture_D6_vs_D0" = "Human ADM culture: day 6 - day 0\n14 donor pairs; TMM/voom weights retained")
models$Cohort <- factor(models$Cohort, levels = names(facet_labels), labels = facet_labels)
estimated <- models[models$Status == "Estimated", ]
missing <- models[models$Status != "Estimated", ]
p <- ggplot2::ggplot(estimated, ggplot2::aes(log2FC, Model_Label)) +
  ggplot2::geom_vline(xintercept = 0, color = "#868B91", linetype = 2, linewidth = .35) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .15) +
  ggplot2::geom_point(ggplot2::aes(color = Variant), size = 2.4) +
  ggplot2::geom_text(data = missing, ggplot2::aes(x = .25, label = "Insufficient marker coverage"),
                     color = "#777777", size = 3, hjust = 0) +
  ggplot2::facet_wrap(~Cohort, ncol = 2, drop = FALSE) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::scale_color_manual(values = c(Unadjusted_for_composition = "#333333", Plus_T_cells = "#168273",
                                        Plus_Fibroblasts = "#C04B58", Plus_both = "#457BB5")) +
  ggplot2::labs(title = "STAT1 sensitivity to immune and stromal expression proxies",
                subtitle = "Original blocking/clinical covariates retained; estimates are not pooled",
                x = "Log2 expression difference (pointwise 95% CI)", y = NULL,
                caption = "MCP-counter scores are expression proxies, not cell fractions.\nSpatial fibroblast coverage: 2 of 8 markers; adjustment not estimated.") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 9),
                 plot.caption = ggplot2::element_text(hjust = 0),
                 panel.spacing = grid::unit(1.0, "lines"))
export(p, "STAT1_composition_sensitivity", 11, 8)

roi <- read.csv(file.path(out, "GSE208536_published_case_crosswalk.csv"))
counts <- read.csv(file.path(out, "GSE208536_published_case_roi_counts.csv"))
missing_states <- unlist(lapply(c("Normal", "ADM", "PDAC"), function(state) {
  cases <- counts$Published_Case[counts[[state]] == 0]
  if (length(cases)) paste(paste(cases, collapse = ", "), "lacks", state) else character()
}))
coverage <- read.csv(file.path(out, "MCPcounter_marker_coverage.csv"))
fmt <- function(x, digits = 3) ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
pv <- function(x) ifelse(is.na(x), NA_character_, format.pval(x, digits = 3, eps = 1e-6))
md_table <- function(x) {
  x[] <- lapply(x, function(v) ifelse(is.na(v), "Not estimated", as.character(v)))
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
table <- data.frame(Dataset = models$Dataset, Comparison = models$Comparison,
                    Model = unname(labels[models$Variant]), Effect = fmt(models$log2FC),
                    CI = ifelse(is.na(models$log2FC), NA_character_,
                                paste(fmt(models$CI_low), "to", fmt(models$CI_high))),
                    P = pv(models$P), Status = models$Status, check.names = FALSE)
report <- c("# Publication Audit: Metadata, Composition And Inflammatory Controls", "",
            "Updated 2026-10-05. These are post hoc sensitivity analyses. The earlier primary results,",
            "frozen families and PLCO predictions are unchanged. This report does not certify journal readiness.", "",
            "## Spatial Case Identities", "",
            "All eight distinct age/sex/grade/stage profiles match unique cases in Table 1 of the",
            "[source article](https://doi.org/10.3389/fimmu.2022.961457). This corroborates the clinical-profile",
            "crosswalk; it is not a directly deposited ROI-to-patient ID. A mislabeled profile would remain",
            "undetectable from this comparison. The unusual ROI counts are retained, not repaired by guessing.", "",
            md_table(counts[, c("Published_Case", "Normal", "ADM", "PDAC")]), "",
            "The 48 ROIs yield 22 observed case/state means. Six inferred cases have all three states.",
            paste0(paste(missing_states, collapse = "; "), ". The source's general selection description does"),
            "not imply that every deposited case has exactly six ROIs. The original grouping is preserved.", "",
            "## STAT1 Composition Sensitivity", "",
            "[MCP-counter](https://doi.org/10.1186/s13059-016-1070-5) implementation and signatures are pinned",
            "to commit b6eac73e91c246fcff0bb1a5c68a816cd588fc48. A named score must have at least three markers",
            "and at least half its published markers to enter a model. This analyst-selected coverage rule",
            "does not validate a truncated signature for a targeted spatial assay. Scores are standardized",
            "within each cohort and are not cell fractions. STAT1 is not a signature marker. Baseline models",
            "retain case/donor blocking or age and sex; ordinary spatial models and full-gene moderated",
            "microarray models match the primary analysis. Covariate adjustment is conditional, not causal.", "",
            md_table(table), "",
            "In GSE15471, fibroblast adjustment reduces the STAT1 tumor/normal estimate from +0.940",
            "to +0.010 log2 units (95% CI -0.609 to +0.629). The broader conditional CI includes zero.",
            "This is consistent with strong dependence on stromal-associated expression but does not prove",
            "that fibroblasts explain the signal or identify the expressing compartment. Adjustment may",
            "remove mediated biology and is affected by collinearity.", "",
            "The spatial ADM/Normal estimate is +0.550 after T-cell adjustment (95% CI +0.285 to +0.816).",
            "Only two of eight fibroblast markers are measured, so that spatial sensitivity cannot be fit.",
            "The later PDAC/ADM contrast remains uncertain. There is no confirmed transient STAT1 peak.",
            "The 11-PDAC/6-CP cohort provides no significant STAT1 difference in these adjusted variants.", "",
            "All variants, unavailable models, pointwise CIs, exploratory model-family BH values, residual",
            "degrees of freedom and condition numbers are saved. Do not select a favorable adjusted p-value.", "",
            "![STAT1 sensitivity](../reanalysis/figures/publication/STAT1_composition_sensitivity.png)", "",
            "## Larger Inflammatory Cohort", "",
            "E-MTAB-1791 contains 457 deposited records. The eligibility rule retains 195 explicitly PDAC",
            "specimens and 59 pancreatitis specimens without an associated-tumor clinical-history label.",
            "Nine tumor-associated pancreatitis specimens are excluded. The",
            "[source publication](https://doi.org/10.1002/ijc.31087) reports 59 CP specimens from 58 patients.",
            "The repeated patient cannot be identified from the SDRF. It reports aggregate median ages of",
            "47.1 years in CP and 63.4 years in PDAC, but individual age/sex fields are not deposited in the SDRF.",
            "Group means/medians are not substituted as covariates. The source's 452 quality-controlled",
            "samples differ from the 457 deposited records; the primary tissue counts match the source.", "")
culture <- models[models$Dataset == "GSE179248" & models$Variant == "Plus_Fibroblasts", ]
report <- append(report, c(
  "### Additional Post Hoc Culture Check", "",
  "After inspecting the initial three-cohort sensitivities, the same named proxies were added",
  "to the already analyzed paired human ADM culture cohort. Donor blocking, outcome-independent",
  "count filtering and TMM normalization were retained; voom weights were recomputed for each design.",
  sprintf("Fibroblast-score adjustment changes the culture STAT1 effect from +0.585 to %s log2 units (95%% CI %s to %s; p=%s).",
          fmt(culture$log2FC), fmt(culture$CI_low), fmt(culture$CI_high), pv(culture$P)),
  "The adjusted culture interval includes zero. T-cell coverage fails the fixed eligibility rule,",
  "so T-cell and combined models are not fitted. Expression proxies can change with cell state",
  "as well as composition; this is not evidence of actual fibroblast accumulation in culture.", ""),
  after = match("## Larger Inflammatory Cohort", report) - 1)

cohort_file <- file.path(out, "E-MTAB-1791_replication_summary.csv")
if (file.exists(cohort_file)) {
  summary <- read.csv(cohort_file)
  family <- read.csv(file.path(out, "E-MTAB-1791_frozen_spatial_family.csv"))
  cp_stat1 <- read.csv(file.path(out, "E-MTAB-1791_STAT1_composition_sensitivity.csv"))
  stat1_family <- family[family$Gene == "STAT1", ]
  cp_stat1$BH_frozen_base_family <- NA_real_
  base_rows <- which(cp_stat1$Variant == "Unadjusted_for_composition")
  cp_stat1$BH_frozen_base_family[base_rows] <-
    stat1_family$BH_frozen_family[match(cp_stat1$Model[base_rows], stat1_family$Model)]
  cp_models <- read.csv(file.path(out, "E-MTAB-1791_models.csv"))
  omit <- read.csv(file.path(out, "E-MTAB-1791_omit_one_CP_STAT1.csv"))
  envelope <- read.csv(file.path(out, "E-MTAB-1791_omit_one_CP_envelope.csv"))
  chip_stat1 <- cp_stat1[cp_stat1$Model == "Chip_adjusted" &
                         cp_stat1$Variant == "Unadjusted_for_composition", ]
  chip_envelope <- envelope[envelope$Model == "Chip_adjusted" & envelope$Gene == "STAT1", ]
  base_designs <- cp_models[cp_models$Variant == "Unadjusted_for_composition", ]
  design_table <- base_designs[, c("Model", "Status", "Chips", "Mixed_Diagnosis_Chips",
                                   "Design_Columns", "Design_Rank", "Residual_DF")]
  omission_table <- do.call(rbind, lapply(split(omit, omit$Model), function(s) {
    estimated <- s[s$Status == "Estimated", ]
    e <- envelope[envelope$Model == s$Model[1] & envelope$Gene == "STAT1", ]
    data.frame(Model = s$Model[1], Estimable_Omissions = nrow(estimated),
               STAT1_Effect_Range = if (nrow(estimated))
                 paste(fmt(min(estimated$log2FC)), "to", fmt(max(estimated$log2FC))) else NA_character_,
               STAT1_AllGene_BH_Range = if (nrow(estimated))
                 paste(pv(min(estimated$BH_all_genes)), "to", pv(max(estimated$BH_all_genes))) else NA_character_,
               STAT1_Frozen_BH_Range = paste(pv(e$Min_BH_frozen), "to", pv(e$Max_BH_frozen)),
               STAT1_Frozen_Significant_All_Omissions = e$Significant_All_Estimable_Omissions)
  }))
  family$Model <- factor(family$Model, levels = c("Unadjusted", "Chip_adjusted"),
                         labels = c("Unadjusted", "Chip fixed effects"))
  family$Replication <- family$BH_frozen_family < .05 & family$Same_direction
  p <- ggplot2::ggplot(family, ggplot2::aes(Spatial_log2FC, log2FC, color = Replication)) +
    ggplot2::geom_hline(yintercept = 0, color = "#AAB1B8", linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, color = "#AAB1B8", linewidth = .3) +
    ggplot2::geom_point(alpha = .65, size = 1.5) +
    ggplot2::geom_point(data = family[family$Gene == "STAT1", ], shape = 21, fill = "white",
                        color = "#333333", size = 3, stroke = 1) +
    ggplot2::facet_wrap(~Model, nrow = 1) +
    ggplot2::scale_color_manual(values = c("TRUE" = "#168273", "FALSE" = "#9B9DA1"),
                                labels = c("FALSE" = "Other frozen genes", "TRUE" = "Family BH < .05, same direction")) +
    ggplot2::labs(title = "Frozen spatial family in a larger inflammatory comparison",
                  subtitle = "195 PDAC / 59 CP specimens; unverified repeat and unavailable individual age/sex",
                  x = "Spatial PDAC - Normal (log2 difference)", y = "E-MTAB-1791 PDAC - CP (log2 difference)",
                  color = NULL, caption = "Different tissue contrasts, not a pooled effect. White outlined point: STAT1.") +
    ggplot2::theme_classic(base_size = 10) + ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
  export(p, "inflammatory_cohort_frozen_family", 10, 5.4)
  report <- c(report, "The source normalized log2 column is parsed for each eligible array. Source byte counts,",
              "ZIP CRCs where used, MD5s, probe order and annotations are checked. Entrez IDs map to",
              "unambiguous current symbols; highest overall mean selects duplicate probes without using labels.",
              "Prepared probe/gene matrices and selected mappings are retained. Models use robust, trend-moderated",
              "limma. Unadjusted and estimable chip-fixed-effect contrasts are reported with all-gene BH and",
              "separate BH across the measurable subset of the frozen 506 spatial genes.", "",
              "### Design And Frozen-Family Results", "",
              md_table(design_table), "", md_table(summary), "",
              "Gene-direction counts and correlations are descriptive. The spatial PDAC/Normal and",
              "bulk PDAC/CP comparisons have different controls and do not test the same tissue contrast.", "",
              "### STAT1 In The Larger Cohort", "",
              md_table(data.frame(Model = cp_stat1$Model, Variant = cp_stat1$Variant,
                                   log2FC = fmt(cp_stat1$log2FC), CI_low = fmt(cp_stat1$CI_low),
                                   CI_high = fmt(cp_stat1$CI_high), P = pv(cp_stat1$P),
                                   BH_all_genes = pv(cp_stat1$BH_all_genes),
                                   BH_frozen_base_family = pv(cp_stat1$BH_frozen_base_family))), "",
              "The frozen-family BH column is evaluated for the two base models only. Composition",
              "variants retain all-gene BH; their smaller exploratory model-family BH is saved separately.",
              "These corrections answer different multiplicity questions and must not be interchanged.", "",
              sprintf("The chip-adjusted base STAT1 effect is +%s log2 units (95%% CI %s to %s), with frozen-family BH=%s but all-gene BH=%s.",
                      fmt(chip_stat1$log2FC), fmt(chip_stat1$CI_low), fmt(chip_stat1$CI_high),
                      pv(chip_stat1$BH_frozen_base_family), pv(chip_stat1$BH_all_genes)),
              "It passes the narrower frozen-family threshold, not the all-gene threshold. The composition",
              "variants must not be selected merely because their adjusted p-values are more favorable.", "",
              "Every possible CP specimen is omitted in turn for each estimable base model. The saved",
              "gene-level envelope reports effect and family-BH extrema, direction stability and the number",
              "of estimable omissions. This is a sensitivity analysis for the unidentified repeat, not donor",
              "blocking or verified independence. It does not remove the substantial age difference.", "",
              md_table(omission_table), "",
              sprintf("Across all 59 chip-adjusted CP omissions, STAT1 remains positive (%s to %s), but frozen-family BH ranges from %s to %s and crosses 0.05.",
                      fmt(chip_envelope$Min_log2FC), fmt(chip_envelope$Max_log2FC),
                      pv(chip_envelope$Min_BH_frozen), pv(chip_envelope$Max_BH_frozen)),
              "Thus its direction is stable in this check, while threshold-level significance is not.", "",
              "![Frozen family](../reanalysis/figures/publication/inflammatory_cohort_frozen_family.png)", "",
              "A PDAC/CP association in this cohort is not proof of malignancy specificity, because individual",
              "clinical adjustment and a verified donor map are unavailable. Bulk composition scores also",
              "cannot establish cellular origin or eliminate epithelial/acinar lineage effects.", "")
} else {
  report <- c(report, "**Expression preparation and analysis are incomplete.** The eligibility audit is complete,",
              "but this cohort must not be counted as analyzed validation evidence. No partial cohort effects",
              "are tested or reported. The preparation script is resumable and checks all 254 specimens.", "")
}
report <- c(report, "## Consequences For Submission", "",
            "The strongest current framing is a metadata-aware reproducibility and transportability study,",
            "with separated tissue and circulating-marker evidence streams. The focused",
            "[prior-art assessment](publication_prior_art.md) documents overlap with the source analyses and",
            "existing tissue meta-analyses. It does not prove novelty.", "",
            "Submission still requires a coherent manuscript that replaces the historical claims, attribution",
            "to every source study, reconciled supplements, author review of the computational results and",
            "journal-specific reporting/declarations. Direct ROI/patient IDs, individual clinical adjustment",
            "for the larger cohort and untouched prediagnostic validation data remain unavailable. These",
            "limits must be disclosed and the claims narrowed; they cannot be solved by more significant tests.", "",
            "## Reproduction", "",
            "See [the dated sensitivity specification](publication_sensitivity_protocol.md), scripts 18-23,",
            "the source/array manifests in data/publication, and reanalysis/results/publication. All new figures",
            "are supplied in PNG and PDF. The old PLCO scores and thresholds are not retuned.")
report <- c(report, "", "## Output Index", "",
            "| Artifact | Contents |", "| --- | --- |",
            "| [Spatial crosswalk](../reanalysis/results/publication/GSE208536_published_case_crosswalk.csv) | All 48 ROIs, inferred profiles and corroborating published case labels |",
            "| [Composition models](../reanalysis/results/publication/STAT1_composition_sensitivity.csv) | All 24 named variants, including unavailable adjustments |",
            "| [Marker coverage](../reanalysis/results/publication/MCPcounter_marker_coverage.csv) | All ten populations in four existing tissue/culture cohorts |",
            "| [Source inputs](../data/publication/publication_download_manifest.csv) | Source URLs, byte counts and checksums |",
            "| [Eligibility audit](../reanalysis/results/publication/E-MTAB-1791_eligibility_mapping.csv) | All 457 records and exact inclusion/exclusion roles |")
if (file.exists(cohort_file)) report <- c(report,
            "| [Array manifest](../data/publication/E-MTAB-1791_array_manifest.csv) | All 254 source arrays and checked preparation records |",
            "| [Gene/probe map](../data/publication/E-MTAB-1791_selected_probe_map.csv) | Outcome-blind selected probes and symbol mappings |",
            "| [Larger-cohort models](../reanalysis/results/publication/E-MTAB-1791_models.csv) | All eight designs, estimability, residual degrees of freedom and condition numbers |",
            "| [All-gene contrasts](../reanalysis/results/publication/E-MTAB-1791_all_gene_contrasts.csv) | Unadjusted and chip-adjusted effects, pointwise CIs and all-gene BH |",
            "| [Frozen family](../reanalysis/results/publication/E-MTAB-1791_frozen_spatial_family.csv) | Measurable frozen genes, directions and separate family BH |",
            "| [Omission envelope](../reanalysis/results/publication/E-MTAB-1791_omit_one_CP_envelope.csv) | All possible single-CP omissions and gene-level sensitivity extrema |",
            "| [Output checks](../reanalysis/provenance/publication_output_checks.txt) | Numerical/structural checks; not certification of biological truth |")
writeLines(report, "docs/publication_audit_methods_results.md")
cat("Publication audit report and available figure pairs saved.\n")
