#!/usr/bin/env Rscript

source("scripts/reanalysis/external_helpers.R")
out <- "reanalysis/results/publication"
read_result <- function(name) read.csv(file.path(out, name))
roi <- read_result("GSE208536_published_case_crosswalk.csv")
counts <- read_result("GSE208536_published_case_roi_counts.csv")
stopifnot(nrow(roi) == 48, length(unique(roi$Published_Case)) == 8,
          all(vapply(split(roi$Original_Profile, roi$Published_Case), function(x) length(unique(x)) == 1, logical(1))),
          all(vapply(split(roi$Published_Case, roi$Original_Profile), function(x) length(unique(x)) == 1, logical(1))),
          sum(counts$Normal) == 16, sum(counts$ADM) == 16, sum(counts$PDAC) == 16,
          sum(counts$Normal > 0 & counts$ADM > 0 & counts$PDAC > 0) == 6)
report <- paste(readLines("docs/publication_audit_methods_results.md"), collapse = "\n")
for (state in c("Normal", "ADM", "PDAC")) {
  missing <- counts$Published_Case[counts[[state]] == 0]
  if (length(missing))
    stopifnot(grepl(paste(paste(missing, collapse = ", "), "lacks", state), report, fixed = TRUE))
}
composition <- read_result("STAT1_composition_sensitivity.csv")
stopifnot(nrow(composition) == 24,
          sum(composition$Status[composition$Dataset != "GSE179248"] == "Estimated") == 14,
          all(composition$CI_low[!is.na(composition$log2FC)] <= composition$log2FC[!is.na(composition$log2FC)]),
          all(composition$CI_high[!is.na(composition$log2FC)] >= composition$log2FC[!is.na(composition$log2FC)]),
          all(composition$Residual_DF[composition$Status == "Estimated"] > 0),
          isTRUE(all.equal(composition$BH_exploratory_model_family, p.adjust(composition$P, "BH"))))
spatial <- read.csv("reanalysis/results/spatial/GSE208536_patient_blocked_all_targets.csv")
spatial <- spatial[spatial$Gene == "STAT1", ]
bulk <- read.csv("reanalysis/results/bulk/GSE143754_gene_level_contrasts.csv")
bulk <- bulk[bulk$Gene == "STAT1", ]
tissue <- read.csv(file.path(external_dir, "tissue_all_gene_contrasts.csv"))
original <- c(GSE208536_ADM_vs_Normal = spatial$Normal_to_ADM_log2,
              GSE208536_PDAC_vs_ADM = spatial$ADM_to_PDAC_log2,
              GSE208536_PDAC_vs_Normal = spatial$Normal_to_PDAC_log2,
              GSE15471_PDAC_vs_Normal = tissue$log2FC[tissue$Dataset == "GSE15471" & tissue$Gene == "STAT1"],
              GSE143754_PDAC_vs_CP = bulk$PDAC_vs_CP_log2FC,
              GSE179248_Culture_D6_vs_D0 = tissue$log2FC[tissue$Dataset == "GSE179248" & tissue$Gene == "STAT1"])
base <- composition[composition$Variant == "Unadjusted_for_composition", ]
stopifnot(max(abs(base$log2FC - original[paste(base$Dataset, base$Comparison, sep = "_")])) < 1e-10)
coverage <- read_result("MCPcounter_marker_coverage.csv")
stopifnot(nrow(coverage) == 40,
          all(coverage$Present_Markers <= coverage$Published_Markers),
          all(coverage$Eligible_Adjustment == (coverage$Present_Markers >= 3 & coverage$Coverage >= .5)))
marker <- read.delim("data/publication/MCPcounter_genes.txt", check.names = FALSE)
stopifnot(!"STAT1" %in% marker[["HUGO symbols"]])
manifest <- read.csv("data/publication/publication_download_manifest.csv")
stopifnot(identical(unname(tools::md5sum(file.path("data/publication", manifest$File))), manifest$MD5))
if ("--composition-only" %in% commandArgs(trailingOnly = TRUE)) {
  writeLines(c("Publication composition checks: PASS", "48 ROIs and eight published clinical-profile crosswalks",
                "Original model variants reproduce earlier estimates; insufficient-coverage models retained",
                "Baseline STAT1 effects reproduce previous results within 1e-10",
                "Coverage, CIs, degrees of freedom, exploratory BH and source hashes checked"),
              "reanalysis/provenance/publication_composition_checks.txt")
  cat("Publication composition checks: PASS\n")
  quit(save = "no", status = 0)
}

meta <- read_result("E-MTAB-1791_eligibility_mapping.csv")
stopifnot(nrow(meta) == 457, sum(meta$Role == "PDAC") == 195, sum(meta$Role == "CP") == 59,
          sum(meta$Role == "Pancreatitis_associated_context") == 9,
          all(meta$Primary_Eligible == meta$Role %in% c("CP", "PDAC")))
arrays <- read.csv("data/publication/E-MTAB-1791_array_manifest.csv")
parser <- read.csv("data/publication/E-MTAB-1791_parser_equivalence.csv")
stopifnot(nrow(parser) == 1, parser$Exact_Base_Parser_Equivalence,
          parser$Source_MD5 == arrays$Source_MD5[match(parser$File, arrays$File)],
          parser$Probe_Count == 49577)
stopifnot(nrow(arrays) == 254, !anyDuplicated(arrays$Specimen), !anyNA(arrays$Source_MD5),
          all(nchar(arrays$Source_MD5) == 32), all(arrays$Probe_Count == 49577),
          setequal(arrays$Specimen, meta$Specimen[meta$Primary_Eligible]))
gene <- read.csv(gzfile("data/publication/E-MTAB-1791_primary_normalized_genes.csv.gz"),
                 check.names = FALSE, colClasses = c("character", rep("numeric", 254)))
map <- read.csv("data/publication/E-MTAB-1791_selected_probe_map.csv")
stopifnot(ncol(gene) == 255, !anyDuplicated(gene$Gene), all(is.finite(as.matrix(gene[, -1]))),
          identical(gene$Gene, map$Gene), setequal(names(gene)[-1], arrays$Specimen), "STAT1" %in% gene$Gene)
probes <- read.csv(gzfile("data/publication/E-MTAB-1791_primary_normalized_probes.csv.gz"),
                   check.names = FALSE, colClasses = c("character", rep("numeric", 254)))
probe_index <- match(map$ProbeID, probes$ProbeID)
stopifnot(nrow(probes) == 49577, !anyNA(probe_index), !anyDuplicated(probes$ProbeID),
          identical(names(probes)[-1], names(gene)[-1]),
          isTRUE(all.equal(unname(as.matrix(gene[, -1])),
                           unname(as.matrix(probes[probe_index, -1])), tolerance = 1e-12)))
rm(probes)
invisible(gc())
stats <- read_result("E-MTAB-1791_all_gene_contrasts.csv")
family <- read_result("E-MTAB-1791_frozen_spatial_family.csv")
reference <- read.csv(file.path(external_dir, "frozen_spatial_reference.csv"))
for (s in split(stats, stats$Model)) {
  stopifnot(nrow(s) == nrow(gene), !anyDuplicated(s$Gene), all(s$CI_low <= s$log2FC),
            all(s$CI_high >= s$log2FC), isTRUE(all.equal(s$BH_all_genes, p.adjust(s$P, "BH"))))
}
for (s in split(family, family$Model)) {
  stopifnot(all(s$Gene %in% reference$Gene),
            isTRUE(all.equal(s$BH_frozen_family, p.adjust(s$P, "BH"))),
            identical(s$Same_direction, sign(s$log2FC) == sign(s$Spatial_log2FC)))
  full <- stats[stats$Model == s$Model[1], ]
  index <- match(s$Gene, full$Gene)
  stopifnot(!anyNA(index), isTRUE(all.equal(s$log2FC, full$log2FC[index])),
            isTRUE(all.equal(s$P, full$P[index])))
}
models <- read_result("E-MTAB-1791_models.csv")
stopifnot(nrow(models) == 8, all(models$PDAC_Specimens == 195), all(models$CP_Specimens == 59),
          all(!models$Deposited_Patient_ID_Available), all(!models$Age_Sex_Adjustment_Available),
          all(models$Design_Columns[models$Status == "Estimated"] == models$Design_Rank[models$Status == "Estimated"]))
omissions <- read_result("E-MTAB-1791_omit_one_CP_STAT1.csv")
stopifnot(all(table(omissions$Model) == 59), all(omissions$Omitted_CP %in% meta$Specimen[meta$Role == "CP"]))
analysis_meta <- meta[match(names(gene)[-1], meta$Specimen), ]
analysis_meta$State <- factor(analysis_meta$Role, levels = c("CP", "PDAC"))
analysis_meta$Chip <- factor(analysis_meta$Chip)
stat1_expression <- as.numeric(unlist(gene[gene$Gene == "STAT1", -1], use.names = FALSE))
composition_scores <- read_result("E-MTAB-1791_composition_scores.csv")
composition_scores <- composition_scores[match(analysis_meta$Specimen, composition_scores$Specimen), ]
stopifnot(identical(composition_scores$Specimen, analysis_meta$Specimen))
analysis_meta$T_Score <- composition_scores$T_Score
analysis_meta$Fibroblast_Score <- composition_scores$Fibroblast_Score
stat1_models <- read_result("E-MTAB-1791_STAT1_composition_sensitivity.csv")
stopifnot(nrow(stat1_models) == sum(models$Status == "Estimated"),
          all(stat1_models$CI_low <= stat1_models$log2FC),
          all(stat1_models$CI_high >= stat1_models$log2FC),
          isTRUE(all.equal(stat1_models$BH_exploratory_model_family, p.adjust(stat1_models$P, "BH"))))
base_terms <- list(Unadjusted = "State", Chip_adjusted = c("Chip", "State"))
variant_terms <- list(Unadjusted_for_composition = character(), Plus_T_cells = "T_Score",
                      Plus_Fibroblasts = "Fibroblast_Score", Plus_both = c("T_Score", "Fibroblast_Score"))
for (i in seq_len(nrow(stat1_models))) {
  row <- stat1_models[i, ]
  terms <- c(base_terms[[row$Model]], variant_terms[[row$Variant]])
  design <- model.matrix(reformulate(terms), analysis_meta)
  estimate <- lm.fit(design, stat1_expression)$coefficients["StatePDAC"]
  stopifnot(abs(row$log2FC - estimate) < 1e-10)
}
for (model in names(base_terms)) {
  saved <- stats[stats$Model == model & stats$Gene == "STAT1", ]
  if (!nrow(saved)) next
  design <- model.matrix(reformulate(base_terms[[model]]), analysis_meta)
  estimate <- lm.fit(design, stat1_expression)$coefficients["StatePDAC"]
  stopifnot(abs(saved$log2FC - estimate) < 1e-10)
  for (j in which(analysis_meta$State == "CP")) {
    row <- omissions[omissions$Model == model & omissions$Omitted_CP == analysis_meta$Specimen[j], ]
    stopifnot(nrow(row) == 1)
    if (row$Status != "Estimated") next
    design <- model.matrix(reformulate(base_terms[[model]]), droplevels(analysis_meta[-j, ]))
    estimate <- lm.fit(design, stat1_expression[-j])$coefficients["StatePDAC"]
    stopifnot(abs(row$log2FC - estimate) < 1e-10)
  }
}
omitted_family <- read_result("E-MTAB-1791_omit_one_CP_frozen_family.csv")
for (s in split(omitted_family, paste(omitted_family$Model, omitted_family$Omitted_CP))) {
  stopifnot(!anyDuplicated(s$Gene), isTRUE(all.equal(s$BH_frozen_family, p.adjust(s$P, "BH"))))
}
envelope <- read_result("E-MTAB-1791_omit_one_CP_envelope.csv")
stopifnot(all(envelope$Min_log2FC <= envelope$Max_log2FC), all(envelope$Min_BH_frozen <= envelope$Max_BH_frozen),
          all(envelope$Estimable_Omissions <= 59))
for (i in seq_len(nrow(envelope))) {
  row <- envelope[i, ]
  s <- omitted_family[omitted_family$Model == row$Model & omitted_family$Gene == row$Gene, ]
  original <- family[family$Model == row$Model & family$Gene == row$Gene, ]
  stopifnot(nrow(s) == row$Estimable_Omissions,
            abs(min(s$log2FC) - row$Min_log2FC) < 1e-12,
            abs(max(s$log2FC) - row$Max_log2FC) < 1e-12,
            abs(min(s$BH_frozen_family) - row$Min_BH_frozen) < 1e-12,
            abs(max(s$BH_frozen_family) - row$Max_BH_frozen) < 1e-12,
            row$Significant_All_Estimable_Omissions == all(s$BH_frozen_family < .05),
            row$Direction_Stable == all(sign(s$log2FC) == sign(original$log2FC)))
}
summary <- read_result("E-MTAB-1791_replication_summary.csv")
stopifnot(nrow(summary) == length(unique(family$Model)))
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  s <- family[family$Model == row$Model, ]
  a <- stats[stats$Model == row$Model, ]
  e <- envelope[envelope$Model == row$Model, ]
  stopifnot(row$All_Genes_Tested == nrow(a),
            row$All_Gene_BH_Significant == sum(a$BH_all_genes < .05),
            row$Measurable_Frozen == nrow(s),
            row$Frozen_BH_Significant == sum(s$BH_frozen_family < .05),
            row$Same_Direction_As_Spatial == sum(s$Same_direction),
            row$Spatial_And_CP_Concordant_Significant ==
              sum(s$Spatial_BH < .05 & s$BH_all_genes < .05 & s$Same_direction),
            row$Frozen_Significant_All_Omissions == sum(e$Significant_All_Estimable_Omissions),
            abs(row$Spearman_Effect_Correlation -
                  cor(s$log2FC, s$Spatial_log2FC, method = "spearman")) < 1e-12)
}
for (path in list.files("scripts/reanalysis", pattern = "\\.R$", full.names = TRUE)) invisible(parse(path))
figures <- list.files("reanalysis/figures/publication", pattern = "\\.pdf$", full.names = TRUE)
stopifnot(length(figures) == 2, all(file.exists(sub("\\.pdf$", ".png", figures))))
writeLines(c("Publication output checks: PASS", "Eight unique published clinical profiles crosswalked, not direct IDs",
              "48 ROIs / 22 observed case-state means; no relabeling of unexpected counts",
              "Original STAT1 base effects reproduced; failed composition models retained",
              "195 PDAC and 59 primary CP specimens; 9 tumor-associated CP specimens excluded",
              "254 source arrays checked; frozen-family BH, CIs and omission envelopes checked",
              "All eight STAT1 model effects and all base-model omissions independently recomputed by ordinary least squares",
              "Frozen-family summary counts and effect correlations independently recomputed",
              "No claim of 254 independent patients or individual clinical adjustment",
              "All reanalysis R scripts parse; two PNG/PDF figure pairs exist"),
            "reanalysis/provenance/publication_output_checks.txt")
cat("Publication analysis checks: PASS\n")
