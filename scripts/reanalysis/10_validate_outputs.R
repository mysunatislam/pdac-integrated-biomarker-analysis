#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
root <- "reanalysis"
spatial <- read.csv(file.path(root, "results/spatial/GSE208536_patient_blocked_all_targets.csv"))
bulk <- read.csv(file.path(root, "results/bulk/GSE143754_gene_level_contrasts.csv"))
integrated <- read.csv(file.path(root, "results/integrated_tissue/CP_aware_spatial_associated_genes.csv"))
spatial_overlap <- read.csv(file.path(root, "results/integrated_tissue/jointly_measurable_spatial_associated_genes.csv"))
plasma <- read.csv(file.path(root, "results/plasma/diagnostic_to_PLCO_performance.csv"))
scores <- read.csv(file.path(root, "results/plasma/PLCO_locked_model_scores.csv"))
diagnostic_scores <- read.csv(file.path(root, "results/plasma/diagnostic_model_scores.csv"))
ev <- read.csv(file.path(root, "results/ev/GSE304572_EV_differential.csv"))
stopifnot(nrow(spatial) == 1825, !anyDuplicated(spatial$Gene), !anyDuplicated(bulk$Gene),
          nrow(plasma) == 4, nrow(scores) == 96, sum(scores$Outcome) == 48,
          nrow(diagnostic_scores) == 203, sum(diagnostic_scores$Outcome) == 121,
          !anyDuplicated(scores$Sample), !anyDuplicated(diagnostic_scores$Sample),
          !any(scores$Sample %in% diagnostic_scores$Sample),
          all(integrated$PDAC_vs_CP_BH_FDR < 0.05),
          all(abs(integrated$PDAC_vs_CP_log2FC) >= 0.5),
          all(integrated$Patient_Blocked_F_BH_FDR_all_targets < 0.05))
significant_steps <- with(spatial, Normal_to_ADM_BH_FDR_all_targets < 0.05 &
                                  ADM_to_PDAC_BH_FDR_all_targets < 0.05)
stopifnot(identical(spatial$ADM_peak,
                   significant_steps & spatial$Normal_to_ADM_log2 > 0 & spatial$ADM_to_PDAC_log2 < 0),
          identical(spatial$ADM_trough,
                   significant_steps & spatial$Normal_to_ADM_log2 < 0 & spatial$ADM_to_PDAC_log2 > 0),
          !any(spatial$ADM_peak & spatial$ADM_trough))
stopifnot(identical(spatial_overlap$Bulk_Spatial_NormalPDAC_Direction_Concordant,
                   sign(spatial_overlap$PDAC_vs_Normal_log2FC) == sign(spatial_overlap$Normal_to_PDAC_log2)),
          identical(spatial_overlap$BulkCP_SpatialNormalPDAC_Direction_Concordant,
                   sign(spatial_overlap$PDAC_vs_CP_log2FC) == sign(spatial_overlap$Normal_to_PDAC_log2)))
for (table in list(spatial, bulk, ev)) {
  fdr_columns <- grep("FDR", names(table), value = TRUE)
  for (column in fdr_columns) {
    stopifnot(all(is.finite(table[[column]])), all(table[[column]] >= 0 & table[[column]] <= 1))
  }
}
for (column in c("PLCO_AUC", "PLCO_Sensitivity", "PLCO_Specificity", "PLCO_Brier")) {
  stopifnot(all(is.finite(plasma[[column]])), all(plasma[[column]] >= 0 & plasma[[column]] <= 1))
}
for (table in list(scores, diagnostic_scores)) {
  prediction_columns <- setdiff(names(table), c("Sample", "Outcome"))
  for (column in prediction_columns) {
    stopifnot(all(is.finite(table[[column]])), all(table[[column]] >= 0 & table[[column]] <= 1))
  }
}
stopifnot(all(plasma$PLCO_AUC_CI_low <= plasma$PLCO_AUC),
          all(plasma$PLCO_AUC_CI_high >= plasma$PLCO_AUC),
          all(plasma$Model %in% names(scores)))
pngs <- list.files(file.path(root, "figures"), pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
pdfs <- sub("\\.png$", ".pdf", pngs)
all_pdfs <- list.files(file.path(root, "figures"), pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
stopifnot(length(pngs) >= 10, all(file.exists(pdfs)),
          setequal(pdfs, all_pdfs),
          all(file.info(c(pngs, pdfs))$size > 1000))
for (path in pdfs) {
  connection <- file(path, "rb")
  header <- rawToChar(readBin(connection, "raw", n = 5))
  close(connection)
  stopifnot(identical(header, "%PDF-"))
}
dir.create(file.path(root, "provenance"), recursive = TRUE, showWarnings = FALSE)
writeLines(c("Output checks: PASS", paste("Spatial targets:", nrow(spatial)),
              paste("Bulk genes:", nrow(bulk)),
              paste("Adjusted CP-aware spatial overlap:", nrow(integrated)),
              paste("Jointly measurable spatial-associated genes:", nrow(spatial_overlap)),
              paste("Plasma models:", nrow(plasma)), paste("PLCO samples:", nrow(scores)),
              paste("PNG/PDF figure pairs:", length(pngs))),
            file.path(root, "provenance/output_checks.txt"))
cat("Output checks passed:", length(pngs), "PNG/PDF pairs.\n")
