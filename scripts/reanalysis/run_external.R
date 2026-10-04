#!/usr/bin/env Rscript
steps <- c("12_download_external.R", "13_audit_external.R", "14_external_tissue.R",
           "15_external_serum.R", "16_external_figures_report.R", "17_validate_external.R",
           "09_provenance.R", "10_validate_outputs.R")
for (step in steps) {
  cat("\nRunning", step, "\n")
  source(file.path("scripts/reanalysis", step), local = new.env(parent = globalenv()))
}
cat("Independent public-cohort extension complete.\n")
