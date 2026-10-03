#!/usr/bin/env Rscript

stopifnot(file.exists("data/source/GSE208536_Processed_data.xlsx"))
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
scripts <- c("00_audit_inputs.R", "01_spatial_patient_level.R", "02_ev_discovery.R",
              "03_plasma_plco_validation.R", "04_bulk_cp_aware.R", "05_secondary_ev_cohort.R",
              "06_integrated_tissue.R", "07_panel_enrichment.R", "08_qc_and_calibration.R",
              "09_provenance.R", "10_validate_outputs.R", "11_write_report.R")
for (script in scripts) {
  cat("\nRunning", script, "\n")
  source(file.path("scripts", "reanalysis", script), local = new.env(parent = globalenv()))
}
