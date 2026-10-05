#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
args <- commandArgs(trailingOnly = TRUE)
run <- function(file, extra = character()) {
  status <- system2(rscript, c(shQuote(file.path("scripts/reanalysis", file)), shQuote(extra)))
  if (status != 0) stop("Failed: ", file)
}
stopifnot(file.exists("reanalysis/results/external/frozen_spatial_reference.csv"))
run("18_publication_inputs.R")
run("19_spatial_composition_audit.R")
prepared <- c("data/publication/E-MTAB-1791_primary_normalized_genes.csv.gz",
               "data/publication/E-MTAB-1791_primary_normalized_probes.csv.gz")
if ("--download" %in% args || !all(file.exists(prepared)))
  run("20_prepare_inflammatory_cohort.R", grep("^--python=", args, value = TRUE))
run("21_inflammatory_cohort_analysis.R")
run("22_publication_figures_report.R")
run("23_validate_publication.R")
run("09_provenance.R")
cat("Publication audit and inflammatory-cohort workflow complete.\n")
