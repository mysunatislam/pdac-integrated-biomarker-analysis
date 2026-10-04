#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
dir.create(".R-library", showWarnings = FALSE)
.libPaths(c(normalizePath(".R-library"), .libPaths()))
cran <- c("BiocManager", "renv", "readxl", "ggplot2", "pROC", "nlme", "scales", "statmod", "glmnet")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = ".R-library")
bioc <- c("limma", "edgeR", "AnnotationDbi", "hta20transcriptcluster.db",
          "hgu133plus2.db", "org.Hs.eg.db", "GO.db")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  BiocManager::install(missing, lib = ".R-library", ask = FALSE, update = FALSE)
}
