#!/usr/bin/env Rscript

if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
out_dir <- "reanalysis/provenance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dependencies <- c("readxl", "ggplot2", "pROC", "nlme", "scales", "limma",
                   "statmod", "AnnotationDbi", "hta20transcriptcluster.db", "org.Hs.eg.db",
                   "GO.db", "renv", "BiocManager", "edgeR", "hgu133plus2.db", "glmnet",
                   "curl", "jsonlite", "remotes", "MCPcounter", "data.table")
stopifnot(all(vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)))
session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_lines, file.path(out_dir, "sessionInfo.txt"))
versions <- data.frame(Package = dependencies, Version = vapply(dependencies,
                         function(p) as.character(packageVersion(p)), character(1)))
write.csv(versions, file.path(out_dir, "package_versions.csv"), row.names = FALSE)
external_files <- if (dir.exists("data/external")) list.files("data/external", full.names = TRUE) else character()
external_files <- external_files[!grepl("family.soft.gz$|\\.part$", external_files)]
publication_files <- if (dir.exists("data/publication"))
  list.files("data/publication", full.names = TRUE) else character()
publication_files <- publication_files[!file.info(publication_files)$isdir]
inputs <- c(list.files("data/source", full.names = TRUE), external_files, publication_files,
             "results/tables/spatial/valid_genes.txt")
hashes <- tools::md5sum(inputs)
write.csv(data.frame(File = names(hashes), MD5 = unname(hashes)),
          file.path(out_dir, "input_md5.csv"), row.names = FALSE)
renv::snapshot(project = getwd(), library = .libPaths(), packages = dependencies,
                lockfile = "renv.lock", prompt = FALSE, force = TRUE)
cat("Session information, input checksums, and dependency lockfile saved.\n")
