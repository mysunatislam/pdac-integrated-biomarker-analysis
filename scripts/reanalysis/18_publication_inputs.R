#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 600)
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
out_dir <- "data/publication"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sha <- "b6eac73e91c246fcff0bb1a5c68a816cd588fc48"
sources <- data.frame(
  File = c("MCPcounter_genes.txt", "MCPcounter_signature_license.txt",
           "E-MTAB-1791_metadata.json", "E-MTAB-1791.sdrf.txt"),
  URL = c(paste0("https://raw.githubusercontent.com/ebecht/MCPcounter/", sha, "/Signatures/genes.txt"),
          paste0("https://raw.githubusercontent.com/ebecht/MCPcounter/", sha, "/Signatures/License.txt"),
          "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-1791",
          "https://www.ebi.ac.uk/biostudies/files/E-MTAB-1791/E-MTAB-1791.sdrf.txt")
)
for (i in seq_len(nrow(sources))) {
  path <- file.path(out_dir, sources$File[i])
  if (!file.exists(path)) curl::curl_download(sources$URL[i], path, quiet = FALSE)
}
sources$MD5 <- unname(tools::md5sum(file.path(out_dir, sources$File)))
write.csv(sources, file.path(out_dir, "publication_download_manifest.csv"), row.names = FALSE)
stopifnot(jsonlite::fromJSON(file.path(out_dir, "E-MTAB-1791_metadata.json"))$accno == "E-MTAB-1791")
sdrf <- read.delim(file.path(out_dir, "E-MTAB-1791.sdrf.txt"), check.names = FALSE)
stopifnot(nrow(sdrf) == 457, !anyDuplicated(sdrf[["Source Name"]]),
          !anyDuplicated(sdrf[["Assay Name"]]))
meta <- data.frame(
  Specimen = sdrf[["Source Name"]], Assay = sdrf[["Assay Name"]],
  Disease = sdrf[["Characteristics[disease]"]],
  Clinical_History = sdrf[["Characteristics[clinical history]"]],
  File = sdrf[["Derived Array Data File"]],
  Archive_URL = sub("^ftp://", "https://", sdrf[["Comment [Derived ArrayExpress FTP file]"]])
)
meta$Chip <- sub("_[A-F]$", "", meta$Assay)
meta$Role <- "Excluded_other_context"
meta$Role[meta$Disease == "pancreatic ductal adenocarcinoma" &
            meta$Clinical_History == "pancreatic ductal adenocarcinoma tissue"] <- "PDAC"
pancreatitis <- meta$Disease == "pancreatitis"
meta$Role[pancreatitis] <- "Pancreatitis_associated_context"
meta$Role[pancreatitis & meta$Clinical_History ==
            "pancreatitis tissue from pancreas with pancreatitis"] <- "CP"
meta$Primary_Eligible <- meta$Role %in% c("PDAC", "CP")
stopifnot(sum(meta$Role == "PDAC") == 195, sum(meta$Role == "CP") == 59,
          sum(meta$Role == "Pancreatitis_associated_context") == 9)
dir.create("reanalysis/results/publication", recursive = TRUE, showWarnings = FALSE)
write.csv(meta, "reanalysis/results/publication/E-MTAB-1791_eligibility_mapping.csv", row.names = FALSE)
counts <- as.data.frame(table(meta$Disease, meta$Clinical_History))
counts <- counts[counts$Freq > 0, ]
names(counts) <- c("Disease", "Clinical_History", "Specimens")
write.csv(counts, "reanalysis/results/publication/E-MTAB-1791_context_counts.csv", row.names = FALSE)
print(counts, row.names = FALSE)
cat("Publication input hashes and cohort metadata saved.\n")
