#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 900)
out_dir <- "data/external"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("reanalysis/results/external", recursive = TRUE, showWarnings = FALSE)
registry <- data.frame(
  Dataset = c("GSE91035", "GSE15471", "GSE179248", "GSE179248",
              "GSE59856", "GSE85589", "GSE91035", "GSE101462"),
  File = c("GSE91035_series_matrix.txt.gz", "GSE15471_series_matrix.txt.gz",
           "GSE179248_series_matrix.txt.gz", "GSE179248_CountReads_D6_0_processed.xlsx",
           "GSE59856_series_matrix.txt.gz", "GSE85589_series_matrix.txt.gz",
           "GSE91035_normalized_data_with_gene_symbol.txt.gz",
           "GSE101462_series_matrix.txt.gz"),
  Directory = c("matrix", "matrix", "matrix", "suppl", "matrix", "matrix", "suppl", "matrix")
)
registry$URL <- vapply(seq_len(nrow(registry)), function(i) {
  group <- sub("[0-9]{3}$", "nnn", registry$Dataset[i])
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/%s/%s",
          group, registry$Dataset[i], registry$Directory[i], registry$File[i])
}, character(1))
for (i in seq_len(nrow(registry))) {
  path <- file.path(out_dir, registry$File[i])
  if (!file.exists(path)) {
    cat("Downloading", registry$File[i], "\n")
    temporary <- paste0(path, ".part")
    download.file(registry$URL[i], temporary, mode = "wb", method = "libcurl")
    stopifnot(file.info(temporary)$size > 100)
    stopifnot(file.rename(temporary, path))
  }
}
registry$MD5 <- unname(tools::md5sum(file.path(out_dir, registry$File)))
for (id in c("GPL18941", "GPL22763")) {
  file <- paste0(id, "_family.soft.gz")
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s",
                 sub("[0-9]{3}$", "nnn", id), id, file)
  platform_path <- file.path(out_dir, file)
  if (!file.exists(platform_path)) {
    download.file(url, paste0(platform_path, ".part"), mode = "wb", method = "libcurl")
    stopifnot(file.rename(paste0(platform_path, ".part"), platform_path))
  }
  platform_connection <- gzfile(platform_path, "rt")
  platform_lines <- character()
  repeat {
    block <- readLines(platform_connection, n = 2000, warn = FALSE)
    if (!length(block)) stop("Platform table missing")
    platform_lines <- c(platform_lines, block)
    if (any(startsWith(block, "!platform_table_end"))) break
  }
  close(platform_connection)
  begin <- which(startsWith(platform_lines, "!platform_table_begin"))
  end <- which(startsWith(platform_lines, "!platform_table_end"))[1]
  stopifnot(length(begin) == 1)
  platform <- read.delim(text = paste(platform_lines[seq.int(begin + 1, end - 1)], collapse = "\n"),
                        quote = "\"", check.names = FALSE)
  annotation_file <- paste0(id, "_annotation.csv")
  annotation_path <- file.path(out_dir, annotation_file)
  if (!file.exists(annotation_path)) {
    write.csv(platform, annotation_path, row.names = FALSE)
  } else {
    deposited <- read.csv(annotation_path, check.names = FALSE)
    stopifnot(isTRUE(all.equal(platform, deposited, check.attributes = FALSE)))
  }
  registry <- rbind(registry, data.frame(
    Dataset = id, File = c(file, annotation_file),
    Directory = c("soft (download cache, not committed)", "derived platform table"),
    URL = url, MD5 = unname(tools::md5sum(c(platform_path, file.path(out_dir, annotation_file))))
  ))
}
write.csv(registry, "reanalysis/results/external/download_manifest.csv", row.names = FALSE)

reference <- read.csv("reanalysis/results/spatial/GSE208536_patient_blocked_all_targets.csv")
reference <- reference[reference$Patient_Blocked_F_BH_FDR_all_targets < 0.05, ]
stopifnot(nrow(reference) == 506, !anyDuplicated(reference$Gene))
reference_path <- "reanalysis/results/external/frozen_spatial_reference.csv"
if (file.exists(reference_path)) {
  stopifnot(isTRUE(all.equal(read.csv(reference_path), reference, check.attributes = FALSE)))
} else {
  write.csv(reference, reference_path, row.names = FALSE)
}
cat("Public source files and the 506-gene spatial reference recorded.\n")
