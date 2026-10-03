#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
source_dir <- file.path("data", "source")
output_dir <- file.path("reanalysis", "results", "input_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

series_fields <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  parts <- strsplit(lines[startsWith(lines, "!Sample_")], "\t", fixed = TRUE)
  keys <- vapply(parts, function(x) substring(x[[1]], 2), character(1))
  values <- lapply(parts, function(x) gsub('^"|"$', "", x[-1]))
  names(values) <- make.unique(keys)
  values
}

sample_table <- function(fields) {
  n <- length(fields[["Sample_geo_accession"]])
  stopifnot(n > 0, all(lengths(fields) == n))
  out <- data.frame(GSM = fields[["Sample_geo_accession"]], Title = fields[["Sample_title"]])
  for (name in setdiff(names(fields), c("Sample_geo_accession", "Sample_title"))) {
    out[[name]] <- fields[[name]]
  }
  out
}

tma_key <- function(x) {
  digits <- regmatches(x, gregexpr("[0-9]+", x))
  vapply(digits, function(d) {
    if (length(d) < 2) return(NA_character_)
    sprintf("TMA%d_%03d", as.integer(d[[1]]), as.integer(d[[2]]))
  }, character(1))
}

geo_path <- function(accession) file.path(source_dir, paste0(accession, "_series_matrix.txt.gz"))

spatial <- sample_table(series_fields(geo_path("GSE208536")))
spatial_xlsx <- readxl::read_excel(file.path(source_dir, "GSE208536_Processed_data.xlsx"), sheet = "Sheet1")
stopifnot(nrow(spatial) == 48, ncol(spatial_xlsx) == 49)
description_columns <- grep("^Sample_description", names(spatial), value = TRUE)
stopifnot(length(description_columns) >= 2)
spatial$TMA_Key <- tma_key(spatial[[description_columns[[2]]]])
expr_columns <- names(spatial_xlsx)[-1]
expr_keys <- tma_key(expr_columns)
stopifnot(!anyNA(spatial$TMA_Key), !anyNA(expr_keys),
          !anyDuplicated(spatial$TMA_Key), !anyDuplicated(expr_keys),
          setequal(spatial$TMA_Key, expr_keys))
spatial <- spatial[match(expr_keys, spatial$TMA_Key), ]
stopifnot(identical(spatial$TMA_Key, expr_keys))
spatial$Expression_Column <- expr_columns
characteristic_columns <- grep("^Sample_characteristics_ch1", names(spatial), value = TRUE)
stopifnot(length(characteristic_columns) == 5)
spatial$Clinical_Profile <- do.call(paste, c(spatial[characteristic_columns[-1]], sep = " | "))
spatial$Condition <- ifelse(grepl("^Adjacent normal", spatial$Title, ignore.case = TRUE), "Normal",
                            ifelse(grepl("^ADM", spatial$Title, ignore.case = TRUE), "ADM", "PDAC"))
stopifnot(all(table(spatial$Condition) == 16))
write.csv(spatial[, c("GSM", "Title", characteristic_columns, "TMA_Key",
                      "Expression_Column", "Clinical_Profile", "Condition")],
          file.path(output_dir, "GSE208536_roi_mapping.csv"), row.names = FALSE)
profile_counts <- as.data.frame.matrix(table(spatial$Clinical_Profile, spatial$Condition))
profile_counts$Clinical_Profile <- rownames(profile_counts)
profile_counts <- profile_counts[, c("Clinical_Profile", "Normal", "ADM", "PDAC")]
write.csv(profile_counts, file.path(output_dir, "GSE208536_clinical_profile_counts.csv"), row.names = FALSE)

read_count_book <- function(path) {
  book <- suppressMessages(readxl::read_excel(path, sheet = "Raw", .name_repair = "minimal"))
  sample_ids <- as.character(unlist(book[1, -1], use.names = FALSE))
  features <- as.character(book[[1]][-(1:2)])
  stopifnot(!anyNA(sample_ids), !anyDuplicated(sample_ids), !anyDuplicated(features))
  list(sample_ids = sample_ids, features = features, n_features = length(features))
}

diagnostic <- read_count_book(file.path(source_dir, "GSE259327_Raw_counts_PDAC_203samples_new.xlsx"))
plco <- read_count_book(file.path(source_dir, "GSE259327_Raw_counts_PLCO_96samples_new.xlsx"))
plasma_meta <- sample_table(series_fields(geo_path("GSE259327")))
stopifnot(length(diagnostic$sample_ids) == 203, length(plco$sample_ids) == 96,
          nrow(plasma_meta) == 299,
          setequal(c(diagnostic$sample_ids, plco$sample_ids), plasma_meta$Title))
plasma_meta$Source_Book <- ifelse(plasma_meta$Title %in% diagnostic$sample_ids, "Diagnostic", "PLCO")
write.csv(plasma_meta[, c("GSM", "Title", grep("^Sample_characteristics_ch1",
                                             names(plasma_meta), value = TRUE), "Source_Book")],
          file.path(output_dir, "GSE259327_sample_metadata.csv"), row.names = FALSE)
ev_meta <- sample_table(series_fields(geo_path("GSE304572")))
write.csv(ev_meta[, c("GSM", "Title", grep("^Sample_characteristics_ch1",
                                          names(ev_meta), value = TRUE))],
          file.path(output_dir, "GSE304572_sample_metadata.csv"), row.names = FALSE)
secondary_meta <- sample_table(series_fields(geo_path("GSE268771")))
write.csv(secondary_meta[, c("GSM", "Title", grep("^Sample_characteristics_ch1",
                                                 names(secondary_meta), value = TRUE))],
          file.path(output_dir, "GSE268771_sample_metadata.csv"), row.names = FALSE)

cohorts <- data.frame(
  Dataset = c("GSE143754", "GSE208536", "GSE259327 diagnostic", "GSE259327 PLCO", "GSE304572", "GSE268771"),
  Samples = c(length(series_fields(geo_path("GSE143754"))[["Sample_geo_accession"]]),
              nrow(spatial), length(diagnostic$sample_ids), length(plco$sample_ids),
              length(series_fields(geo_path("GSE304572"))[["Sample_geo_accession"]]),
              length(series_fields(geo_path("GSE268771"))[["Sample_geo_accession"]]))
)
write.csv(cohorts, file.path(output_dir, "cohort_sizes.csv"), row.names = FALSE)
cat("Spatial ROI map: 48 of 48 matched by TMA ID\n")
cat("Distinct spatial clinical profiles:", nrow(profile_counts), "\n")
print(profile_counts, row.names = FALSE)
cat("Diagnostic and PLCO count books matched to GEO titles\n")
cat("Diagnostic features:", diagnostic$n_features, "PLCO features:", plco$n_features,
    "shared:", length(intersect(diagnostic$features, plco$features)), "\n")
print(cohorts, row.names = FALSE)
