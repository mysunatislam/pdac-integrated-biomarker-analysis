#!/usr/bin/env Rscript
source("scripts/reanalysis/external_helpers.R")
ids <- c("GSE91035", "GSE15471", "GSE179248", "GSE59856", "GSE85589", "GSE101462")
for (id in ids) {
  z <- read_external_geo(id)
  write.csv(z$meta, file.path(external_dir, paste0(id, "_deposited_metadata.csv")), row.names = FALSE)
  cat("\n", id, "samples", nrow(z$meta), "expression", dim(z$x), "\n")
  cat("Sample titles:\n")
  print(head(z$meta$Title, 12))
  cols <- grep("^Sample_characteristics_ch1|^Sample_source_name_ch1|^Sample_data_processing",
               names(z$meta), value = TRUE)
  for (col in cols) {
    cat(col, "\n")
    print(table(z$meta[[col]]), max = 25)
  }
  if (!is.null(z$x)) {
    print(summary(as.vector(z$x)))
    cat("Missing values:", sum(is.na(z$x)), "all-missing rows:",
        sum(rowSums(!is.na(z$x)) == 0), "\n")
    print(head(rownames(z$x), 12))
  }
}
cat("\nGSE91035 deposited annotation columns:\n")
print(read.delim(gzfile("data/external/GSE91035_normalized_data_with_gene_symbol.txt.gz"),
                 nrows = 3, check.names = FALSE)[, 1:8])
book_path <- "data/external/GSE179248_CountReads_D6_0_processed.xlsx"
cat("\nADM workbook sheets:\n")
print(readxl::excel_sheets(book_path))
for (sheet in readxl::excel_sheets(book_path)) {
  b <- readxl::read_excel(book_path, sheet = sheet, n_max = 4, .name_repair = "minimal")
  cat(sheet, "\n")
  print(names(b))
  print(as.data.frame(b[, seq_len(min(8, ncol(b)))]))
}
