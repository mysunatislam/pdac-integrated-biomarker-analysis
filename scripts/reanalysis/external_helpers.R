options(stringsAsFactors = FALSE)
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
external_dir <- "reanalysis/results/external"
dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)
external_path <- function(id) file.path("data/external", paste0(id, "_series_matrix.txt.gz"))

read_external_geo <- function(id, expression = TRUE, path = external_path(id)) {
  lines <- readLines(gzfile(path), warn = FALSE)
  parts <- strsplit(lines[startsWith(lines, "!Sample_")], "\t", fixed = TRUE)
  fields <- lapply(parts, function(x) gsub('^"|"$', "", x[-1]))
  names(fields) <- make.unique(vapply(parts, function(x) substring(x[1], 2), character(1)))
  stopifnot(all(lengths(fields) == length(fields$Sample_geo_accession)))
  meta <- as.data.frame(fields, check.names = FALSE)
  names(meta)[names(meta) == "Sample_geo_accession"] <- "GSM"
  names(meta)[names(meta) == "Sample_title"] <- "Title"
  result <- list(meta = meta, lines = lines[startsWith(lines, "!Series_")])
  if (expression) {
    start <- which(startsWith(lines, "!series_matrix_table_begin"))
    end <- which(startsWith(lines, "!series_matrix_table_end"))
    stopifnot(length(start) == 1, length(end) == 1)
    tab <- read.delim(text = paste(lines[seq.int(start + 1, end - 1)], collapse = "\n"),
                      check.names = FALSE, quote = "\"", na.strings = c("null", "NA"))
    if (nrow(tab) > 0 && ncol(tab) > 1) {
      x <- as.matrix(tab[, -1, drop = FALSE])
      storage.mode(x) <- "double"
      rownames(x) <- as.character(tab[[1]])
      stopifnot(!anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)), !any(is.infinite(x)))
      meta <- meta[match(colnames(x), meta$GSM), , drop = FALSE]
      stopifnot(identical(meta$GSM, colnames(x)))
      result$meta <- meta
      result$x <- x
    }
  }
  result
}

characteristic <- function(meta, prefix, required = TRUE) {
  cols <- grep("^Sample_characteristics_ch1", names(meta), value = TRUE)
  hit <- cols[vapply(meta[cols], function(x) any(startsWith(tolower(x), tolower(prefix))), logical(1))]
  if (!length(hit) && !required) return(rep(NA_character_, nrow(meta)))
  stopifnot(length(hit) == 1)
  value <- meta[[hit]]
  result <- substring(value, nchar(prefix) + 1)
  result[!startsWith(tolower(value), tolower(prefix))] <- NA_character_
  trimws(result)
}

collapse_probes <- function(x, symbols) {
  stopifnot(nrow(x) == length(symbols))
  valid <- !is.na(symbols) & nzchar(symbols) & !grepl("[/,;[:space:]]", symbols)
  x <- x[valid, , drop = FALSE]
  symbols <- symbols[valid]
  pick <- order(symbols, -rowMeans(x), rownames(x))
  pick <- pick[!duplicated(symbols[pick])]
  probes <- rownames(x)[pick]
  x <- x[pick, , drop = FALSE]
  rownames(x) <- symbols[pick]
  list(x = x, probes = probes)
}

limma_result <- function(fit, contrast, dataset, comparison) {
  effect <- fit$coefficients[, contrast]
  se <- fit$stdev.unscaled[, contrast] * sqrt(fit$s2.post)
  margin <- qt(0.975, fit$df.total) * se
  data.frame(Dataset = dataset, Comparison = comparison, Gene = rownames(fit),
             log2FC = effect, CI_low = effect - margin, CI_high = effect + margin,
             SE = se, P = fit$p.value[, contrast],
             BH_all_genes = p.adjust(fit$p.value[, contrast], "BH"))
}

attach_reference <- function(stats, spatial_contrast) {
  reference <- read.csv(file.path(external_dir, "frozen_spatial_reference.csv"))
  out <- merge(stats, reference, by = "Gene", sort = FALSE)
  out$BH_frozen_family <- p.adjust(out$P, "BH")
  out$Spatial_log2FC <- out[[paste0(spatial_contrast, "_log2")]]
  out$Spatial_BH <- out[[paste0(spatial_contrast, "_BH_FDR_all_targets")]]
  out$Same_direction <- sign(out$log2FC) == sign(out$Spatial_log2FC)
  out$Spatial_Contrast <- spatial_contrast
  out
}
