#!/usr/bin/env Rscript

source("scripts/reanalysis/external_helpers.R")
options(timeout = 600, digits = 17)
args <- commandArgs(trailingOnly = TRUE)
python_arg <- grep("^--python=", args, value = TRUE)
python <- if (length(python_arg)) sub("^--python=", "", python_arg[1]) else Sys.which("python")
out <- "data/publication"
meta <- read.csv("reanalysis/results/publication/E-MTAB-1791_eligibility_mapping.csv")
meta <- meta[meta$Primary_Eligible, ]
stopifnot(nrow(meta) == 254, sum(meta$Role == "PDAC") == 195,
          sum(meta$Role == "CP") == 59, !anyDuplicated(meta$Specimen))
default_cache <- if (.Platform$OS.type == "windows")
  file.path(Sys.getenv("LOCALAPPDATA"), "Codex", "PDAC-public-data-cache") else
  file.path(path.expand("~"), ".cache", "pdac-public-data")
cache_root <- Sys.getenv("PDAC_DATA_CACHE", default_cache)
raw_dir <- file.path(cache_root, "emtab1791_arrays")
cache_dir <- file.path(cache_root, "emtab1791_cache")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
manifest_file <- file.path(out, "E-MTAB-1791_array_manifest.csv")
manifest <- if (file.exists(manifest_file)) read.csv(manifest_file) else data.frame()
if (nrow(manifest) && !"Archive_URL" %in% names(manifest)) {
  manifest$Archive_URL <- meta$Archive_URL[match(manifest$Specimen, meta$Specimen)]
  manifest$Retrieval <- "Uncompressed canonical file"
}
url_root <- "https://www.ebi.ac.uk/biostudies/files/E-MTAB-1791/"
file_records <- list()
collect_files <- function(node) {
  if (!is.list(node)) return(invisible(NULL))
  if (identical(node$type, "file") && !is.null(node$path) && !is.null(node$size)) {
    file_records[[length(file_records) + 1]] <<- data.frame(File = node$path, Bytes = node$size)
  } else for (item in node) collect_files(item)
}
collect_files(jsonlite::fromJSON(file.path(out, "E-MTAB-1791_metadata.json"), simplifyVector = FALSE))
sizes <- unique(do.call(rbind, file_records))
stopifnot(!anyDuplicated(sizes$File))
expected_bytes <- sizes$Bytes[match(meta$File, sizes$File)]
stopifnot(!anyNA(expected_bytes))

verified_raw <- function(path, row) {
  record_path <- paste0(path, ".verified.json")
  if (!file.exists(path) || !file.exists(record_path)) return(FALSE)
  tryCatch({
    record <- jsonlite::read_json(record_path, simplifyVector = TRUE)
    isTRUE(record$File == meta$File[row] && record$Bytes == expected_bytes[row] &&
             file.info(path)$size == expected_bytes[row] &&
             record$MD5 == unname(tools::md5sum(path)))
  }, error = function(e) FALSE)
}

parse_array <- function(path, reference = FALSE) {
  header <- names(read.delim(path, nrows = 0, check.names = FALSE))
  normalized <- grep("\\.mean\\.qnorm\\.log$", header, value = TRUE)
  detection <- grep("\\.(pvalue|p)$", header, value = TRUE)
  stopifnot(length(normalized) == 1, length(detection) == 1)
  columns <- c("ProbeID", "PROBE_ID", "ENTREZ_GENE_ID", "SYMBOL", normalized, detection)
  stopifnot(all(columns %in% header))
  classes <- setNames(rep("NULL", length(header)), header)
  classes[c("ProbeID", "PROBE_ID", "ENTREZ_GENE_ID", "SYMBOL")] <- "character"
  classes[c(normalized, detection)] <- "numeric"
  if (reference) {
    tab <- read.delim(path, check.names = FALSE, colClasses = classes, quote = "\"",
                      comment.char = "", na.strings = c("", "NA"))
  } else {
    # R's numeric conversion and missing-string rule preserve the reference parser exactly.
    tab <- data.table::fread(file = path, sep = "\t", select = columns,
                            colClasses = list(character = columns), data.table = FALSE,
                            nThread = 2, na.strings = c("", "NA"), strip.white = FALSE,
                            showProgress = FALSE)
    tab[[normalized]] <- as.numeric(tab[[normalized]])
    tab[[detection]] <- as.numeric(tab[[detection]])
    for (column in columns[1:4])
      tab[[column]][tab[[column]] %in% c("", "NA")] <- NA_character_
  }
  stopifnot(nrow(tab) > 40000, !anyDuplicated(tab$ProbeID),
            all(is.finite(tab[[normalized]])),
            all(is.na(tab[[detection]]) | (tab[[detection]] >= 0 & tab[[detection]] <= 1)))
  tab <- tab[order(tab$ProbeID), ]
  rownames(tab) <- NULL
  list(annotation = tab[, c("ProbeID", "PROBE_ID", "ENTREZ_GENE_ID", "SYMBOL")],
       expression = tab[[normalized]], detection = tab[[detection]],
       normalized_column = normalized, detection_column = detection)
}

# Bound temporary storage while allowing independent archive members to transfer.
batch_size <- 8
parser_verified <- FALSE
for (start in seq.int(1, nrow(meta), by = batch_size)) {
  ii <- seq.int(start, min(start + batch_size - 1, nrow(meta)))
  paths <- file.path(raw_dir, meta$File[ii])
  caches <- file.path(cache_dir, paste0(meta$Specimen[ii], ".rds"))
  cached_row <- if (nrow(manifest)) match(meta$Specimen[ii], manifest$Specimen) else rep(NA_integer_, length(ii))
  cached_md5 <- unname(tools::md5sum(caches))
  needs <- !file.exists(caches) | is.na(cached_row) | is.na(cached_md5)
  if (nrow(manifest)) needs <- needs | (!is.na(cached_row) & cached_md5 != manifest$Cache_MD5[cached_row])
  needs[is.na(needs)] <- TRUE
  retrieval <- rep("Uncompressed canonical file", length(ii))
  if (any(needs)) {
    complete <- vapply(seq_along(ii), function(j) verified_raw(paths[j], ii[j]), logical(1))
    for (j in which(complete)) retrieval[j] <-
      jsonlite::read_json(paste0(paths[j], ".verified.json"))$Retrieval
    jj <- which(needs & !complete)
    if (length(jj) && nzchar(python)) {
      status <- system2(python, c(shQuote("scripts/reanalysis/fetch_inflammatory_arrays.py"),
                                  "--output-dir", shQuote(raw_dir), shQuote(meta$Specimen[ii[jj]])))
      if (status != 0) stop("ZIP retrieval failed; install the pinned Python requirements or retry")
      stopifnot(all(vapply(jj, function(j) verified_raw(paths[j], ii[j]), logical(1))))
      retrieval[jj] <- "SDRF ZIP member (CRC checked)"
      jj <- integer()
    }
    fallback_jj <- jj
    for (attempt in seq_len(3)) {
      if (!length(jj)) break
      status <- try(curl::multi_download(paste0(url_root, meta$File[ii[jj]]), paths[jj],
                                        resume = FALSE, progress = TRUE, multiplex = FALSE,
                                        timeout = 180, connecttimeout = 30, http_version = 2), silent = TRUE)
      correct_size <- file.exists(paths[jj]) & !is.na(file.info(paths[jj])$size) &
        file.info(paths[jj])$size == expected_bytes[ii[jj]]
      failed <- if (inherits(status, "try-error")) seq_along(jj) else
        which(!correct_size | !status$status_code %in% c(200, 206))
      if (!length(failed)) break
      jj <- jj[failed]
      Sys.sleep(2 * attempt)
      if (attempt == 3) stop("Download failed for: ", paste(meta$Specimen[ii[jj]], collapse = ", "))
    }
    for (j in fallback_jj) {
      jsonlite::write_json(list(File = meta$File[ii[j]], Bytes = file.info(paths[j])$size,
                                MD5 = unname(tools::md5sum(paths[j])),
                                Retrieval = "Uncompressed canonical file"),
                           paste0(paths[j], ".verified.json"), auto_unbox = TRUE)
      stopifnot(verified_raw(paths[j], ii[j]))
    }
    for (j in which(needs)) {
      row <- ii[j]
      stopifnot(file.info(paths[j])$size == expected_bytes[row])
      parsed <- parse_array(paths[j])
      md5 <- unname(tools::md5sum(paths[j]))
      if (!parser_verified) {
        stopifnot(identical(parsed, parse_array(paths[j], reference = TRUE)))
        write.csv(data.frame(File = meta$File[row], Source_MD5 = md5,
                             Probe_Count = length(parsed$expression),
                             Data_Table_Version = as.character(packageVersion("data.table")),
                             Exact_Base_Parser_Equivalence = TRUE),
                  file.path(out, "E-MTAB-1791_parser_equivalence.csv"), row.names = FALSE)
        parser_verified <- TRUE
        cat("Fast/reference parser equivalence: PASS\n")
      }
      old <- if (nrow(manifest)) match(meta$Specimen[row], manifest$Specimen) else NA_integer_
      if (!is.na(old)) stopifnot(identical(md5, manifest$Source_MD5[old]))
      saveRDS(parsed, caches[j], compress = "gzip")
      entry <- data.frame(Specimen = meta$Specimen[row], Assay = meta$Assay[row],
                          File = meta$File[row], URL = paste0(url_root, meta$File[row]),
                          Source_Bytes = file.info(paths[j])$size, Source_MD5 = md5,
                          Probe_Count = length(parsed$expression),
                          Normalized_Column = parsed$normalized_column,
                          Cache_MD5 = unname(tools::md5sum(caches[j])),
                          Archive_URL = meta$Archive_URL[row], Retrieval = retrieval[j])
      if (!is.na(old)) manifest <- manifest[-old, , drop = FALSE]
      manifest <- rbind(manifest, entry)
      write.csv(manifest, manifest_file, row.names = FALSE)
      stopifnot(identical(dirname(normalizePath(paths[j], winslash = "/")),
                          normalizePath(raw_dir, winslash = "/")))
      unlink(paths[j], recursive = FALSE)
      unlink(paste0(paths[j], ".verified.json"), recursive = FALSE)
    }
  }
  cat(sprintf("Validated compact arrays: %d/%d\n", min(start + batch_size - 1, nrow(meta)), nrow(meta)))
  flush.console()
}

manifest <- manifest[match(meta$Specimen, manifest$Specimen), ]
stopifnot(!anyNA(manifest$Specimen), identical(meta$Specimen, manifest$Specimen))
cache_paths <- file.path(cache_dir, paste0(meta$Specimen, ".rds"))
stopifnot(identical(unname(tools::md5sum(cache_paths)), manifest$Cache_MD5))
first <- readRDS(cache_paths[1])
annotation <- first$annotation
x <- matrix(NA_real_, nrow(annotation), nrow(meta),
            dimnames = list(annotation$ProbeID, meta$Specimen))
qc <- list()
for (i in seq_len(nrow(meta))) {
  array <- readRDS(cache_paths[i])
  stopifnot(identical(array$annotation, annotation))
  x[, i] <- array$expression
  qc[[i]] <- data.frame(Specimen = meta$Specimen[i], Probes = nrow(annotation),
                       Median_log2 = median(array$expression),
                       IQR_log2 = IQR(array$expression),
                       Fraction_detection_P_below_0_01 = mean(array$detection < 0.01, na.rm = TRUE))
}
stopifnot(all(is.finite(x)))
write.csv(annotation, file.path(out, "E-MTAB-1791_probe_annotation.csv"), row.names = FALSE)
probe_file <- file.path(out, "E-MTAB-1791_primary_normalized_probes.csv.gz")
con <- gzfile(probe_file, "wt")
write.csv(data.frame(ProbeID = rownames(x), x, check.names = FALSE), con, row.names = FALSE)
close(con)

db <- org.Hs.eg.db::org.Hs.eg.db
entrez <- annotation$ENTREZ_GENE_ID
known <- intersect(entrez, AnnotationDbi::keys(db, keytype = "ENTREZID"))
symbols <- AnnotationDbi::mapIds(db, keys = known, keytype = "ENTREZID",
                               column = "SYMBOL", multiVals = "asNA")
symbol <- unname(symbols[entrez])
gene <- collapse_probes(x, symbol)
probe_index <- match(gene$probes, annotation$ProbeID)
gene_map <- data.frame(Gene = rownames(gene$x), ProbeID = gene$probes,
                       ILMN_Probe = annotation$PROBE_ID[probe_index],
                       Entrez_ID = entrez[probe_index],
                       Deposited_Symbol = annotation$SYMBOL[probe_index])
write.csv(gene_map, file.path(out, "E-MTAB-1791_selected_probe_map.csv"), row.names = FALSE)
con <- gzfile(file.path(out, "E-MTAB-1791_primary_normalized_genes.csv.gz"), "wt")
write.csv(data.frame(Gene = rownames(gene$x), gene$x, check.names = FALSE), con, row.names = FALSE)
close(con)
write.csv(do.call(rbind, qc), "reanalysis/results/publication/E-MTAB-1791_array_qc.csv", row.names = FALSE)
write.csv(manifest, manifest_file, row.names = FALSE)
stopifnot("STAT1" %in% rownames(gene$x), all(file.info(c(probe_file,
             file.path(out, "E-MTAB-1791_primary_normalized_genes.csv.gz")))$size < 100 * 1024^2))
cat(sprintf("Prepared %d probes and %d genes for 195 PDAC / 59 CP specimens.\n",
            nrow(x), nrow(gene$x)))
