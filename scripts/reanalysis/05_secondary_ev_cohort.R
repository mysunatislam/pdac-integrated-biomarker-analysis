#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source("scripts/reanalysis/00_audit_inputs.R", local = TRUE)
out_dir <- file.path("reanalysis", "results", "secondary_ev")
fig_dir <- file.path("reanalysis", "figures", "secondary_ev")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rpm <- read.delim(gzfile(file.path(source_dir, "GSE268771_miRNA_expression_RPM.txt.gz")),
                  check.names = FALSE)
feature <- sub("^hsa-", "", as.character(rpm[[1]]))
expression <- as.matrix(rpm[, -1])
storage.mode(expression) <- "double"
rownames(expression) <- feature
stopifnot(!anyDuplicated(feature), all(is.finite(expression)), all(expression >= 0))

metadata <- secondary_meta
metadata$Expression_ID <- sub(".*, ", "", metadata$Title)
stopifnot(!anyDuplicated(metadata$Expression_ID),
          setequal(metadata$Expression_ID, colnames(expression)))
metadata <- metadata[match(colnames(expression), metadata$Expression_ID), ]
stopifnot(identical(metadata$Expression_ID, colnames(expression)))
metadata$Diagnosis <- sub("^diagnosis: ", "", metadata[["Sample_characteristics_ch1.1"]])
metadata$Stage <- sub("^Stage: ", "", metadata[["Sample_characteristics_ch1.2"]])
metadata$Palliative <- sub("^treatment \\(palliative\\): ", "",
                           metadata[["Sample_characteristics_ch1.3"]])
stopifnot(all(as.integer(table(metadata$Diagnosis)[c("PDAC", "benign",
                                                    "healthy control")]) ==
              c(51L, 12L, 3L)))
write.csv(metadata[, c("GSM", "Expression_ID", "Diagnosis", "Stage", "Palliative")],
          file.path(out_dir, "GSE268771_sample_map.csv"), row.names = FALSE)

value <- log2(expression + 1)
pdac <- metadata$Diagnosis == "PDAC"
benign <- metadata$Diagnosis == "benign"
all_controls <- !pdac
test_feature <- function(miRNA) {
  x <- value[miRNA, ]
  test_p <- function(a, b) {
    if (length(unique(c(a, b))) < 2) return(1)
    suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  }
  p_all <- test_p(x[pdac], x[all_controls])
  p_benign <- test_p(x[pdac], x[benign])
  data.frame(
    miRNA = miRNA,
    PDAC_minus_all_controls_log2RPM = median(x[pdac]) - median(x[all_controls]),
    PDAC_minus_benign_log2RPM = median(x[pdac]) - median(x[benign]),
    PDAC_vs_all_controls_p = p_all,
    PDAC_vs_benign_p = p_benign
  )
}
stats <- do.call(rbind, lapply(rownames(value), test_feature))
stats$PDAC_vs_all_controls_BH_FDR <- p.adjust(stats$PDAC_vs_all_controls_p, "BH")
stats$PDAC_vs_benign_BH_FDR <- p.adjust(stats$PDAC_vs_benign_p, "BH")
write.csv(stats, file.path(out_dir, "GSE268771_exploratory_differential.csv"),
          row.names = FALSE)

historical6 <- c("miR-184", "miR-107", "miR-216b-5p", "miR-20a-5p",
                 "miR-134-3p", "miR-215-3p")
published3 <- c("let-7i-5p", "miR-130a-3p", "miR-221-3p")
markers <- c(historical6, published3)
availability <- data.frame(miRNA = markers, Available = markers %in% rownames(value))
write.csv(availability, file.path(out_dir, "signature_feature_availability.csv"),
          row.names = FALSE)
available <- intersect(markers, rownames(value))
marker_stats <- stats[match(available, stats$miRNA), ]
marker_stats$AUC_PDAC_vs_all_controls <- vapply(available, function(miRNA) {
  as.numeric(pROC::auc(pROC::roc(as.integer(pdac), value[miRNA, ],
                                levels = c(0, 1), direction = "<", quiet = TRUE)))
}, numeric(1))
write.csv(marker_stats, file.path(out_dir, "signature_marker_descriptive.csv"),
          row.names = FALSE)

composition <- data.frame(
  Group = c("Healthy", "Benign", "PDAC I", "PDAC II", "PDAC III", "PDAC IV"),
  Samples = c(sum(metadata$Diagnosis == "healthy control"),
              sum(benign), sapply(as.character(1:4), function(stage) {
                sum(pdac & metadata$Stage == stage)
              }))
)
write.csv(composition, file.path(out_dir, "cohort_composition.csv"), row.names = FALSE)
draw_composition <- function() {
  barplot(composition$Samples, names.arg = composition$Group, las = 2,
          col = c("#237A77", "#6A9370", "#D3A44C", "#C86B31",
                  "#AD4C4C", "#873D62"),
          ylab = "Number of participants", main = "GSE268771 cohort composition",
          ylim = c(0, max(composition$Samples) + 5))
}
png(file.path(fig_dir, "cohort_composition.png"), width = 1800,
    height = 1300, res = 300)
par(mar = c(6, 4, 3, 1))
draw_composition()
dev.off()
pdf(file.path(fig_dir, "cohort_composition.pdf"), width = 6, height = 4.3)
par(mar = c(6, 4, 3, 1))
draw_composition()
dev.off()

print(composition, row.names = FALSE)
cat("Measured miRNAs:", nrow(stats), "\n")
cat("Markers available:", length(available), "of", length(markers), "\n")
cat("PDAC versus benign FDR < 0.05:",
    sum(stats$PDAC_vs_benign_BH_FDR < 0.05, na.rm = TRUE), "\n")
