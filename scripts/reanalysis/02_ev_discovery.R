#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source("scripts/reanalysis/00_audit_inputs.R", local = TRUE)
out_dir <- file.path("reanalysis", "results", "ev")
fig_dir <- file.path("reanalysis", "figures", "ev")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

book <- readxl::read_excel(file.path(source_dir, "GSE304572_processed_data.xlsx"),
                           sheet = "counts")
features <- as.character(book[[1]])
counts <- as.matrix(book[, -1])
storage.mode(counts) <- "double"
rownames(counts) <- features
stopifnot(!anyDuplicated(features), all(is.finite(counts)), all(counts >= 0))

metadata <- ev_meta
metadata$Expression_ID <- sub("-.*$", "", sub(".*Patient: ", "", metadata$Title))
stopifnot(!anyDuplicated(metadata$Expression_ID),
          setequal(metadata$Expression_ID, colnames(counts)))
metadata <- metadata[match(colnames(counts), metadata$Expression_ID), ]
stopifnot(identical(metadata$Expression_ID, colnames(counts)))
metadata$Group <- ifelse(grepl("^PDAC", metadata$Title, ignore.case = TRUE), "PDAC",
                         ifelse(grepl("^chronic pancreatitis", metadata$Title,
                                      ignore.case = TRUE), "CP", "IPMN"))
stopifnot(all(as.integer(table(metadata$Group)[c("CP", "IPMN", "PDAC")]) ==
              c(10L, 10L, 65L)))
write.csv(metadata[, c("GSM", "Title", "Expression_ID", "Group")],
          file.path(out_dir, "GSE304572_sample_map.csv"), row.names = FALSE)

is_mirna <- grepl("^(hsa-)?(miR-|let-)", features, ignore.case = TRUE)
counts <- counts[is_mirna, , drop = FALSE]
col_totals <- colSums(counts)
stopifnot(all(col_totals > 0))
log_cpm <- log2(sweep(counts, 2, col_totals, "/") * 1e6 + 1)
detected <- rowSums(log_cpm >= log2(2)) >= 10
log_cpm <- log_cpm[detected, , drop = FALSE]
case <- metadata$Group == "PDAC"
design <- model.matrix(~case)
voom <- limma::voom(counts[detected, , drop = FALSE], design, lib.size = col_totals)
voom_fit <- limma::eBayes(limma::lmFit(voom, design), robust = TRUE)
voom_table <- limma::topTable(voom_fit, coef = 2, number = Inf, sort.by = "P")
voom_table$miRNA <- rownames(voom_table)
write.csv(voom_table, file.path(out_dir, "GSE304572_voom_sensitivity.csv"), row.names = FALSE)

test_gene <- function(feature) {
  x <- log_cpm[feature, ]
  data.frame(
    miRNA = feature,
    Median_PDAC = median(x[case]),
    Median_Benign = median(x[!case]),
    Median_Difference_log2CPM = median(x[case]) - median(x[!case]),
    Wilcoxon_p = suppressWarnings(wilcox.test(x[case], x[!case], exact = FALSE)$p.value),
    Detectable_Samples = sum(x >= log2(2))
  )
}
stats <- do.call(rbind, lapply(rownames(log_cpm), test_gene))
stats$BH_FDR <- p.adjust(stats$Wilcoxon_p, method = "BH")
stats <- stats[order(stats$BH_FDR, -abs(stats$Median_Difference_log2CPM)), ]
write.csv(stats, file.path(out_dir, "GSE304572_EV_differential.csv"), row.names = FALSE)

plasma_features <- intersect(diagnostic$features, plco$features)
eligible <- stats[stats$miRNA %in% plasma_features, ]
candidate20 <- head(eligible, 20)
stopifnot(nrow(candidate20) == 20)
write.csv(candidate20, file.path(out_dir, "GSE304572_EV_candidate20.csv"), row.names = FALSE)

volcano <- ggplot2::ggplot(stats, ggplot2::aes(Median_Difference_log2CPM,
                                               -log10(pmax(BH_FDR, .Machine$double.xmin)))) +
  ggplot2::geom_point(ggplot2::aes(color = BH_FDR < 0.05), size = 1.3, alpha = 0.8) +
  ggplot2::scale_color_manual(values = c("FALSE" = "#66727C", "TRUE" = "#B7404B"),
                              guide = "none") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#66727C") +
  ggplot2::labs(x = "Median difference, PDAC minus benign (log2 CPM)",
                y = "-log10(BH FDR)", title = "Plasma EV miRNA discovery",
                subtitle = "GSE304572: 65 PDAC, 10 CP, 10 IPMN") +
  ggplot2::theme_classic(base_size = 11)
ggplot2::ggsave(file.path(fig_dir, "EV_differential.png"), volcano,
                width = 6, height = 4.3, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "EV_differential.pdf"), volcano,
                width = 6, height = 4.3)

cat("Measured miRNAs:", sum(is_mirna), "Detectable:", nrow(stats), "\n")
cat("EV miRNAs at FDR < 0.05:", sum(stats$BH_FDR < 0.05), "\n")
cat("Voom sensitivity miRNAs at FDR < 0.05:", sum(voom_table$adj.P.Val < 0.05), "\n")
cat("Eligible for diagnostic and PLCO:", nrow(eligible), "\n")
print(candidate20[, c("miRNA", "Median_Difference_log2CPM", "BH_FDR")],
      row.names = FALSE)
