#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
out_dir <- file.path("reanalysis", "results", "bulk")
fig_dir <- file.path("reanalysis", "figures", "bulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

source("scripts/reanalysis/00_audit_inputs.R", local = TRUE)
path <- file.path(source_dir, "GSE143754_series_matrix.txt.gz")
lines <- readLines(gzfile(path), warn = FALSE)
begin <- which(startsWith(lines, "!series_matrix_table_begin"))
end <- which(startsWith(lines, "!series_matrix_table_end"))
stopifnot(length(begin) == 1, length(end) == 1, end > begin)
matrix_table <- read.delim(
  text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
  check.names = FALSE, quote = "\"", stringsAsFactors = FALSE
)
probe <- as.character(matrix_table[[1]])
expression <- as.matrix(matrix_table[, -1])
storage.mode(expression) <- "double"
rownames(expression) <- probe
stopifnot(!anyDuplicated(probe), all(is.finite(expression)))

metadata <- sample_table(series_fields(path))
metadata <- metadata[match(colnames(expression), metadata$GSM), ]
stopifnot(identical(metadata$GSM, colnames(expression)))
characteristics <- grep("^Sample_characteristics_ch1", names(metadata), value = TRUE)
disease_column <- characteristics[vapply(metadata[characteristics], function(x) {
  any(startsWith(x, "disease state:"))
}, logical(1))]
stopifnot(length(disease_column) == 1)
group <- sub("^disease state: ", "", metadata[[disease_column]])
group <- factor(group, levels = c("Adjacent Normal", "Chronic Pancreatitis", "Tumor"))
stopifnot(identical(as.integer(table(group)), c(9L, 6L, 11L)))
metadata$Group <- as.character(group)
write.csv(metadata[, c("GSM", "Title", "Group", characteristics)],
          file.path(out_dir, "GSE143754_sample_metadata.csv"), row.names = FALSE)

annotation_db <- hta20transcriptcluster.db::hta20transcriptcluster.db
known_probes <- AnnotationDbi::keys(annotation_db, keytype = "PROBEID")
mapped <- AnnotationDbi::mapIds(
  annotation_db, keys = intersect(probe, known_probes), column = "SYMBOL",
  keytype = "PROBEID", multiVals = "asNA"
)
symbol <- unname(mapped[probe])
probe_map <- data.frame(Probe = probe, Gene = symbol,
                        Annotation = paste0("hta20transcriptcluster.db ",
                                             packageVersion("hta20transcriptcluster.db")))
write.csv(probe_map, file.path(out_dir, "GSE143754_probe_annotation.csv"), row.names = FALSE)
valid <- !is.na(symbol) & nzchar(symbol)
expression <- expression[valid, , drop = FALSE]
symbol <- symbol[valid]
probe <- probe[valid]
mean_expression <- rowMeans(expression)
pick <- order(symbol, -mean_expression, probe)
keep <- pick[!duplicated(symbol[pick])]
expression <- expression[keep, , drop = FALSE]
symbol <- symbol[keep]
probe <- probe[keep]
stopifnot(!anyDuplicated(symbol))

age_column <- characteristics[vapply(metadata[characteristics], function(x) {
  all(startsWith(x, "age:"))
}, logical(1))]
sex_column <- characteristics[vapply(metadata[characteristics], function(x) {
  all(startsWith(x, "Sex:"))
}, logical(1))]
stopifnot(length(age_column) == 1, length(sex_column) == 1)
age <- as.numeric(sub("^age: ", "", metadata[[age_column]]))
sex <- factor(sub("^Sex: ", "", metadata[[sex_column]]))
stopifnot(all(is.finite(age)), nlevels(sex) == 2)
design <- model.matrix(~0 + group + age + sex)
colnames(design)[1:3] <- c("Normal", "CP", "PDAC")
stopifnot(qr(design)$rank == ncol(design))

# Exclude adjacent-normal samples from primary inference because pairing is unknown.
primary_index <- group %in% c("Chronic Pancreatitis", "Tumor")
primary_metadata <- data.frame(
  Group = factor(ifelse(group[primary_index] == "Tumor", "PDAC", "CP"),
                 levels = c("CP", "PDAC")),
  Age = age[primary_index], Sex = droplevels(sex[primary_index])
)
primary_design <- model.matrix(~0 + Group + Age + Sex, primary_metadata)
colnames(primary_design)[1:2] <- c("CP", "PDAC")
stopifnot(sum(primary_index) == 17, qr(primary_design)$rank == ncol(primary_design))
primary_contrast <- limma::makeContrasts(PDAC_vs_CP = PDAC - CP, levels = primary_design)
primary_fit <- limma::lmFit(expression[, primary_index, drop = FALSE], primary_design)
primary_fit <- limma::eBayes(limma::contrasts.fit(primary_fit, primary_contrast),
                            trend = TRUE, robust = TRUE)

unadjusted_design <- model.matrix(~0 + Group, primary_metadata)
colnames(unadjusted_design) <- c("CP", "PDAC")
unadjusted_fit <- limma::lmFit(expression[, primary_index, drop = FALSE], unadjusted_design)
unadjusted_contrasts <- limma::makeContrasts(PDAC_vs_CP = PDAC - CP,
                                             levels = unadjusted_design)
unadjusted_fit <- limma::eBayes(limma::contrasts.fit(unadjusted_fit, unadjusted_contrasts),
                               trend = TRUE, robust = TRUE)
fit <- limma::lmFit(expression, design)
contrasts <- limma::makeContrasts(
  PDAC_vs_CP = PDAC - CP, PDAC_vs_Normal = PDAC - Normal,
  CP_vs_Normal = CP - Normal, levels = design
)
fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend = TRUE, robust = TRUE)
contrast <- function(name, model = fit) {
  effect <- model$coefficients[, name]
  se <- model$stdev.unscaled[, name] * sqrt(model$s2.post)
  p <- model$p.value[, name]
  margin <- qt(0.975, model$df.total) * se
  list(effect = effect, p = p, fdr = p.adjust(p, method = "BH"),
        low = effect - margin, high = effect + margin)
}
pdac_cp <- contrast("PDAC_vs_CP", primary_fit)
pdac_cp_all_groups <- contrast("PDAC_vs_CP")
pdac_normal <- contrast("PDAC_vs_Normal")
cp_normal <- contrast("CP_vs_Normal")
stats <- data.frame(
  Gene = symbol, Probe = probe, Mean_Log2_Expression = mean_expression[keep],
  PDAC_vs_CP_log2FC = pdac_cp$effect, PDAC_vs_CP_p = pdac_cp$p,
  PDAC_vs_CP_BH_FDR = pdac_cp$fdr,
  PDAC_vs_CP_CI_low = pdac_cp$low, PDAC_vs_CP_CI_high = pdac_cp$high,
  PDAC_vs_Normal_log2FC = pdac_normal$effect,
  PDAC_vs_Normal_p = pdac_normal$p, PDAC_vs_Normal_BH_FDR = pdac_normal$fdr,
  PDAC_vs_Normal_CI_low = pdac_normal$low, PDAC_vs_Normal_CI_high = pdac_normal$high,
  CP_vs_Normal_log2FC = cp_normal$effect,
  CP_vs_Normal_p = cp_normal$p, CP_vs_Normal_BH_FDR = cp_normal$fdr,
  CP_vs_Normal_CI_low = cp_normal$low, CP_vs_Normal_CI_high = cp_normal$high
)
stats$Unadjusted_PDAC_vs_CP_log2FC <- as.numeric(unadjusted_fit$coefficients[, 1])
stats$Unadjusted_PDAC_vs_CP_BH_FDR <- p.adjust(unadjusted_fit$p.value[, 1], "BH")
stats$AllGroup_Adjusted_PDAC_vs_CP_log2FC <- pdac_cp_all_groups$effect
stats$AllGroup_Adjusted_PDAC_vs_CP_BH_FDR <- pdac_cp_all_groups$fdr
stats <- stats[order(stats$PDAC_vs_CP_BH_FDR, -abs(stats$PDAC_vs_CP_log2FC)), ]
write.csv(stats, file.path(out_dir, "GSE143754_gene_level_contrasts.csv"), row.names = FALSE)

cp_aware <- stats[stats$PDAC_vs_CP_BH_FDR < 0.05 &
                    abs(stats$PDAC_vs_CP_log2FC) >= 0.5, ]
write.csv(cp_aware, file.path(out_dir, "GSE143754_PDAC_vs_CP_candidates.csv"),
          row.names = FALSE)
spatial_features <- read.csv(file.path("reanalysis", "results", "spatial",
                                       "GSE208536_patient_blocked_all_targets.csv"))
overlap <- merge(cp_aware, spatial_features, by = "Gene")
write.csv(overlap, file.path(out_dir, "CP_aware_genes_measured_in_spatial_panel.csv"),
          row.names = FALSE)

volcano <- ggplot2::ggplot(stats, ggplot2::aes(
  PDAC_vs_CP_log2FC, -log10(pmax(PDAC_vs_CP_BH_FDR, .Machine$double.xmin)))) +
  ggplot2::geom_point(ggplot2::aes(color = PDAC_vs_CP_BH_FDR < 0.05),
                      size = 0.9, alpha = 0.7) +
  ggplot2::scale_color_manual(values = c("FALSE" = "#66727C", "TRUE" = "#B7404B"),
                              guide = "none") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#66727C") +
  ggplot2::labs(x = "PDAC minus chronic pancreatitis (log2 expression)",
                y = "-log10(BH FDR)", title = "PDAC versus chronic pancreatitis",
                subtitle = "GSE143754: 11 PDAC and 6 chronic pancreatitis tissues") +
  ggplot2::theme_classic(base_size = 11)
ggplot2::ggsave(file.path(fig_dir, "PDAC_vs_CP_volcano.png"), volcano,
                width = 6, height = 4.3, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "PDAC_vs_CP_volcano.pdf"), volcano,
                width = 6, height = 4.3)

cat("Annotated genes:", nrow(stats), "\n")
cat("PDAC vs CP at FDR < 0.05:", sum(stats$PDAC_vs_CP_BH_FDR < 0.05), "\n")
cat("PDAC vs CP with |log2FC| >= 0.5:", nrow(cp_aware), "\n")
cat("CP-aware candidates measured in GeoMx panel:", nrow(overlap), "\n")
cat("Primary PDAC-versus-CP model: 17 specimens, adjusted for age and sex.\n")
cat("All-group adjusted PDAC vs CP at FDR < 0.05:",
    sum(stats$AllGroup_Adjusted_PDAC_vs_CP_BH_FDR < 0.05), "\n")
