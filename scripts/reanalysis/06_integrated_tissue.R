#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source("scripts/reanalysis/00_audit_inputs.R", local = TRUE)
out_dir <- file.path("reanalysis", "results", "integrated_tissue")
fig_dir <- file.path("reanalysis", "figures", "integrated_tissue")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

bulk <- read.csv(file.path("reanalysis", "results", "bulk",
                           "GSE143754_gene_level_contrasts.csv"))
spatial <- read.csv(file.path("reanalysis", "results", "spatial",
                              "GSE208536_patient_blocked_all_targets.csv"))
integrated <- merge(bulk, spatial, by = "Gene")
integrated$CP_Aware <- integrated$PDAC_vs_CP_BH_FDR < 0.05 &
  abs(integrated$PDAC_vs_CP_log2FC) >= 0.5
integrated$CP_Aware_Unadjusted <- integrated$Unadjusted_PDAC_vs_CP_BH_FDR < 0.05 &
  abs(integrated$Unadjusted_PDAC_vs_CP_log2FC) >= 0.5
integrated$Spatial_Associated <- integrated$Patient_Blocked_F_BH_FDR_all_targets < 0.05
integrated$Normal_to_ADM_Associated <-
  integrated$Normal_to_ADM_BH_FDR_all_targets < 0.05
integrated$ADM_to_PDAC_Associated <-
  integrated$ADM_to_PDAC_BH_FDR_all_targets < 0.05

confirmed <- integrated[integrated$Spatial_Associated, ]
confirmed$Stage_Pattern <- ifelse(
  confirmed$Normal_to_ADM_Associated & confirmed$ADM_to_PDAC_Associated,
  "Both contrasts",
  ifelse(confirmed$Normal_to_ADM_Associated, "Normal-to-ADM",
         ifelse(confirmed$ADM_to_PDAC_Associated, "ADM-to-PDAC",
                "Omnibus only"))
)
confirmed <- confirmed[order(confirmed$Patient_Blocked_F_BH_FDR_all_targets), ]
confirmed$Bulk_Spatial_NormalPDAC_Direction_Concordant <-
  sign(confirmed$PDAC_vs_Normal_log2FC) == sign(confirmed$Normal_to_PDAC_log2)
confirmed$BulkCP_SpatialNormalPDAC_Direction_Concordant <-
  sign(confirmed$PDAC_vs_CP_log2FC) == sign(confirmed$Normal_to_PDAC_log2)
complete <- read.csv("reanalysis/results/spatial/GSE208536_complete_profiles_all_targets.csv")
loo <- read.csv("reanalysis/results/spatial/GSE208536_leave_one_profile_out.csv")
confirmed$Complete_Profile_Omnibus_FDR <-
  complete$Patient_Blocked_F_BH_FDR_all_targets[match(confirmed$Gene, complete$Gene)]
confirmed$LOO_Omnibus_FDR_Max <- vapply(confirmed$Gene, function(gene) {
  max(loo$Patient_Blocked_F_BH_FDR_all_targets[loo$Gene == gene])
}, numeric(1))
confirmed$LOO_Early_Direction_Consistent <- vapply(confirmed$Gene, function(gene) {
  all(sign(loo$Normal_to_ADM_log2[loo$Gene == gene]) ==
        sign(confirmed$Normal_to_ADM_log2[confirmed$Gene == gene]))
}, logical(1))
confirmed$LOO_Late_Direction_Consistent <- vapply(confirmed$Gene, function(gene) {
  all(sign(loo$ADM_to_PDAC_log2[loo$Gene == gene]) ==
        sign(confirmed$ADM_to_PDAC_log2[confirmed$Gene == gene]))
}, logical(1))
write.csv(confirmed[confirmed$CP_Aware, ], file.path(out_dir, "CP_aware_spatial_associated_genes.csv"),
          row.names = FALSE)
write.csv(confirmed, file.path(out_dir, "jointly_measurable_spatial_associated_genes.csv"),
          row.names = FALSE)
write.csv(confirmed[confirmed$CP_Aware_Unadjusted, ],
          file.path(out_dir, "unadjusted_CP_overlap_exploratory.csv"), row.names = FALSE)
write.csv(integrated, file.path(out_dir, "all_jointly_measurable_genes.csv"), row.names = FALSE)
summary <- data.frame(
  Metric = c("GSE143754 annotated genes", "Age-sex-adjusted CP-aware genes",
             "Age-sex-adjusted CP-aware genes on GeoMx panel",
             "Age-sex-adjusted CP-aware spatial overlap",
             "Jointly measurable spatial-associated genes",
             "Spatial overlap: Normal-to-ADM associated", "Spatial overlap: ADM-to-PDAC associated",
             "Spatial overlap: both contrasts", "Unadjusted CP-aware genes",
             "Unadjusted CP-aware spatial overlap"),
  Count = c(nrow(bulk), sum(bulk$PDAC_vs_CP_BH_FDR < 0.05 &
                             abs(bulk$PDAC_vs_CP_log2FC) >= 0.5),
            sum(integrated$CP_Aware), sum(confirmed$CP_Aware), nrow(confirmed),
            sum(confirmed$Normal_to_ADM_Associated),
            sum(confirmed$ADM_to_PDAC_Associated),
            sum(confirmed$Normal_to_ADM_Associated &
                  confirmed$ADM_to_PDAC_Associated),
            sum(bulk$Unadjusted_PDAC_vs_CP_BH_FDR < 0.05 &
                  abs(bulk$Unadjusted_PDAC_vs_CP_log2FC) >= 0.5),
            sum(confirmed$CP_Aware_Unadjusted))
)
write.csv(summary, file.path(out_dir, "analysis_flow_counts.csv"), row.names = FALSE)

top <- head(confirmed$Gene, 20)
map <- read.csv(file.path("reanalysis", "results", "input_audit",
                          "GSE208536_roi_mapping.csv"))
map$Patient <- sprintf("P%02d", match(map$Clinical_Profile, unique(map$Clinical_Profile)))
map$Condition <- factor(map$Condition, levels = c("Normal", "ADM", "PDAC"))
raw <- as.matrix(spatial_xlsx[, -1])
storage.mode(raw) <- "double"
rownames(raw) <- as.character(spatial_xlsx[[1]])
stopifnot(all(top %in% rownames(raw)), identical(colnames(raw), map$Expression_Column))
roi_groups <- split(seq_len(nrow(map)), interaction(map$Patient, map$Condition, drop = TRUE))
heat <- vapply(roi_groups, function(ii) {
  rowMeans(log2(raw[top, ii, drop = FALSE] + 1))
}, numeric(length(top)))
rownames(heat) <- top
labels <- do.call(rbind, lapply(roi_groups, function(ii) {
  data.frame(Patient = map$Patient[ii[[1]]],
             Condition = as.character(map$Condition[ii[[1]]]))
}))
order_index <- order(labels$Patient, match(labels$Condition, c("Normal", "ADM", "PDAC")))
heat <- heat[, order_index, drop = FALSE]
labels <- labels[order_index, ]
column_labels <- paste(labels$Patient, labels$Condition, sep = " ")
standardized <- t(scale(t(heat)))
standardized[!is.finite(standardized)] <- 0
long <- data.frame(
  Gene = rep(rownames(standardized), times = ncol(standardized)),
  Sample = rep(column_labels, each = nrow(standardized)),
  Z = as.vector(standardized)
)
long$Gene <- factor(long$Gene, levels = rev(top))
long$Sample <- factor(long$Sample, levels = column_labels)
heatmap <- ggplot2::ggplot(long, ggplot2::aes(Sample, Gene, fill = Z)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "#237A77", mid = "#F6F7F6",
                                 high = "#B7404B", midpoint = 0,
                                 limits = c(-2.5, 2.5), oob = scales::squish,
                                 name = "Gene z-score") +
  ggplot2::labs(x = "Inferred patient profile and tissue state", y = NULL,
                title = "Spatial tissue-state associations",
                subtitle = "Top 20 jointly measurable genes ranked by spatial omnibus FDR") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(fig_dir, "spatial_tissue_state_heatmap.png"), heatmap,
                width = 10, height = 7, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "spatial_tissue_state_heatmap.pdf"), heatmap,
                width = 10, height = 7)

print(summary, row.names = FALSE)
print(table(confirmed$Stage_Pattern))
print(confirmed[confirmed$Gene %in% c("STAT1", "ITGA2", "TOP2A", "ACTB"),
                c("Gene", "PDAC_vs_CP_log2FC", "PDAC_vs_CP_BH_FDR",
                  "Patient_Blocked_F_BH_FDR_all_targets",
                  "Normal_to_ADM_log2", "ADM_to_PDAC_log2", "Stage_Pattern")],
      row.names = FALSE)
