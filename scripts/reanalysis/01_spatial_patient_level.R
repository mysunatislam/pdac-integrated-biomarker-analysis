#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source("scripts/reanalysis/00_audit_inputs.R", local = TRUE)

out_dir <- file.path("reanalysis", "results", "spatial")
fig_dir <- file.path("reanalysis", "figures", "spatial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

map <- spatial
map$Patient <- sprintf("P%02d", match(map$Clinical_Profile, unique(map$Clinical_Profile)))
map$Condition <- factor(map$Condition, levels = c("Normal", "ADM", "PDAC"))
stopifnot(length(unique(map$Patient)) == 8, all(table(map$Condition) == 16))

target <- as.character(spatial_xlsx[[1]])
raw_expr <- as.matrix(spatial_xlsx[, -1])
storage.mode(raw_expr) <- "double"
stopifnot(!anyDuplicated(target), all(is.finite(raw_expr)), all(raw_expr >= 0),
          identical(colnames(raw_expr), map$Expression_Column))
rownames(raw_expr) <- target
expr <- log2(raw_expr + 1)

roi_groups <- split(seq_len(nrow(map)), interaction(map$Patient, map$Condition, drop = TRUE))
patient_means <- vapply(roi_groups, function(ii) rowMeans(expr[, ii, drop = FALSE]),
                        numeric(nrow(expr)))
rownames(patient_means) <- target
stage_meta <- do.call(rbind, lapply(roi_groups, function(ii) {
  data.frame(Patient = map$Patient[ii[[1]]], Condition = as.character(map$Condition[ii[[1]]]),
             ROI_Count = length(ii))
}))
stage_meta$Patient <- factor(stage_meta$Patient)
stage_meta$Condition <- factor(stage_meta$Condition, levels = c("Normal", "ADM", "PDAC"))
stopifnot(ncol(patient_means) == nrow(stage_meta))
write.csv(stage_meta, file.path(out_dir, "patient_stage_roi_counts.csv"), row.names = FALSE)

contrast_stats <- function(fit, contrast) {
  beta <- coef(fit)
  covariance <- vcov(fit)
  l <- setNames(rep(0, length(beta)), names(beta))
  l[names(contrast)] <- contrast
  estimate <- sum(l * beta)
  se <- sqrt(as.numeric(t(l) %*% covariance %*% l))
  p <- 2 * pt(-abs(estimate / se), df.residual(fit))
  c(estimate = estimate, SE = se, p = p,
    CI_low = estimate - qt(0.975, df.residual(fit)) * se,
    CI_high = estimate + qt(0.975, df.residual(fit)) * se)
}

patient_test <- function(gene, ii = seq_len(nrow(stage_meta))) {
  d <- droplevels(stage_meta[ii, ])
  d$Expression <- as.numeric(patient_means[gene, ii])
  fit <- lm(Expression ~ Patient + Condition, data = d)
  an <- anova(fit)
  early <- contrast_stats(fit, c(ConditionADM = 1))
  late <- contrast_stats(fit, c(ConditionADM = -1, ConditionPDAC = 1))
  overall <- contrast_stats(fit, c(ConditionPDAC = 1))
  data.frame(
    Gene = gene, Patient_Blocked_F_p = an["Condition", "Pr(>F)"],
    Normal_to_ADM_log2 = early["estimate"], Normal_to_ADM_p = early["p"],
    Normal_to_ADM_CI_low = early["CI_low"], Normal_to_ADM_CI_high = early["CI_high"],
    ADM_to_PDAC_log2 = late["estimate"], ADM_to_PDAC_p = late["p"],
    ADM_to_PDAC_CI_low = late["CI_low"], ADM_to_PDAC_CI_high = late["CI_high"],
    Normal_to_PDAC_log2 = overall["estimate"], Normal_to_PDAC_p = overall["p"]
  )
}

stats <- do.call(rbind, lapply(target, patient_test))
for (column in c("Patient_Blocked_F_p", "Normal_to_ADM_p", "ADM_to_PDAC_p",
                 "Normal_to_PDAC_p")) {
  stats[[sub("_p$", "_BH_FDR_all_targets", column)]] <- p.adjust(stats[[column]], method = "BH")
}
stats$ADM_peak <- with(stats, Normal_to_ADM_BH_FDR_all_targets < 0.05 &
                         ADM_to_PDAC_BH_FDR_all_targets < 0.05 &
                         Normal_to_ADM_log2 > 0 & ADM_to_PDAC_log2 < 0)
stats$ADM_trough <- with(stats, Normal_to_ADM_BH_FDR_all_targets < 0.05 &
                           ADM_to_PDAC_BH_FDR_all_targets < 0.05 &
                           Normal_to_ADM_log2 < 0 & ADM_to_PDAC_log2 > 0)
stats <- stats[order(stats$Patient_Blocked_F_p), ]
write.csv(stats, file.path(out_dir, "GSE208536_patient_blocked_all_targets.csv"), row.names = FALSE)

complete_patients <- names(which(table(stage_meta$Patient) == 3))
complete_index <- which(stage_meta$Patient %in% complete_patients)
complete_stats <- do.call(rbind, lapply(target, patient_test, ii = complete_index))
complete_stats$Patient_Blocked_F_BH_FDR_all_targets <-
  p.adjust(complete_stats$Patient_Blocked_F_p, "BH")
write.csv(complete_stats, file.path(out_dir, "GSE208536_complete_profiles_all_targets.csv"),
          row.names = FALSE)
loo_stats <- do.call(rbind, lapply(levels(stage_meta$Patient), function(patient) {
  ii <- which(stage_meta$Patient != patient)
  result <- do.call(rbind, lapply(target, patient_test, ii = ii))
  result$Left_Out_Profile <- patient
  result$Patient_Blocked_F_BH_FDR_all_targets <-
    p.adjust(result$Patient_Blocked_F_p, "BH")
  result
}))
write.csv(loo_stats, file.path(out_dir, "GSE208536_leave_one_profile_out.csv"),
          row.names = FALSE)

historical80 <- scan("results/tables/spatial/valid_genes.txt", what = character(), quiet = TRUE)
stopifnot(length(historical80) == 80, all(historical80 %in% target))
stats80 <- stats[match(historical80, stats$Gene), ]

mixed_test <- function(gene) {
  d <- map[, c("Patient", "Condition")]
  d$Expression <- as.numeric(expr[gene, ])
  fit <- tryCatch(
    nlme::lme(Expression ~ Condition, random = ~1 | Patient, data = d,
              method = "REML", control = nlme::lmeControl(msMaxIter = 100)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  as.numeric(anova(fit)["Condition", "p-value"])
}
stats80$ROI_Mixed_Model_p <- vapply(stats80$Gene, mixed_test, numeric(1))
stats80$ROI_Mixed_Model_BH_FDR_80 <- p.adjust(stats80$ROI_Mixed_Model_p, method = "BH")
write.csv(stats80, file.path(out_dir, "GSE208536_historical80_patient_aware.csv"), row.names = FALSE)

if ("STAT1" %in% target) {
  stat1 <- stage_meta
  stat1$Expression <- as.numeric(patient_means["STAT1", ])
  write.csv(stat1, file.path(out_dir, "STAT1_patient_stage_means.csv"), row.names = FALSE)
  plot_stat1 <- ggplot2::ggplot(stat1, ggplot2::aes(Condition, Expression, group = Patient)) +
    ggplot2::geom_line(color = "#7A8490", linewidth = 0.55, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(color = Condition), size = 2.6) +
    ggplot2::scale_color_manual(values = c(Normal = "#237A77", ADM = "#C86B31",
                                            PDAC = "#873D62"), guide = "none") +
    ggplot2::labs(x = NULL, y = "Mean log2(normalized expression + 1)",
                  title = "STAT1 across tissue states",
                  subtitle = "ROI means within eight inferred patient profiles") +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(file.path(fig_dir, "STAT1_patient_stage.png"), plot_stat1,
                  width = 5.5, height = 4, dpi = 300)
  ggplot2::ggsave(file.path(fig_dir, "STAT1_patient_stage.pdf"), plot_stat1,
                  width = 5.5, height = 4)
}

cat("Measured targets:", nrow(stats), "\n")
cat("All-target FDR < 0.05:", sum(stats$Patient_Blocked_F_BH_FDR_all_targets < 0.05), "\n")
cat("Historical 80 at all-target FDR < 0.05:",
    sum(stats80$Patient_Blocked_F_BH_FDR_all_targets < 0.05), "\n")
cat("Historical 80 at mixed-model FDR < 0.05:",
    sum(stats80$ROI_Mixed_Model_BH_FDR_80 < 0.05, na.rm = TRUE), "\n")
cat("Mixed-model fit failures:", sum(is.na(stats80$ROI_Mixed_Model_p)), "\n")
print(stats[stats$Gene == "STAT1", ], row.names = FALSE)
