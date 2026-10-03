#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
fig_dir <- "reanalysis/figures/qc"
out_dir <- "reanalysis/results/qc"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(plot, stem, width = 7, height = 5) {
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".", ext)), plot,
                    width = width, height = height, dpi = 300)
  }
}
train_scores <- read.csv("reanalysis/results/plasma/diagnostic_model_scores.csv")
test_scores <- read.csv("reanalysis/results/plasma/PLCO_locked_model_scores.csv")
performance <- read.csv("reanalysis/results/plasma/diagnostic_to_PLCO_performance.csv")
qc <- read.csv("reanalysis/results/plasma/sample_assay_qc.csv")
qc$Group <- factor(qc$Outcome, levels = c(0, 1), labels = c("Control", "PDAC"))
assay_plot <- ggplot2::ggplot(qc, ggplot2::aes(Cohort, log10(Total_Counts), fill = Group)) +
  ggplot2::geom_boxplot(outlier.size = 0.6) +
  ggplot2::scale_fill_manual(values = c(Control = "#237A77", PDAC = "#B7404B")) +
  ggplot2::labs(x = NULL, y = "Log10 deposited total counts", fill = NULL,
                title = "Plasma assay count distributions") +
  ggplot2::theme_classic(base_size = 11)
save_plot(assay_plot, "plasma_total_counts")
calibration_bins <- do.call(rbind, lapply(performance$Model, function(model) {
  score <- test_scores[[model]]
  bins <- cut(rank(score, ties.method = "first"), breaks = seq(0, 96, length.out = 7))
  do.call(rbind, lapply(split(seq_along(score), bins), function(ii) {
    interval <- binom.test(sum(test_scores$Outcome[ii]), length(ii))$conf.int
    data.frame(Model = gsub("_", " ", model), N = length(ii),
               Mean_Score = mean(score[ii]), Observed_Fraction = mean(test_scores$Outcome[ii]),
               Observed_CI_low = interval[[1]], Observed_CI_high = interval[[2]])
  }))
}))
write.csv(calibration_bins, file.path(out_dir, "PLCO_calibration_bins.csv"), row.names = FALSE)
calibration_plot <- ggplot2::ggplot(calibration_bins,
                                    ggplot2::aes(Mean_Score, Observed_Fraction)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#7A8490") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Observed_CI_low, ymax = Observed_CI_high),
                         width = 0.015, color = "#237A77") +
  ggplot2::geom_point(size = 2, color = "#B7404B") +
  ggplot2::facet_wrap(~Model, ncol = 2) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "Mean model score in rank bin", y = "Observed case fraction",
                title = "PLCO calibration in the deposited case-control sample",
                subtitle = "Six equal-size bins; intervals are exact binomial 95% CIs") +
  ggplot2::theme_classic(base_size = 10)
save_plot(calibration_plot, "PLCO_calibration", 7, 6)
clinical <- read.csv("reanalysis/results/bulk/GSE143754_sample_metadata.csv")
clinical$Age <- as.numeric(sub("^age: ", "", clinical$Sample_characteristics_ch1))
clinical$Sex <- sub("^Sex: ", "", clinical$Sample_characteristics_ch1.1)
age_summary <- aggregate(Age ~ Group, clinical,
                          function(x) c(N = length(x), Median = median(x), Min = min(x), Max = max(x)))
write.csv(age_summary, file.path(out_dir, "bulk_age_summary.csv"), row.names = FALSE)
age_plot <- ggplot2::ggplot(clinical, ggplot2::aes(Group, Age, color = Sex)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 203),
                      size = 2.4) +
  ggplot2::scale_color_manual(values = c(Female = "#873D62", Male = "#237A77")) +
  ggplot2::labs(x = NULL, y = "Age (years)", title = "Bulk tissue cohort demographics") +
  ggplot2::theme_classic(base_size = 11)
save_plot(age_plot, "bulk_demographics")
cat("QC and calibration figures generated.\n")
