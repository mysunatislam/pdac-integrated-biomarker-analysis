#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source("scripts/reanalysis/00_audit_inputs.R", local = TRUE)
out_dir <- file.path("reanalysis", "results", "plasma")
fig_dir <- file.path("reanalysis", "figures", "plasma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

plot_only <- "--figures-only" %in% commandArgs(trailingOnly = TRUE)
if (!plot_only) {
read_assay <- function(path) {
  b <- suppressMessages(readxl::read_excel(path, sheet = "Raw", .name_repair = "minimal"))
  ids <- as.character(unlist(b[1, -1], use.names = FALSE))
  totals <- as.numeric(unlist(b[2, -1], use.names = FALSE))
  features <- as.character(b[[1]][-(1:2)])
  raw <- as.matrix(b[-c(1, 2), -1])
  storage.mode(raw) <- "double"
  rownames(raw) <- features
  colnames(raw) <- ids
  stopifnot(all(is.finite(raw)), all(raw >= 0), all(is.finite(totals)),
            all(totals > 0), !anyDuplicated(ids), !anyDuplicated(features))
  x <- t(log2(sweep(raw, 2, totals, "/") * 1e6 + 1))
  list(X = x, ids = ids, totals = totals)
}

train <- read_assay(file.path(source_dir, "GSE259327_Raw_counts_PDAC_203samples_new.xlsx"))
test <- read_assay(file.path(source_dir, "GSE259327_Raw_counts_PLCO_96samples_new.xlsx"))
stopifnot(identical(colnames(train$X), colnames(test$X)))
metadata <- plasma_meta[match(c(train$ids, test$ids), plasma_meta$Title), ]
stopifnot(identical(metadata$Title, c(train$ids, test$ids)))
label <- metadata[["Sample_characteristics_ch1.2"]]
y_train <- as.integer(label[seq_along(train$ids)] == "disease: PDAC")
y_test <- as.integer(label[length(train$ids) + seq_along(test$ids)] ==
                       "disease: PLCO_Pre-PDAC")
stopifnot(sum(y_train) == 121, sum(y_test) == 48)

historical6 <- c("miR-184", "miR-107", "miR-216b-5p", "miR-20a-5p",
                 "miR-134-3p", "miR-215-3p")
published3 <- c("let-7i-5p", "miR-130a-3p", "miR-221-3p")
candidate20 <- read.csv(file.path("reanalysis", "results", "ev",
                                  "GSE304572_EV_candidate20.csv"))$miRNA
stopifnot(length(candidate20) == 20,
          all(c(historical6, published3, candidate20) %in% colnames(train$X)))

fit_panel <- function(x, y, features) {
  d <- as.data.frame(x[, features, drop = FALSE], check.names = FALSE)
  names(d) <- paste0("x", seq_along(features))
  d$Outcome <- y
  fit <- suppressWarnings(glm(Outcome ~ ., data = d, family = binomial(),
                              control = glm.control(maxit = 100)))
  stopifnot(fit$converged, all(is.finite(coef(fit))))
  fit
}

predict_panel <- function(fit, x, features) {
  d <- as.data.frame(x[, features, drop = FALSE], check.names = FALSE)
  names(d) <- paste0("x", seq_along(features))
  as.numeric(predict(fit, newdata = d, type = "response"))
}

stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (class in 0:1) {
    ids <- sample(which(y == class))
    fold[ids] <- rep(seq_len(k), length.out = length(ids))
  }
  fold
}

select_features <- function(x, y, candidates, n, filter_detection = FALSE) {
  if (filter_detection) {
    detected <- colSums(x[, candidates, drop = FALSE] >= log2(2)) >=
      max(5, ceiling(0.1 * nrow(x)))
    candidates <- candidates[detected]
  }
  p <- vapply(candidates, function(feature) {
    if (length(unique(x[, feature])) < 2) return(1)
    suppressWarnings(wilcox.test(x[y == 1, feature], x[y == 0, feature],
                                 exact = FALSE)$p.value)
  }, numeric(1))
  stopifnot(length(candidates) >= n, all(is.finite(p)))
  head(candidates[order(p, candidates)], n)
}

auc_value <- function(y, p) {
  as.numeric(pROC::auc(pROC::roc(y, p, levels = c(0, 1),
                                direction = "<", quiet = TRUE)))
}

crossval_fixed <- function(x, y, features, seed = 20261001) {
  fold <- stratified_folds(y, 10, seed)
  pred <- rep(NA_real_, length(y))
  for (i in seq_len(10)) {
    fit <- fit_panel(x[fold != i, , drop = FALSE], y[fold != i], features)
    pred[fold == i] <- predict_panel(fit, x[fold == i, , drop = FALSE], features)
  }
  stopifnot(all(is.finite(pred)))
  pred
}

inner_cv_auc <- function(x, y, candidates, n, seed, filter_detection = FALSE) {
  fold <- stratified_folds(y, 5, seed)
  pred <- rep(NA_real_, length(y))
  for (i in seq_len(5)) {
    features <- select_features(x[fold != i, , drop = FALSE],
                                y[fold != i], candidates, n, filter_detection)
    fit <- fit_panel(x[fold != i, , drop = FALSE], y[fold != i], features)
    pred[fold == i] <- predict_panel(fit, x[fold == i, , drop = FALSE], features)
  }
  stopifnot(all(is.finite(pred)))
  auc_value(y, pred)
}

candidate_sizes <- c(3L, 6L, 9L)
nested_selection <- function(candidates, prefix, filter_detection = FALSE) {
  outer_fold <- stratified_folds(y_train, 10, 20261002)
  oof <- rep(NA_real_, length(y_train))
  outer_selection <- vector("list", 10)
  for (i in seq_len(10)) {
    x_fit <- train$X[outer_fold != i, , drop = FALSE]
    y_fit <- y_train[outer_fold != i]
    inner_auc <- vapply(candidate_sizes, function(n) {
      inner_cv_auc(x_fit, y_fit, candidates, n, 20261002 + i, filter_detection)
    }, numeric(1))
    n <- candidate_sizes[which.max(inner_auc)]
    features <- select_features(x_fit, y_fit, candidates, n, filter_detection)
    model <- fit_panel(x_fit, y_fit, features)
    oof[outer_fold == i] <- predict_panel(
      model, train$X[outer_fold == i, , drop = FALSE], features
    )
    outer_selection[[i]] <- data.frame(Outer_Fold = i, Feature_Count = n,
                                       Inner_AUC = max(inner_auc),
                                       Features = paste(features, collapse = ";"))
    cat(prefix, "outer fold", i, "complete\n")
  }
  stopifnot(all(is.finite(oof)))
  write.csv(do.call(rbind, outer_selection),
            file.path(out_dir, paste0(prefix, "_pipeline_outer_folds.csv")), row.names = FALSE)
  full_inner_auc <- vapply(candidate_sizes, function(n) {
    inner_cv_auc(train$X, y_train, candidates, n, 20261025, filter_detection)
  }, numeric(1))
  full_n <- candidate_sizes[which.max(full_inner_auc)]
  features <- select_features(train$X, y_train, candidates, full_n, filter_detection)
  write.csv(data.frame(Feature_Count = candidate_sizes, Inner_CV_AUC = full_inner_auc),
            file.path(out_dir, paste0(prefix, "_pipeline_feature_count_selection.csv")),
            row.names = FALSE)
  list(features = features, oof = oof)
}
ev_pipeline <- nested_selection(candidate20, "EV")
all_mirnas <- colnames(train$X)[grepl("^(miR-|let-)", colnames(train$X))]
all_pipeline <- nested_selection(all_mirnas, "Unrestricted", filter_detection = TRUE)

models <- list(
  Historical_six = list(features = historical6,
                        oof = crossval_fixed(train$X, y_train, historical6)),
  Published_three_refit = list(features = published3,
                               oof = crossval_fixed(train$X, y_train, published3)),
  EV_ranked_nested = ev_pipeline,
  Unrestricted_nested = all_pipeline
)
for (name in names(models)) {
  fit <- fit_panel(train$X, y_train, models[[name]]$features)
  models[[name]]$fit <- fit
  models[[name]]$train_score <- predict_panel(fit, train$X, models[[name]]$features)
  models[[name]]$plco_score <- predict_panel(fit, test$X, models[[name]]$features)
}

clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
score_summary <- function(name, item) {
  roc_train <- pROC::roc(y_train, item$train_score, levels = c(0, 1),
                         direction = "<", quiet = TRUE)
  roc_test <- pROC::roc(y_test, item$plco_score, levels = c(0, 1),
                        direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_test, method = "delong"))
  cv_roc <- pROC::roc(y_train, item$oof, levels = c(0, 1), direction = "<", quiet = TRUE)
  cv_ci <- as.numeric(pROC::ci.auc(cv_roc, method = "delong"))
  threshold <- as.numeric(pROC::coords(roc_train, x = "best",
                                       best.method = "youden",
                                       ret = "threshold"))[[1]]
  positive <- item$plco_score >= threshold
  sens_ci <- binom.test(sum(positive[y_test == 1]), sum(y_test == 1))$conf.int
  spec_ci <- binom.test(sum(!positive[y_test == 0]), sum(y_test == 0))$conf.int
  calibration <- glm(y_test ~ qlogis(clip(item$plco_score)),
                     family = binomial())
  data.frame(
    Model = name, Features = paste(item$features, collapse = ";"),
    Apparent_AUC = as.numeric(pROC::auc(roc_train)),
    Diagnostic_CV_AUC = auc_value(y_train, item$oof),
    Diagnostic_CV_Conditional_CI_low = cv_ci[[1]],
    Diagnostic_CV_Conditional_CI_high = cv_ci[[3]],
    PLCO_AUC = as.numeric(pROC::auc(roc_test)),
    PLCO_AUC_CI_low = ci[[1]], PLCO_AUC_CI_high = ci[[3]],
    Training_Youden_Threshold = threshold,
    PLCO_Sensitivity = mean(positive[y_test == 1]),
    PLCO_Sensitivity_CI_low = sens_ci[[1]], PLCO_Sensitivity_CI_high = sens_ci[[2]],
    PLCO_Specificity = mean(!positive[y_test == 0]),
    PLCO_Specificity_CI_low = spec_ci[[1]], PLCO_Specificity_CI_high = spec_ci[[2]],
    PLCO_Brier = mean((item$plco_score - y_test)^2),
    PLCO_Calibration_Intercept = unname(coef(calibration)[[1]]),
    PLCO_Calibration_Slope = unname(coef(calibration)[[2]])
  )
}
summary <- do.call(rbind, lapply(names(models), function(name) {
  score_summary(name, models[[name]])
}))
write.csv(summary, file.path(out_dir, "diagnostic_to_PLCO_performance.csv"),
          row.names = FALSE)

scores <- data.frame(Sample = test$ids, Outcome = y_test)
for (name in names(models)) scores[[name]] <- models[[name]]$plco_score
write.csv(scores, file.path(out_dir, "PLCO_locked_model_scores.csv"), row.names = FALSE)

diagnostic_scores <- data.frame(Sample = train$ids, Outcome = y_train)
for (name in names(models)) {
  diagnostic_scores[[paste0(name, "_Fitted")]] <- models[[name]]$train_score
  diagnostic_scores[[paste0(name, "_OOF")]] <- models[[name]]$oof
}
write.csv(diagnostic_scores, file.path(out_dir, "diagnostic_model_scores.csv"), row.names = FALSE)
coefficients <- do.call(rbind, lapply(names(models), function(name) {
  data.frame(Model = name, Feature = c("(Intercept)", models[[name]]$features),
             Estimate = as.numeric(coef(models[[name]]$fit)))
}))
write.csv(coefficients, file.path(out_dir, "model_coefficients.csv"), row.names = FALSE)
qc <- data.frame(Sample = c(train$ids, test$ids),
                  Cohort = c(rep("Diagnostic", length(train$ids)), rep("PLCO", length(test$ids))),
                  Outcome = c(y_train, y_test), Total_Counts = c(train$totals, test$totals))
write.csv(qc, file.path(out_dir, "sample_assay_qc.csv"), row.names = FALSE)
test_rocs <- lapply(models, function(item) {
  pROC::roc(y_test, item$plco_score, levels = c(0, 1), direction = "<", quiet = TRUE)
})
pairs <- combn(names(models), 2, simplify = FALSE)
comparisons <- do.call(rbind, lapply(pairs, function(pair) {
  comparison <- pROC::roc.test(test_rocs[[pair[[1]]]], test_rocs[[pair[[2]]]],
                               method = "delong", paired = TRUE)
  data.frame(Model_1 = pair[[1]], Model_2 = pair[[2]],
             PLCO_Delta_AUC = as.numeric(pROC::auc(test_rocs[[pair[[1]]]]) -
                                          pROC::auc(test_rocs[[pair[[2]]]])),
             Paired_DeLong_p = comparison$p.value)
}))
comparisons$BH_FDR <- p.adjust(comparisons$Paired_DeLong_p, "BH")
write.csv(comparisons,
          file.path(out_dir, "PLCO_paired_model_comparison.csv"), row.names = FALSE)
} else {
  summary <- read.csv(file.path(out_dir, "diagnostic_to_PLCO_performance.csv"))
  scores <- read.csv(file.path(out_dir, "PLCO_locked_model_scores.csv"))
  diagnostic_scores <- read.csv(file.path(out_dir, "diagnostic_model_scores.csv"))
  qc <- read.csv(file.path(out_dir, "sample_assay_qc.csv"))
  y_test <- scores$Outcome
  models <- setNames(lapply(summary$Model, function(name) {
    list(train_score = diagnostic_scores[[paste0(name, "_Fitted")]],
         plco_score = scores[[name]])
  }), summary$Model)
  test_rocs <- lapply(models, function(item) {
    pROC::roc(y_test, item$plco_score, levels = c(0, 1), direction = "<", quiet = TRUE)
  })
}

display_names <- c(Historical_six = "Historical six", Published_three_refit = "Published three (refit)",
                    EV_ranked_nested = "EV-ranked selection", Unrestricted_nested = "Unrestricted selection")
colors <- c("#B7404B", "#237A77", "#364B72", "#A16B20")
roc_data <- do.call(rbind, lapply(names(test_rocs), function(name) {
  r <- test_rocs[[name]]
  data.frame(Model = name, FPR = 1 - r$specificities, TPR = r$sensitivities)
}))
roc_data$Model <- factor(roc_data$Model, levels = names(models))
legend_labels <- paste0(display_names[names(models)], ": AUC ",
                         sprintf("%.3f", summary$PLCO_AUC))
roc_plot <- ggplot2::ggplot(roc_data, ggplot2::aes(FPR, TPR, color = Model)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "#A0A0A0", linetype = 2) +
  ggplot2::geom_path(linewidth = 0.8) +
  ggplot2::scale_color_manual(values = colors, labels = legend_labels) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = "False-positive rate", y = "True-positive rate", color = NULL,
                title = "Prediagnostic PLCO evaluation",
                subtitle = "48 cases and 48 controls; models fitted on diagnostic data") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
save_plot <- function(plot, stem, width = 6, height = 6) {
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".png")), plot,
                  width = width, height = height, dpi = 300)
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height)
}
save_plot(roc_plot, "PLCO_ROC_models", 6, 7)
score_long <- do.call(rbind, lapply(names(models), function(name) {
  data.frame(Model = display_names[[name]], Cohort = qc$Cohort,
             Outcome = factor(qc$Outcome, levels = c(0, 1), labels = c("Control", "PDAC")),
             Score = c(models[[name]]$train_score, models[[name]]$plco_score))
}))
score_distribution <- ggplot2::ggplot(score_long,
                                       ggplot2::aes(Cohort, log10(pmax(Score, 1e-10)), fill = Outcome)) +
  ggplot2::geom_boxplot(outlier.size = 0.7, linewidth = 0.4) +
  ggplot2::facet_wrap(~Model, ncol = 2) +
  ggplot2::scale_fill_manual(values = c(Control = "#237A77", PDAC = "#B7404B")) +
  ggplot2::labs(x = NULL, y = "Log10 model probability", fill = NULL,
                title = "Model score shift across cohorts") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
save_plot(score_distribution, "score_distributions", 7, 5.5)
score_ranges <- aggregate(Score ~ Model + Cohort + Outcome, score_long,
                           function(x) c(Min = min(x), Median = median(x), Max = max(x)))
write.csv(score_ranges, file.path(out_dir, "score_distribution_summary.csv"), row.names = FALSE)

if (!plot_only) {
  cat("Final EV-ranked features:", paste(ev_pipeline$features, collapse = ", "), "\n")
  cat("Final unrestricted features:", paste(all_pipeline$features, collapse = ", "), "\n")
  print(comparisons, row.names = FALSE)
}
print(summary, row.names = FALSE)
