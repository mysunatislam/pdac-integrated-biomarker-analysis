#!/usr/bin/env Rscript
source("scripts/reanalysis/external_helpers.R")
historical6 <- c("miR-184", "miR-107", "miR-216b-5p", "miR-20a-5p", "miR-134-3p", "miR-215-3p")
published3 <- c("let-7i-5p", "miR-130a-3p", "miR-221-3p")
candidate20 <- read.csv("reanalysis/results/ev/GSE304572_EV_candidate20.csv")$miRNA
markers <- sort(unique(c(historical6, published3, candidate20)))

# Lock marker orientation to the previously analyzed diagnostic plasma cohort.
b <- readxl::read_excel("data/source/GSE259327_Raw_counts_PDAC_203samples_new.xlsx",
                        sheet = "Raw", .name_repair = "minimal")
ids <- as.character(unlist(b[1, -1], use.names = FALSE))
totals <- as.numeric(unlist(b[2, -1], use.names = FALSE))
features <- as.character(b[[1]][-(1:2)])
counts <- as.matrix(b[-c(1, 2), -1])
storage.mode(counts) <- "double"
expression <- log2(sweep(counts, 2, totals, "/") * 1e6 + 1)
meta <- read.csv("reanalysis/results/input_audit/GSE259327_sample_metadata.csv")
meta <- meta[match(ids, meta$Title), ]
case <- meta$Sample_characteristics_ch1.2 == "disease: PDAC"
stopifnot(sum(case) == 121, identical(meta$Title, ids), all(markers %in% features))
difference <- rowMeans(expression[, case]) - rowMeans(expression[, !case])
names(difference) <- features
marker_reference <- data.frame(miRNA = markers, Diagnostic_Plasma_Mean_Difference = difference[markers],
                               Fixed_Direction = sign(difference[markers]),
                               Historical6 = markers %in% historical6,
                               Published3 = markers %in% published3,
                               EV_Candidate20 = markers %in% candidate20)
stopifnot(all(marker_reference$Fixed_Direction != 0))
write.csv(marker_reference, file.path(external_dir, "frozen_miRNA_reference.csv"), row.names = FALSE)

annotation <- read.csv("data/external/GPL18941_annotation.csv")
stopifnot(!anyDuplicated(annotation$ID))
all_results <- list()
roc_results <- list()
coverage <- list()
serum <- list()
roc_stats <- function(y, score) {
  r <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- pROC::ci.auc(r, method = "delong")
  c(AUC = as.numeric(pROC::auc(r)), CI_low = as.numeric(ci[1]), CI_high = as.numeric(ci[3]))
}
for (id in c("GSE59856", "GSE85589")) {
  cat(id, "age/sex-adjusted serum comparisons\n")
  z <- read_external_geo(id)
  mature_id <- sub("_st$", "", rownames(z$x))
  symbol <- sub("^hsa-", "", annotation$miRNA_ID_LIST[match(mature_id, annotation$ID)])
  unambiguous <- grepl("^MIMAT[0-9]+$", mature_id) & !is.na(symbol) &
    !grepl("[/,;[:space:]]", symbol)
  finite_rows <- rowSums(is.na(z$x)) == 0
  keep <- unambiguous & finite_rows
  gene <- collapse_probes(z$x[keep, , drop = FALSE], symbol[keep])
  x <- gene$x
  coverage[[id]] <- data.frame(Dataset = id, miRNA = markers,
                              Measurable_Complete = markers %in% rownames(x))
  write.csv(data.frame(Probe = rownames(z$x), Mature_ID = mature_id, miRNA = symbol,
                      Unambiguous = unambiguous, Complete_All_Samples = finite_rows),
            file.path(external_dir, paste0(id, "_miRNA_annotation_audit.csv")), row.names = FALSE)
  meta <- z$meta
  meta$Age <- suppressWarnings(as.numeric(characteristic(meta, "age:")))
  meta$Sex <- factor(characteristic(meta, "gender:"))
  if (id == "GSE59856") {
    meta$Disease <- characteristic(meta, "disease state:")
    meta$Stage <- characteristic(meta, "tumor stage:")
    meta$Comparator <- ifelse(meta$Disease == "pancreatic cancer", "PC",
                         ifelse(meta$Disease == "healthy control", "Healthy",
                           ifelse(meta$Disease == "benign pancreatic or biliary tract diseases",
                                  "Benign", "Other_cancer")))
    stopifnot(identical(as.integer(table(factor(meta$Comparator,
                      levels = c("PC", "Healthy", "Benign", "Other_cancer")))), c(100L, 150L, 21L, 300L)))
  } else {
    source <- trimws(meta$Sample_source_name_ch1)
    meta$Disease <- sub("^serum miRNA \\((.*)\\).*$", "\\1", source)
    meta$CA19_9 <- characteristic(meta, "ca19-9:")
    meta$Comparator <- ifelse(meta$Disease == "pancreatic cancer", "PC",
                         ifelse(meta$Disease == "healthy control", "Healthy",
                           ifelse(meta$Disease == "cholelithiasis", "Benign", "Other_cancer")))
    stopifnot(identical(as.integer(table(factor(meta$Comparator,
                      levels = c("PC", "Healthy", "Benign", "Other_cancer")))), c(88L, 19L, 10L, 115L)))
  }
  stopifnot(nlevels(meta$Sex) == 2, all(is.finite(x)))
  meta$Age_Available <- is.finite(meta$Age)
  write.csv(meta, file.path(external_dir, paste0(id, "_serum_analysis_mapping.csv")), row.names = FALSE)
  serum[[id]] <- list(x = x, meta = meta)
  for (control in c("Healthy", "Benign", "Other_cancer")) {
    available <- meta$Comparator %in% c("PC", control)
    keep <- available & meta$Age_Available
    d <- meta[keep, ]
    d$State <- factor(d$Comparator, levels = c(control, "PC"))
    design <- model.matrix(~State + Age + Sex, d)
    stopifnot(qr(design)$rank == ncol(design))
    fit <- limma::eBayes(limma::lmFit(x[, keep, drop = FALSE], design), trend = TRUE, robust = TRUE)
    stats <- limma_result(fit, "StatePC", id, paste0("PC_vs_", control))
    stats$Adjusted_Cases <- sum(d$State == "PC")
    stats$Adjusted_Controls <- sum(d$State == control)
    stats$Missing_Age_Excluded <- sum(available & !meta$Age_Available)
    names(stats)[names(stats) == "Gene"] <- "miRNA"
    all_results[[length(all_results) + 1]] <- stats
    keep <- available
    d <- meta[keep, ]
    for (marker in intersect(markers, rownames(x))) {
      signed <- x[marker, keep] * marker_reference$Fixed_Direction[match(marker, markers)]
      r <- roc_stats(as.integer(d$Comparator == "PC"), signed)
      roc_results[[length(roc_results) + 1]] <- data.frame(
        Dataset = id, Comparison = paste0("PC_vs_", control), Subgroup = "All",
        miRNA = marker, Cases = sum(d$Comparator == "PC"), Controls = sum(d$Comparator == control),
        AUC = r["AUC"], CI_low = r["CI_low"], CI_high = r["CI_high"])
    }
  }
  # Predefined clinical subgroups remain diagnostic case-control comparisons.
  subgroup <- if (id == "GSE59856") meta$Comparator == "PC" & meta$Stage %in% c("pStage IIA", "pStage IIB") else
    meta$Comparator == "PC" & meta$CA19_9 == "<37"
  for (control in c("Healthy", "Benign")) {
    eligible_control <- meta$Comparator == control
    if (id == "GSE85589") eligible_control <- eligible_control & meta$CA19_9 == "<37"
    keep <- subgroup | eligible_control
    d <- meta[keep, ]
    for (marker in intersect(markers, rownames(x))) {
      signed <- x[marker, keep] * marker_reference$Fixed_Direction[match(marker, markers)]
      r <- roc_stats(as.integer(d$Comparator == "PC"), signed)
      roc_results[[length(roc_results) + 1]] <- data.frame(
        Dataset = id, Comparison = paste0("PC_vs_", control),
        Subgroup = ifelse(id == "GSE59856", "Pathologic_stage_II_non_yp", "CA19_9_less_than_37"),
        miRNA = marker, Cases = sum(d$Comparator == "PC"), Controls = sum(d$Comparator == control),
        AUC = r["AUC"], CI_low = r["CI_low"], CI_high = r["CI_high"])
    }
  }
}
stats <- do.call(rbind, all_results)
family <- merge(stats[stats$miRNA %in% markers, ], marker_reference, by = "miRNA")
family$BH_frozen_within_contrast <- ave(family$P, interaction(family$Dataset, family$Comparison),
                                     FUN = function(p) p.adjust(p, "BH"))
family$BH_frozen_all_comparisons <- ave(family$P, family$Dataset, FUN = function(p) p.adjust(p, "BH"))
family$Same_Direction_As_Plasma <- sign(family$log2FC) == family$Fixed_Direction
write.csv(stats, file.path(external_dir, "serum_all_miRNA_contrasts.csv"), row.names = FALSE)
write.csv(family, file.path(external_dir, "serum_frozen_miRNA_family.csv"), row.names = FALSE)
write.csv(do.call(rbind, roc_results), file.path(external_dir, "serum_fixed_direction_marker_AUC.csv"), row.names = FALSE)
write.csv(do.call(rbind, coverage), file.path(external_dir, "serum_marker_coverage.csv"), row.names = FALSE)

cat("GSE85589: exploratory nested-CV CA19-9 increment analysis\n")
z <- serum[["GSE85589"]]
keep <- z$meta$Comparator %in% c("PC", "Healthy", "Benign") & z$meta$CA19_9 %in% c("<37", ">37")
meta <- z$meta[keep, ]
y <- as.integer(meta$Comparator == "PC")
base <- cbind(Age = meta$Age, Male = as.integer(meta$Sex == "M"), CA19_9_positive = as.integer(meta$CA19_9 == ">37"))
stopifnot(sum(y) == 88, sum(!y) == 27, all(c(historical6, published3) %in% rownames(z$x)))
stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  out <- integer(length(y))
  for (class in 0:1) {
    index <- sample(which(y == class))
    out[index] <- rep(seq_len(k), length.out = length(index))
  }
  out
}
fold <- stratified_folds(y, 10, 20261004)
panels <- list(CA19_9_age_sex = character(), CA19_9_plus_published3 = published3,
               CA19_9_plus_historical6 = historical6)
predictions <- data.frame(GSM = meta$GSM, Case = y, Comparator = meta$Comparator,
                           CA19_9 = meta$CA19_9, Outer_Fold = fold)
cv_results <- list()
tuning <- list()
lambda_grid <- 10^seq(2, -4, length.out = 60)
for (model in names(panels)) {
  x <- cbind(base, t(z$x[panels[[model]], keep, drop = FALSE]))
  score <- numeric(length(y))
  for (i in seq_len(10)) {
    train <- fold != i
    inner <- stratified_folds(y[train], 5, 20261004 + i)
    fit <- glmnet::cv.glmnet(x[train, , drop = FALSE], y[train], family = "binomial",
                            alpha = 0, lambda = lambda_grid, foldid = inner,
                            type.measure = "deviance", standardize = TRUE, parallel = FALSE)
    score[!train] <- as.numeric(predict(fit, newx = x[!train, , drop = FALSE],
                                       s = "lambda.min", type = "response"))
    tuning[[length(tuning) + 1]] <- data.frame(Model = model, Outer_Fold = i,
                                               Lambda = fit$lambda.min, Train_N = sum(train), Test_N = sum(!train))
  }
  stopifnot(all(is.finite(score)), all(score >= 0 & score <= 1))
  predictions[[model]] <- score
  r <- roc_stats(y, score)
  cv_results[[model]] <- data.frame(Model = model, N = length(y), Cases = sum(y), Controls = sum(!y),
                                    AUC = r["AUC"], CI_low = r["CI_low"], CI_high = r["CI_high"],
                                    Brier = mean((score - y)^2))
}
cv <- do.call(rbind, cv_results)
comparison <- do.call(rbind, lapply(names(panels)[-1], function(model) {
  a <- pROC::roc(y, predictions[[model]], direction = "<", quiet = TRUE)
  b <- pROC::roc(y, predictions$CA19_9_age_sex, direction = "<", quiet = TRUE)
  test <- pROC::roc.test(a, b, method = "delong", paired = TRUE)
  data.frame(Model = model, Delta_AUC_vs_Baseline = as.numeric(pROC::auc(a) - pROC::auc(b)),
             Conditional_DeLong_P = test$p.value)
}))
comparison$BH_two_comparisons <- p.adjust(comparison$Conditional_DeLong_P, "BH")
write.csv(cv, file.path(external_dir, "GSE85589_CA19_9_nested_CV_performance.csv"), row.names = FALSE)
write.csv(predictions, file.path(external_dir, "GSE85589_CA19_9_out_of_fold_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, tuning), file.path(external_dir, "GSE85589_CA19_9_inner_lambda.csv"), row.names = FALSE)
write.csv(comparison, file.path(external_dir, "GSE85589_CA19_9_conditional_comparison.csv"), row.names = FALSE)
print(cv, row.names = FALSE)
print(comparison, row.names = FALSE)
cat("Marker coverage:\n")
print(table(do.call(rbind, coverage)$Dataset, do.call(rbind, coverage)$Measurable_Complete))
