#!/usr/bin/env Rscript

source("scripts/reanalysis/external_helpers.R")
out <- "reanalysis/results/publication"
tab <- read.csv(gzfile("data/publication/E-MTAB-1791_primary_normalized_genes.csv.gz"),
                check.names = FALSE, colClasses = c("character", rep("numeric", 254)))
x <- as.matrix(tab[, -1])
rownames(x) <- tab$Gene
meta <- read.csv(file.path(out, "E-MTAB-1791_eligibility_mapping.csv"))
meta <- meta[match(colnames(x), meta$Specimen), ]
meta$State <- factor(meta$Role, levels = c("CP", "PDAC"))
meta$Chip <- factor(meta$Chip)
stopifnot(ncol(x) == 254, all(is.finite(x)), !anyDuplicated(rownames(x)),
          identical(colnames(x), meta$Specimen), all(meta$Primary_Eligible),
          sum(meta$State == "PDAC") == 195, sum(meta$State == "CP") == 59)
chips <- as.data.frame.matrix(table(meta$Chip, meta$State))
chips$Chip <- rownames(chips)
chips$Both_Diagnoses <- chips$CP > 0 & chips$PDAC > 0
write.csv(chips, file.path(out, "E-MTAB-1791_chip_balance.csv"), row.names = FALSE)

markers <- read.delim("data/publication/MCPcounter_genes.txt", check.names = FALSE)
stopifnot(!"STAT1" %in% markers[["HUGO symbols"]])
scores <- MCPcounter::MCPcounter.estimate(x, featuresType = "HUGO_symbols", genes = markers)
groups <- split(markers[["HUGO symbols"]], markers[["Cell population"]])
coverage <- do.call(rbind, lapply(names(groups), function(cell) {
  available <- intersect(unique(groups[[cell]]), rownames(x))
  variance <- if (cell %in% rownames(scores)) sd(scores[cell, ]) else NA_real_
  data.frame(Dataset = "E-MTAB-1791", Cell_Population = cell,
             Published_Markers = length(unique(groups[[cell]])), Present_Markers = length(available),
             Coverage = length(available) / length(unique(groups[[cell]])),
             Eligible_Adjustment = length(available) >= 3 &&
               length(available) / length(unique(groups[[cell]])) >= .5 &&
               is.finite(variance) && variance > 0, Markers = paste(available, collapse = ";"))
}))
for (cell in c("T cells", "Fibroblasts")) {
  eligible <- coverage$Eligible_Adjustment[coverage$Cell_Population == cell]
  if (length(eligible) && eligible)
    meta[[if (cell == "T cells") "T_Score" else "Fibroblast_Score"]] <- as.numeric(scale(scores[cell, ]))
}
write.csv(coverage, file.path(out, "E-MTAB-1791_MCPcounter_marker_coverage.csv"), row.names = FALSE)
write.csv(cbind(meta, as.data.frame(t(scores), check.names = FALSE)),
          file.path(out, "E-MTAB-1791_composition_scores.csv"), row.names = FALSE)

fit_model <- function(expression, metadata, terms) {
  d <- model.matrix(reformulate(terms), droplevels(metadata))
  if (qr(d)$rank != ncol(d) || nrow(d) - qr(d)$rank <= 0)
    return(list(status = "Diagnosis not estimable or no residual degrees of freedom",
                design = d))
  fit <- limma::eBayes(limma::lmFit(expression, d), trend = TRUE, robust = TRUE)
  stats <- limma_result(fit, "StatePDAC", "E-MTAB-1791", "PDAC_vs_primary_CP")
  list(status = "Estimated", stats = stats, design = d)
}
base_models <- list(Unadjusted = "State", Chip_adjusted = c("Chip", "State"))
variants <- list(Unadjusted_for_composition = character(), Plus_T_cells = "T_Score",
                 Plus_Fibroblasts = "Fibroblast_Score", Plus_both = c("T_Score", "Fibroblast_Score"))
all_results <- frozen_results <- model_results <- stat1_results <- list()
for (base in names(base_models)) {
  for (variant in names(variants)) {
    terms <- c(base_models[[base]], variants[[variant]])
    fit <- if (!all(variants[[variant]] %in% names(meta)))
      list(status = "Required score failed marker coverage or variance rule") else fit_model(x, meta, terms)
    cat(base, "/", variant, ":", fit$status, "\n")
    flush.console()
    d <- fit$design
    model_results[[paste(base, variant)]] <- data.frame(
      Model = base, Variant = variant, Status = fit$status, Specimens = ncol(x),
      PDAC_Specimens = 195, CP_Specimens = 59, Deposited_Patient_ID_Available = FALSE,
      Age_Sex_Adjustment_Available = FALSE, Chips = nlevels(meta$Chip),
      Mixed_Diagnosis_Chips = sum(chips$Both_Diagnoses),
      Design_Columns = if (is.null(d)) NA_integer_ else ncol(d),
      Design_Rank = if (is.null(d)) NA_integer_ else qr(d)$rank,
      Residual_DF = if (is.null(d)) NA_integer_ else nrow(d) - qr(d)$rank,
      Condition_Number = if (is.null(d)) NA_real_ else kappa(d))
    if (fit$status != "Estimated") next
    stats <- fit$stats
    stats$Model <- base
    stats$Variant <- variant
    stat1_results[[paste(base, variant)]] <- stats[stats$Gene == "STAT1", ]
    if (variant == "Unadjusted_for_composition") {
      all_results[[base]] <- stats
      frozen_results[[base]] <- attach_reference(stats, "Normal_to_PDAC")
    }
  }
}
all_stats <- do.call(rbind, all_results)
frozen <- do.call(rbind, frozen_results)
models <- do.call(rbind, model_results)
stat1 <- do.call(rbind, stat1_results)
stat1$BH_exploratory_model_family <- p.adjust(stat1$P, "BH")
write.csv(all_stats, file.path(out, "E-MTAB-1791_all_gene_contrasts.csv"), row.names = FALSE)
write.csv(frozen, file.path(out, "E-MTAB-1791_frozen_spatial_family.csv"), row.names = FALSE)
write.csv(models, file.path(out, "E-MTAB-1791_models.csv"), row.names = FALSE)
write.csv(stat1, file.path(out, "E-MTAB-1791_STAT1_composition_sensitivity.csv"), row.names = FALSE)

# One of 59 CP specimens is an unidentified repeat in the source publication.
# Omitting every possible specimen tests sensitivity without inventing a donor map.
omissions <- omit_family <- list()
reference <- read.csv(file.path(external_dir, "frozen_spatial_reference.csv"))
cp_index <- which(meta$State == "CP")
for (base in names(base_models)) {
  if (!base %in% names(all_results)) next
  for (k in seq_along(cp_index)) {
    j <- cp_index[k]
    fit <- fit_model(x[, -j, drop = FALSE], meta[-j, , drop = FALSE], base_models[[base]])
    if (fit$status != "Estimated") {
      omissions[[paste(base, j)]] <- data.frame(Model = base, Omitted_CP = meta$Specimen[j],
                                                Status = fit$status, log2FC = NA_real_,
                                                CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
                                                BH_all_genes = NA_real_, Frozen_Significant = NA_integer_)
      next
    }
    s <- fit$stats[fit$stats$Gene == "STAT1", ]
    family <- fit$stats[fit$stats$Gene %in% reference$Gene, ]
    family$BH_frozen_family <- p.adjust(family$P, "BH")
    family$Model <- base
    family$Omitted_CP <- meta$Specimen[j]
    omit_family[[paste(base, j)]] <- family
    omissions[[paste(base, j)]] <- data.frame(Model = base, Omitted_CP = meta$Specimen[j],
                                             Status = fit$status, log2FC = s$log2FC,
                                             CI_low = s$CI_low, CI_high = s$CI_high, P = s$P,
                                             BH_all_genes = s$BH_all_genes,
                                             Frozen_Significant = sum(family$BH_frozen_family < .05))
    if (k %% 10 == 0) {
      cat(base, ":", k, "of", length(cp_index), "CP omissions completed\n")
      flush.console()
    }
  }
  cat(base, ": completed all 59 CP omissions\n")
}
omit <- do.call(rbind, omissions)
write.csv(omit, file.path(out, "E-MTAB-1791_omit_one_CP_STAT1.csv"), row.names = FALSE)
omitted_family <- do.call(rbind, omit_family)
write.csv(omitted_family, file.path(out, "E-MTAB-1791_omit_one_CP_frozen_family.csv"), row.names = FALSE)
envelope <- do.call(rbind, lapply(split(omitted_family, paste(omitted_family$Model, omitted_family$Gene)), function(s) {
  original <- frozen[frozen$Model == s$Model[1] & frozen$Gene == s$Gene[1], ]
  data.frame(Model = s$Model[1], Gene = s$Gene[1], Estimable_Omissions = nrow(s),
             Min_log2FC = min(s$log2FC), Max_log2FC = max(s$log2FC),
             Min_BH_frozen = min(s$BH_frozen_family), Max_BH_frozen = max(s$BH_frozen_family),
             Significant_All_Estimable_Omissions = all(s$BH_frozen_family < .05),
             Direction_Stable = all(sign(s$log2FC) == sign(original$log2FC)))
}))
write.csv(envelope, file.path(out, "E-MTAB-1791_omit_one_CP_envelope.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(frozen, frozen$Model), function(s) {
  a <- all_stats[all_stats$Model == s$Model[1], ]
  e <- envelope[envelope$Model == s$Model[1], ]
  data.frame(Model = s$Model[1], All_Genes_Tested = nrow(a), All_Gene_BH_Significant = sum(a$BH_all_genes < .05),
             Measurable_Frozen = nrow(s), Frozen_BH_Significant = sum(s$BH_frozen_family < .05),
             Same_Direction_As_Spatial = sum(s$Same_direction),
             Spatial_And_CP_Concordant_Significant = sum(s$Spatial_BH < .05 & s$BH_all_genes < .05 & s$Same_direction),
             Frozen_Significant_All_Omissions = sum(e$Significant_All_Estimable_Omissions),
             Spearman_Effect_Correlation = cor(s$log2FC, s$Spatial_log2FC, method = "spearman"))
}))
write.csv(summary, file.path(out, "E-MTAB-1791_replication_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
print(stat1[, c("Model", "Variant", "log2FC", "CI_low", "CI_high", "P", "BH_all_genes")], row.names = FALSE)
cat("Larger inflammatory cohort association and sensitivity analyses saved.\n")
