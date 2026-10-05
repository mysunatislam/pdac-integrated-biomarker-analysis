#!/usr/bin/env Rscript

source("scripts/reanalysis/external_helpers.R")
out <- "reanalysis/results/publication"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
markers <- read.delim("data/publication/MCPcounter_genes.txt", check.names = FALSE)
stopifnot(!("STAT1" %in% markers[["HUGO symbols"]]),
          utils::packageDescription("MCPcounter")$RemoteSha ==
            "b6eac73e91c246fcff0bb1a5c68a816cd588fc48")

roi <- read.csv("reanalysis/results/input_audit/GSE208536_roi_mapping.csv")
published <- read.csv("data/publication/GSE208536_published_cases.csv")
clinical <- grep("^Sample_characteristics_ch1", names(roi), value = TRUE)[-1]
parts <- lapply(roi[clinical], function(x) toupper(trimws(sub("^[^:]+:", "", x))))
key <- do.call(paste, c(parts, sep = "|"))
published_key <- with(published, paste(Age, Sex, Grade, Stage, sep = "|"))
stopifnot(!anyDuplicated(published_key), length(unique(key)) == 8,
          setequal(unique(key), published_key))
roi$Published_Case <- published$Published_Case[match(key, published_key)]
roi$Inference_Status <- "Unique clinical-attribute crosswalk; not a directly deposited patient ID"
roi$Original_Profile <- sprintf("P%02d", match(roi$Clinical_Profile, unique(roi$Clinical_Profile)))
write.csv(roi, file.path(out, "GSE208536_published_case_crosswalk.csv"), row.names = FALSE)
case_counts <- as.data.frame.matrix(table(roi$Published_Case, roi$Condition))
case_counts$Published_Case <- rownames(case_counts)
write.csv(case_counts, file.path(out, "GSE208536_published_case_roi_counts.csv"), row.names = FALSE)

coverage_list <- score_list <- model_list <- list()
score_cohort <- function(x, meta, id) {
  stopifnot(identical(colnames(x), meta$Sample), all(is.finite(x)))
  scores <- MCPcounter::MCPcounter.estimate(x, featuresType = "HUGO_symbols", genes = markers)
  groups <- split(markers[["HUGO symbols"]], markers[["Cell population"]])
  coverage <- do.call(rbind, lapply(names(groups), function(cell) {
    available <- intersect(unique(groups[[cell]]), rownames(x))
    sd_score <- if (cell %in% rownames(scores)) sd(scores[cell, ]) else NA_real_
    data.frame(Dataset = id, Cell_Population = cell, Published_Markers = length(unique(groups[[cell]])),
               Present_Markers = length(available), Coverage = length(available) / length(unique(groups[[cell]])),
               Eligible_Adjustment = length(available) >= 3 &&
                 length(available) / length(unique(groups[[cell]])) >= 0.5 &&
                 is.finite(sd_score) && sd_score > 0, Markers = paste(available, collapse = ";"))
  }))
  coverage_list[[id]] <<- coverage
  score_list[[id]] <<- cbind(Dataset = id, meta, as.data.frame(t(scores), check.names = FALSE))
  for (cell in c("T cells", "Fibroblasts")) {
    name <- if (cell == "T cells") "T_Score" else "Fibroblast_Score"
    eligible <- coverage$Eligible_Adjustment[coverage$Cell_Population == cell]
    if (length(eligible) && eligible) meta[[name]] <- as.numeric(scale(scores[cell, ]))
  }
  meta
}

model_cohort <- function(x, meta, id, base, contrasts, spatial = FALSE, rna_seq = NULL) {
  variants <- list(Unadjusted_for_composition = character(), Plus_T_cells = "T_Score",
                   Plus_Fibroblasts = "Fibroblast_Score", Plus_both = c("T_Score", "Fibroblast_Score"))
  for (variant in names(variants)) {
    covariates <- variants[[variant]]
    reason <- ""
    if (!all(covariates %in% names(meta))) reason <- "Required score failed marker coverage or variance rule"
    formula <- reformulate(c(base, covariates))
    if (!nzchar(reason)) {
      design <- model.matrix(formula, meta)
      if (qr(design)$rank < ncol(design)) reason <- "Rank-deficient design"
      if (nrow(design) - qr(design)$rank <= 0) reason <- "No residual degrees of freedom"
    }
    for (label in names(contrasts)) {
      if (nzchar(reason)) {
        result <- data.frame(Dataset = id, Comparison = label, Variant = variant,
                             Status = reason, Samples = ncol(x), log2FC = NA_real_,
                             SE = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
                             BH_all_genes = NA_real_, Residual_DF = NA_real_, Condition_Number = NA_real_)
      } else {
        l <- setNames(rep(0, ncol(design)), colnames(design))
        l[names(contrasts[[label]])] <- contrasts[[label]]
        stopifnot(all(names(contrasts[[label]]) %in% colnames(design)))
        if (spatial) {
          fit <- lm.fit(design, as.numeric(x["STAT1", ]))
          sigma2 <- sum(fit$residuals^2) / fit$df.residual
          covariance <- chol2inv(chol(crossprod(design))) * sigma2
          beta <- sum(l * fit$coefficients)
          se <- sqrt(as.numeric(t(l) %*% covariance %*% l))
          df <- fit$df.residual
          p <- 2 * pt(-abs(beta / se), df)
          bh <- NA_real_
        } else {
          fit <- if (is.null(rna_seq)) limma::lmFit(x, design) else
            limma::lmFit(limma::voom(rna_seq, design, plot = FALSE), design)
          fit <- limma::eBayes(limma::contrasts.fit(fit, matrix(l, ncol = 1,
                                     dimnames = list(names(l), label))),
                              trend = is.null(rna_seq), robust = TRUE)
          i <- match("STAT1", rownames(x))
          beta <- fit$coefficients[i, 1]
          se <- fit$stdev.unscaled[i, 1] * sqrt(fit$s2.post[i])
          df <- fit$df.total[i]
          p <- fit$p.value[i, 1]
          bh <- p.adjust(fit$p.value[, 1], "BH")[i]
        }
        margin <- qt(0.975, df) * se
        result <- data.frame(Dataset = id, Comparison = label, Variant = variant, Status = "Estimated",
                             Samples = ncol(x), log2FC = beta, SE = se, CI_low = beta - margin,
                             CI_high = beta + margin, P = p, BH_all_genes = bh,
                             Residual_DF = nrow(design) - ncol(design), Condition_Number = kappa(design))
      }
      model_list[[paste(id, variant, label, sep = "_")]] <<- result
    }
  }
}

book <- readxl::read_excel("data/source/GSE208536_Processed_data.xlsx", sheet = "Sheet1")
x <- log2(as.matrix(book[, -1]) + 1)
rownames(x) <- book[[1]]
stopifnot(identical(colnames(x), roi$Expression_Column))
groups <- split(seq_len(nrow(roi)), interaction(roi$Original_Profile, roi$Condition, drop = TRUE))
xs <- vapply(groups, function(ii) rowMeans(x[, ii, drop = FALSE]), numeric(nrow(x)))
rownames(xs) <- rownames(x)
ms <- do.call(rbind, lapply(names(groups), function(name) {
  ii <- groups[[name]]
  data.frame(Sample = name, Profile = roi$Original_Profile[ii[1]], State = roi$Condition[ii[1]],
             Published_Case = roi$Published_Case[ii[1]], ROIs = length(ii))
}))
ms$Profile <- factor(ms$Profile)
ms$State <- factor(ms$State, levels = c("Normal", "ADM", "PDAC"))
ms <- score_cohort(xs, ms, "GSE208536")
model_cohort(xs, ms, "GSE208536", c("Profile", "State"),
             list(ADM_vs_Normal = c(StateADM = 1), PDAC_vs_ADM = c(StatePDAC = 1, StateADM = -1),
                  PDAC_vs_Normal = c(StatePDAC = 1)), spatial = TRUE)

g <- read_external_geo("GSE15471")
probe_map <- read.csv("reanalysis/results/external/GSE15471_selected_probe_map.csv")
stopifnot(all(c("Gene", "Probe") %in% names(probe_map)))
stopifnot(all(probe_map$Probe %in% rownames(g$x)))
xp <- g$x[match(probe_map$Probe, rownames(g$x)), , drop = FALSE]
rownames(xp) <- probe_map$Gene
ma <- read.csv("reanalysis/results/external/GSE15471_analysis_mapping.csv")
ma <- ma[match(colnames(xp), ma$GSM), ]
included <- which(ma$Included)
groups <- split(included, paste(ma$Donor[included], ma$State[included], sep = "_"))
xp <- vapply(groups, function(ii) rowMeans(xp[, ii, drop = FALSE]), numeric(nrow(xp)))
rownames(xp) <- probe_map$Gene
mp <- do.call(rbind, lapply(names(groups), function(name) {
  ii <- groups[[name]]
  data.frame(Sample = name, Donor = ma$Donor[ii[1]], State = ma$State[ii[1]])
}))
mp$Donor <- factor(mp$Donor)
mp$State <- factor(mp$State, levels = c("Normal", "PDAC"))
mp <- score_cohort(xp, mp, "GSE15471")
model_cohort(xp, mp, "GSE15471", c("Donor", "State"), list(PDAC_vs_Normal = c(StatePDAC = 1)))

g <- read_external_geo("GSE143754", path = "data/source/GSE143754_series_matrix.txt.gz")
probe_map <- read.csv("reanalysis/results/bulk/GSE143754_gene_level_contrasts.csv")
xb <- g$x[match(probe_map$Probe, rownames(g$x)), , drop = FALSE]
rownames(xb) <- probe_map$Gene
group <- characteristic(g$meta, "disease state:")
ii <- which(group %in% c("Chronic Pancreatitis", "Tumor"))
xb <- xb[, ii, drop = FALSE]
mb <- data.frame(Sample = colnames(xb),
                 State = factor(ifelse(group[ii] == "Tumor", "PDAC", "CP"), levels = c("CP", "PDAC")),
                 Age = as.numeric(characteristic(g$meta, "age:")[ii]),
                 Sex = factor(characteristic(g$meta, "Sex:")[ii]))
mb <- score_cohort(xb, mb, "GSE143754")
model_cohort(xb, mb, "GSE143754", c("State", "Age", "Sex"), list(PDAC_vs_CP = c(StatePDAC = 1)))

rm(book, x, xs, xp, xb, g)
invisible(gc())
book <- readxl::read_excel("data/external/GSE179248_CountReads_D6_0_processed.xlsx",
                          sheet = "CountReads_D6_0_processed", .name_repair = "minimal")
counts <- as.matrix(book[, -1])
storage.mode(counts) <- "double"
ensembl <- as.character(book[[1]])
annotation <- read.csv("reanalysis/results/external/GSE179248_ensembl_annotation.csv")
symbols <- annotation$Gene[match(ensembl, annotation$ENSEMBL)]
valid <- !is.na(symbols) & nzchar(symbols)
counts <- rowsum(counts[valid, , drop = FALSE], symbols[valid], reorder = TRUE)
mc <- read.csv("reanalysis/results/external/GSE179248_analysis_mapping.csv")
mc <- mc[match(colnames(counts), mc$Expression_Column), ]
stopifnot(identical(mc$Expression_Column, colnames(counts)), nrow(mc) == 28,
          all(table(mc$Donor, mc$Day) == 1))
mc <- data.frame(Sample = colnames(counts), Donor = factor(mc$Donor),
                 State = factor(mc$Day, levels = c("D0", "D6")))
y <- edgeR::DGEList(counts)
keep <- edgeR::filterByExpr(y, group = mc$State)
y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
v <- limma::voom(y, model.matrix(~Donor + State, mc), plot = FALSE)
mc <- score_cohort(v$E, mc, "GSE179248")
model_cohort(v$E, mc, "GSE179248", c("Donor", "State"),
             list(Culture_D6_vs_D0 = c(StateD6 = 1)), rna_seq = y)

write.csv(do.call(rbind, coverage_list), file.path(out, "MCPcounter_marker_coverage.csv"), row.names = FALSE)
for (id in names(score_list)) write.csv(score_list[[id]], file.path(out, paste0(id, "_composition_scores.csv")), row.names = FALSE)
results <- do.call(rbind, model_list)
results$BH_exploratory_model_family <- p.adjust(results$P, "BH")
write.csv(results, file.path(out, "STAT1_composition_sensitivity.csv"), row.names = FALSE)
print(results[, c("Dataset", "Comparison", "Variant", "log2FC", "P", "Status")], row.names = FALSE)
cat("Eight unique published cases matched; composition sensitivity tables saved.\n")
