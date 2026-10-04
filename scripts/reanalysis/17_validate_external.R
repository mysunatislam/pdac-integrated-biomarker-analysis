#!/usr/bin/env Rscript
source("scripts/reanalysis/external_helpers.R")
read_result <- function(name) read.csv(file.path(external_dir, name))
reference <- read_result("frozen_spatial_reference.csv")
tissue <- read_result("tissue_all_gene_contrasts.csv")
family <- read_result("tissue_frozen_spatial_family.csv")
models <- read_result("tissue_models.csv")
stopifnot(nrow(reference) == 506, !anyDuplicated(reference$Gene),
          all(family$Gene %in% reference$Gene),
          !anyDuplicated(paste(tissue$Dataset, tissue$Comparison, tissue$Gene)),
          all(is.finite(tissue$P)), all(tissue$P >= 0 & tissue$P <= 1),
          all(tissue$CI_low <= tissue$log2FC), all(tissue$CI_high >= tissue$log2FC),
          identical(family$Same_direction, sign(family$log2FC) == sign(family$Spatial_log2FC)))
for (s in split(tissue, paste(tissue$Dataset, tissue$Comparison))) {
  stopifnot(isTRUE(all.equal(s$BH_all_genes, p.adjust(s$P, "BH"))))
  model <- models[models$Dataset == s$Dataset[1] & models$Comparison == s$Comparison[1], ]
  stopifnot(nrow(model) == 1, model$Genes_Tested == nrow(s),
            model$All_Gene_BH_Significant == sum(s$BH_all_genes < .05))
}
for (s in split(family, paste(family$Dataset, family$Comparison))) {
  stopifnot(isTRUE(all.equal(s$BH_frozen_family, p.adjust(s$P, "BH"))))
  original <- reference[match(s$Gene, reference$Gene), ]
  stopifnot(identical(s$Spatial_log2FC, original[[paste0(s$Spatial_Contrast[1], "_log2")]]))
}
paired <- read_result("GSE15471_donor_tissue_mapping.csv")
adm <- read_result("GSE179248_analysis_mapping.csv")
stopifnot(nrow(paired) == 70, all(table(paired$Donor, paired$State) == 1),
          !any(paired$Donor == "51721"), sum(paired$Array_Count == 2) == 6,
          nrow(adm) == 28, all(table(adm$Donor, adm$Day) == 1),
          !anyDuplicated(adm$GSM), !anyDuplicated(adm$Expression_Column),
          all(adm$TMM_Factor > 0), all(adm$Filtered_Library_Size > 0))
confound <- read_result("GSE101462_eligibility_decision.csv")
stopifnot(confound$Design_Rank < confound$Design_Columns,
          !"GSE101462" %in% tissue$Dataset, nrow(models) == 5)
mapping <- read_result("GSE91035_platform_annotation_audit.csv")
selected <- read_result("GSE91035_selected_probe_map.csv")
stopifnot(all(selected$Gene %in% AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")),
          !anyDuplicated(selected$Gene), !anyDuplicated(selected$Probe),
          !anyNA(mapping$Probe))

serum <- read_result("serum_all_miRNA_contrasts.csv")
family <- read_result("serum_frozen_miRNA_family.csv")
markers <- read_result("frozen_miRNA_reference.csv")
aucs <- read_result("serum_fixed_direction_marker_AUC.csv")
stopifnot(nrow(markers) == 27, all(family$miRNA %in% markers$miRNA),
          all(is.finite(serum$P)), !anyDuplicated(paste(serum$Dataset, serum$Comparison, serum$miRNA)),
          identical(family$Same_Direction_As_Plasma, sign(family$log2FC) == family$Fixed_Direction),
          all(aucs$AUC >= 0 & aucs$AUC <= 1), all(aucs$CI_low <= aucs$AUC),
          all(aucs$CI_high >= aucs$AUC), all(aucs$Cases > 0), all(aucs$Controls > 0))
for (s in split(serum, paste(serum$Dataset, serum$Comparison))) {
  stopifnot(isTRUE(all.equal(s$BH_all_genes, p.adjust(s$P, "BH"))))
}
for (s in split(family, paste(family$Dataset, family$Comparison))) {
  stopifnot(isTRUE(all.equal(s$BH_frozen_within_contrast, p.adjust(s$P, "BH"))))
}
for (s in split(family, family$Dataset)) {
  stopifnot(isTRUE(all.equal(s$BH_frozen_all_comparisons, p.adjust(s$P, "BH"))))
}
cv <- read_result("GSE85589_CA19_9_nested_CV_performance.csv")
scores <- read_result("GSE85589_CA19_9_out_of_fold_scores.csv")
tuning <- read_result("GSE85589_CA19_9_inner_lambda.csv")
stopifnot(nrow(scores) == 115, sum(scores$Case) == 88, !anyDuplicated(scores$GSM),
          nrow(cv) == 3, nrow(tuning) == 30, all(table(scores$Outer_Fold, scores$Case) > 0))
for (model in cv$Model) {
  row <- cv[cv$Model == model, ]
  auc <- as.numeric(pROC::auc(pROC::roc(scores$Case, scores[[model]], direction = "<", quiet = TRUE)))
  stopifnot(abs(auc - row$AUC) < 1e-12,
            abs(mean((scores[[model]] - scores$Case)^2) - row$Brier) < 1e-12,
            all(is.finite(scores[[model]])))
}
stage <- aucs[aucs$Subgroup == "Pathologic_stage_II_non_yp", ]
ca_low <- aucs[aucs$Subgroup == "CA19_9_less_than_37", ]
stopifnot(nrow(stage) > 0, all(stage$Cases == 17), all(ca_low$Cases == 23),
          all(ca_low$Controls[ca_low$Comparison == "PC_vs_Benign"] == 7),
          all(ca_low$Controls[ca_low$Comparison == "PC_vs_Healthy"] == 19))
manifest <- read_result("download_manifest.csv")
committed <- manifest[!grepl("family.soft.gz$", manifest$File), ]
hashes <- unname(tools::md5sum(file.path("data/external", committed$File)))
stopifnot(identical(hashes, committed$MD5))
figures <- list.files("reanalysis/figures/external", pattern = "\\.pdf$", full.names = TRUE)
stopifnot(length(figures) == 4, all(file.exists(sub("\\.pdf$", ".png", figures))))
for (path in list.files("scripts/reanalysis", pattern = "\\.R$", full.names = TRUE)) parse(path)
writeLines(c("External output checks: PASS", "35 paired tissue donors; 14 paired culture donors",
              "CP n=2 kept exploratory; perfect-confounding cohort not fit",
              "27 locked serum markers; ROC direction fixed to diagnostic plasma",
              "115 serum CV participants; three models, shared outer folds, 30 inner lambda selections",
              "BH, CIs, cohort sizes, subgroup labels, input hashes and script parsing checked"),
            "reanalysis/provenance/external_output_checks.txt")
cat("External analysis checks: PASS\n")
