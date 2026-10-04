#!/usr/bin/env Rscript
source("scripts/reanalysis/external_helpers.R")
tissue_results <- list()
reference_results <- list()
models <- list()
save_tissue <- function(stats, spatial_contrast, cases, controls, donors, role) {
  tissue_results[[length(tissue_results) + 1]] <<- stats
  reference_results[[length(reference_results) + 1]] <<- attach_reference(stats, spatial_contrast)
  models[[length(models) + 1]] <<- data.frame(
    Dataset = stats$Dataset[1], Comparison = stats$Comparison[1], Cases = cases,
    Controls = controls, Paired_Donors = donors, Genes_Tested = nrow(stats),
    All_Gene_BH_Significant = sum(stats$BH_all_genes < .05), Role = role
  )
}

cat("GSE15471: donor-blocked, technical replicates averaged\n")
z <- read_external_geo("GSE15471")
z$meta$Donor <- characteristic(z$meta, "patient:")
z$meta$State <- ifelse(characteristic(z$meta, "sample:") == "tumor", "PDAC", "Normal")
z$meta$Included <- z$meta$Donor != "51721"
z$meta$Exclusion <- ifelse(z$meta$Included, "", "Published failed-QC donor pair")
write.csv(z$meta, file.path(external_dir, "GSE15471_analysis_mapping.csv"), row.names = FALSE)
indices <- which(z$meta$Included)
key <- paste(z$meta$Donor[indices], z$meta$State[indices], sep = "_")
keys <- sort(unique(key))
x <- vapply(keys, function(k) rowMeans(z$x[, indices[key == k], drop = FALSE]),
            numeric(nrow(z$x)))
rownames(x) <- rownames(z$x)
meta <- z$meta[indices[match(keys, key)], c("Donor", "State")]
meta$Array_Count <- as.integer(table(factor(key, levels = keys)))
meta$Key <- keys
stopifnot(nrow(meta) == 70, all(table(meta$Donor, meta$State) == 1),
          length(unique(meta$Donor)) == 35, sum(meta$Array_Count == 2) == 6)
write.csv(meta, file.path(external_dir, "GSE15471_donor_tissue_mapping.csv"), row.names = FALSE)
db <- hgu133plus2.db::hgu133plus2.db
symbols <- AnnotationDbi::mapIds(db, keys = rownames(x), keytype = "PROBEID",
                               column = "SYMBOL", multiVals = "asNA")
gene <- collapse_probes(x, unname(symbols[rownames(x)]))
write.csv(data.frame(Gene = rownames(gene$x), Probe = gene$probes),
          file.path(external_dir, "GSE15471_selected_probe_map.csv"), row.names = FALSE)
meta$State <- factor(meta$State, levels = c("Normal", "PDAC"))
design <- model.matrix(~factor(Donor) + State, meta)
stopifnot(qr(design)$rank == ncol(design))
fit <- limma::eBayes(limma::lmFit(gene$x, design), trend = TRUE, robust = TRUE)
stats <- limma_result(fit, "StatePDAC", "GSE15471", "PDAC_vs_paired_Normal")
save_tissue(stats, "Normal_to_PDAC", 35, 35, 35, "Paired tumor/adjacent-normal association")

cat("GSE91035: independent two-group comparisons, CP n=2\n")
z <- read_external_geo("GSE91035")
b <- read.delim(gzfile("data/external/GSE91035_normalized_data_with_gene_symbol.txt.gz"),
                check.names = FALSE, quote = "\"")
probe_index <- match(rownames(z$x), b$ID_REF)
sample_index <- match(sub(",.*", "", z$meta$Title), names(b))
stopifnot(!anyNA(probe_index), !anyNA(sample_index))
stopifnot(max(abs(as.matrix(b[probe_index, sample_index]) - z$x)) < 1e-6)
platform <- read.csv("data/external/GPL22763_annotation.csv", check.names = FALSE)
stopifnot(!anyDuplicated(platform$ID))
platform <- platform[match(rownames(z$x), platform$ID), ]
stopifnot(identical(platform$ID, rownames(z$x)))
db <- org.Hs.eg.db::org.Hs.eg.db
entrez <- as.character(platform$LOCUSLINK_ID)
known <- intersect(entrez, AnnotationDbi::keys(db, keytype = "ENTREZID"))
mapped <- AnnotationDbi::mapIds(db, keys = known, keytype = "ENTREZID",
                              column = "SYMBOL", multiVals = "asNA")
symbol <- unname(mapped[entrez])
deposited_symbol <- trimws(platform$GENE_SYMBOL)
exact <- deposited_symbol %in% AnnotationDbi::keys(db, keytype = "SYMBOL")
symbol[is.na(symbol) & exact] <- deposited_symbol[is.na(symbol) & exact]
unknown <- is.na(symbol) & !exact
alias_keys <- intersect(deposited_symbol[unknown], AnnotationDbi::keys(db, keytype = "ALIAS"))
aliases <- AnnotationDbi::mapIds(db, keys = alias_keys, keytype = "ALIAS",
                               column = "SYMBOL", multiVals = "asNA")
symbol[unknown] <- unname(aliases[deposited_symbol[unknown]])
write.csv(data.frame(Probe = rownames(z$x), Supplement_Label = b$GeneSymbol[probe_index],
                    Platform_Label = deposited_symbol, Entrez_ID = entrez, Canonical_Gene = symbol),
          file.path(external_dir, "GSE91035_platform_annotation_audit.csv"), row.names = FALSE)
gene <- collapse_probes(z$x, symbol)
z$meta$State <- characteristic(z$meta, "disease state:")
z$meta$State[z$meta$State == "Chronic Pancreatitis"] <- "CP"
write.csv(z$meta, file.path(external_dir, "GSE91035_analysis_mapping.csv"), row.names = FALSE)
write.csv(data.frame(Gene = rownames(gene$x), Probe = gene$probes),
          file.path(external_dir, "GSE91035_selected_probe_map.csv"), row.names = FALSE)
stopifnot(identical(as.integer(table(factor(z$meta$State,
                       levels = c("normal", "benign", "CP", "PDAC")))), c(8L, 15L, 2L, 25L)))
for (control in c("CP", "normal", "benign")) {
  keep <- z$meta$State %in% c(control, "PDAC")
  state <- factor(z$meta$State[keep], levels = c(control, "PDAC"))
  design <- model.matrix(~state)
  fit <- limma::eBayes(limma::lmFit(gene$x[, keep, drop = FALSE], design),
                       trend = TRUE, robust = TRUE)
  comparison <- paste0("PDAC_vs_", control)
  stats <- limma_result(fit, "statePDAC", "GSE91035", comparison)
  role <- switch(control,
                 CP = "Exploratory: only two CP controls, no clinical adjustment",
                 normal = "Exploratory: commercially sourced normal RNA",
                 benign = "Exploratory: adjacent benign tissue, pairing unavailable")
  save_tissue(stats, "Normal_to_PDAC", sum(state == "PDAC"), sum(state == control), NA, role)
}

cat("GSE179248: paired donor RNA-seq culture day 6 versus day 0\n")
z <- read_external_geo("GSE179248", expression = FALSE)
z$meta$Donor <- sub("^donor ", "", characteristic(z$meta, "individual:"))
z$meta$Day <- sub("^ADM_Day", "D", characteristic(z$meta, "time:"))
z$meta$Expression_Column <- paste0("s_", z$meta$Donor, z$meta$Day)
b <- readxl::read_excel("data/external/GSE179248_CountReads_D6_0_processed.xlsx",
                        sheet = "CountReads_D6_0_processed", .name_repair = "minimal")
meta <- z$meta[match(names(b)[-1], z$meta$Expression_Column), ]
stopifnot(identical(meta$Expression_Column, names(b)[-1]),
          all(table(meta$Donor, meta$Day) == 1), nrow(meta) == 28,
          length(unique(meta$Donor)) == 14)
counts <- as.matrix(b[, -1])
storage.mode(counts) <- "double"
ensembl <- as.character(b[[1]])
stopifnot(!anyDuplicated(ensembl), all(is.finite(counts)), all(counts >= 0),
          all(counts == round(counts)))
db <- org.Hs.eg.db::org.Hs.eg.db
known <- intersect(ensembl, AnnotationDbi::keys(db, keytype = "ENSEMBL"))
symbols <- AnnotationDbi::mapIds(db, keys = known, keytype = "ENSEMBL",
                               column = "SYMBOL", multiVals = "asNA")
symbol <- unname(symbols[ensembl])
write.csv(data.frame(ENSEMBL = ensembl, Gene = symbol),
          file.path(external_dir, "GSE179248_ensembl_annotation.csv"), row.names = FALSE)
valid <- !is.na(symbol) & nzchar(symbol)
counts <- rowsum(counts[valid, , drop = FALSE], symbol[valid], reorder = TRUE)
meta$Day <- factor(meta$Day, levels = c("D0", "D6"))
design <- model.matrix(~factor(Donor) + Day, meta)
stopifnot(qr(design)$rank == ncol(design))
y <- edgeR::DGEList(counts)
keep <- edgeR::filterByExpr(y, group = meta$Day)
write.csv(data.frame(Gene = rownames(counts), Pass_filterByExpr = keep),
          file.path(external_dir, "GSE179248_count_filter.csv"), row.names = FALSE)
y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
v <- limma::voom(y, design, plot = FALSE)
fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
stats <- limma_result(fit, "DayD6", "GSE179248", "Culture_D6_vs_D0")
meta$Filtered_Library_Size <- y$samples$lib.size
meta$TMM_Factor <- y$samples$norm.factors
write.csv(meta, file.path(external_dir, "GSE179248_analysis_mapping.csv"), row.names = FALSE)
save_tissue(stats, "Normal_to_ADM", 14, 14, 14, "Paired human culture response, not in-vivo progression")
write.csv(data.frame(GSM = meta$GSM, Expression_Column = meta$Expression_Column,
                    Donor = meta$Donor, Day = meta$Day, STAT1_log2CPM = v$E["STAT1", ]),
          file.path(external_dir, "GSE179248_STAT1_donor_expression.csv"), row.names = FALSE)

# Perfect chip/disease confounding cannot be solved by adding a batch covariate.
z <- read_external_geo("GSE101462", expression = FALSE)
z$meta$State <- characteristic(z$meta, "tissue type:")
z$meta$Preservation <- characteristic(z$meta, "tissue storage:")
z$meta$Chip <- characteristic(z$meta, "sentrix_id:")
keep <- z$meta$State %in% c("PDAC", "pancreatitis") & grepl("FFPE", z$meta$Preservation)
d <- model.matrix(~factor(State) + factor(Chip), z$meta[keep, ])
stopifnot(sum(keep) == 13, qr(d)$rank < ncol(d))
write.csv(z$meta, file.path(external_dir, "GSE101462_excluded_confounded_mapping.csv"), row.names = FALSE)
write.csv(data.frame(Dataset = "GSE101462", FFPE_PDAC = 3, FFPE_Pancreatitis = 10,
                    Design_Columns = ncol(d), Design_Rank = qr(d)$rank,
                    Decision = "No inferential validation: disease and chip perfectly confounded"),
          file.path(external_dir, "GSE101462_eligibility_decision.csv"), row.names = FALSE)

all_stats <- do.call(rbind, tissue_results)
frozen <- do.call(rbind, reference_results)
write.csv(all_stats, file.path(external_dir, "tissue_all_gene_contrasts.csv"), row.names = FALSE)
write.csv(frozen, file.path(external_dir, "tissue_frozen_spatial_family.csv"), row.names = FALSE)
write.csv(do.call(rbind, models), file.path(external_dir, "tissue_models.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(frozen, paste(frozen$Dataset, frozen$Comparison)), function(s) {
  spatial_sig <- s$Spatial_BH < 0.05
  data.frame(Dataset = s$Dataset[1], Comparison = s$Comparison[1], Measurable_Frozen = nrow(s),
             External_BH_All_Significant = sum(s$BH_all_genes < 0.05),
             External_BH_Family_Significant = sum(s$BH_frozen_family < 0.05),
             Same_Direction = sum(s$Same_direction),
             Spatial_Contrast_Significant = sum(spatial_sig),
             Spatial_Significant_Same_Direction = sum(spatial_sig & s$Same_direction),
             Dual_AllGene_Significant_SameDirection = sum(spatial_sig & s$Same_direction & s$BH_all_genes < 0.05),
             Spearman_Effect_Correlation = cor(s$log2FC, s$Spatial_log2FC, method = "spearman"))
}))
write.csv(summary, file.path(external_dir, "tissue_replication_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
print(all_stats[all_stats$Gene == "STAT1", ], row.names = FALSE)
