#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
out_dir <- file.path("reanalysis", "results", "enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bulk <- read.csv("reanalysis/results/bulk/GSE143754_gene_level_contrasts.csv")
spatial <- read.csv("reanalysis/results/spatial/GSE208536_patient_blocked_all_targets.csv")
selected <- read.csv("reanalysis/results/integrated_tissue/jointly_measurable_spatial_associated_genes.csv")$Gene
universe <- intersect(bulk$Gene, spatial$Gene)
db <- org.Hs.eg.db::org.Hs.eg.db
universe <- intersect(universe, AnnotationDbi::keys(db, keytype = "SYMBOL"))
annotation <- AnnotationDbi::select(db, keys = universe, keytype = "SYMBOL",
                                   columns = c("GOALL", "ONTOLOGYALL"))
annotation <- unique(annotation[!is.na(annotation$GOALL) &
                                  annotation$ONTOLOGYALL == "BP", c("SYMBOL", "GOALL")])
annotated_universe <- intersect(universe, unique(annotation$SYMBOL))
selected <- intersect(selected, annotated_universe)
gene_sets <- split(annotation$SYMBOL, annotation$GOALL)
gene_sets <- lapply(gene_sets, unique)
gene_sets <- gene_sets[lengths(gene_sets) >= 10 & lengths(gene_sets) <= 500]
stopifnot(length(selected) > 0, length(annotated_universe) > length(selected))
results <- do.call(rbind, lapply(names(gene_sets), function(id) {
  members <- intersect(gene_sets[[id]], annotated_universe)
  hits <- intersect(selected, members)
  data.frame(GO_ID = id, Selected_Hits = length(hits), Selected_Total = length(selected),
             Background_Hits = length(members), Background_Total = length(annotated_universe),
             Enrichment_Ratio = (length(hits) / length(selected)) /
               (length(members) / length(annotated_universe)),
             Hypergeometric_p = phyper(length(hits) - 1, length(members),
                                       length(annotated_universe) - length(members),
                                       length(selected), lower.tail = FALSE),
             Genes = paste(sort(hits), collapse = ";"))
}))
results$BH_FDR <- p.adjust(results$Hypergeometric_p, "BH")
if (requireNamespace("GO.db", quietly = TRUE)) {
  results$Term <- unname(AnnotationDbi::mapIds(
    GO.db::GO.db, keys = results$GO_ID, keytype = "GOID", column = "TERM", multiVals = "first"))
}
results <- results[order(results$BH_FDR, results$Hypergeometric_p), ]
write.csv(results, file.path(out_dir, "GO_BP_jointly_measurable_background.csv"), row.names = FALSE)
write.csv(data.frame(Gene = annotated_universe, Selected = annotated_universe %in% selected),
          file.path(out_dir, "GO_BP_analysis_universe.csv"), row.names = FALSE)
cat("GO-annotated jointly measurable genes:", length(annotated_universe), "\n")
cat("Selected genes:", length(selected), "Tested GO BP terms:", nrow(results), "\n")
cat("GO BP terms at FDR < 0.05:", sum(results$BH_FDR < 0.05), "\n")
