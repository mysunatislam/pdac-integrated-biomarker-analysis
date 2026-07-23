# Integrated PDAC Biomarker Analysis

This repository contains the cleaned effective scripts, processed public input files, derived tables, networks, and manuscript-ready figures for an integrated pancreatic ductal adenocarcinoma (PDAC) biomarker analysis.

The project combines bulk tissue transcriptomics, GeoMx spatial transcriptomics, protein-protein interaction analysis, extracellular-vesicle miRNA discovery, plasma miRNA model development, secondary cross-study assessment, and curated miRNA-hub-gene evidence.

## Analysis Overview

1. Bulk tissue discovery in `GSE143754` identified PDAC, chronic pancreatitis, and normal pancreas differential-expression patterns.
2. A CP-aware malignant-switch filter retained genes altered in PDAC versus CP while excluding CP-versus-normal inflammatory changes.
3. Spatial validation in `GSE208536` tested malignant-switch genes across Normal, ADM, and PDAC GeoMx regions.
4. STRING/Cytoscape prioritization intersected stage-associated genes with PPI hub genes.
5. EV/plasma miRNA analyses used `GSE304572`, `GSE259327`, and `GSE268771` to assess candidate miRNAs and model transferability.
6. Curated multiMiR evidence summarized miR-107 and miR-20a-5p associations with the 12 shared hub genes.

## Repository Layout

```text
scripts/
  01_bulk_tissue_malignant_switch.R
  02_spatial_stage_ppi_hub_analysis.R
  03_ev_mirna_model_transferability.R
  04_curated_multimir_hub_evidence.R

data/source/
  Processed/public GEO source matrices and metadata used by the recovered scripts.

results/tables/
  Final and supporting CSV/TXT result tables grouped by analysis layer.

results/figures/manuscript/
  Exact figures embedded in the corrected manuscript draft, provided as both PNG and PDF.

results/figures/components/
  Effective component plots used to build or support the manuscript figures, provided as PNG/PDF pairs.

results/networks/
  Portable GraphML network exports.

docs/
  Figure captions and source manifest.
```

## Notes on Curation

Raw GEO archives, `.CEL.gz` files, `.RData` workspaces, RStudio metadata, nested `.git` folders, duplicate files, temporary Office locks, and zipped result bundles were excluded.

The included scripts were recovered from the working project and RStudio source history, then path-normalized from the original local `D:/R_proj/Proj_Panc` layout. They are preserved as effective research scripts rather than refactored into a package. Some sections still require public GEO downloads and installed Bioconductor/R packages to rerun fully.

## Key Outputs

- `results/tables/discovery/Malignant_Switch_628_Genes.csv`
- `results/tables/spatial/Corrected_80_Gene_Spatial_Progression_Statistics.csv`
- `results/tables/enrichment_ppi/Supplementary_Table_12_Shared_StageAssociated_PPI_Hub_Genes.csv`
- `results/tables/mirna/GSE259327_Validation_DEG.csv`
- `results/tables/mirna/AUC_CV_results.csv`
- `results/tables/mirna/Supplementary_Validated_Evidence_miR107_miR20a5p_Shared12.csv`
- `results/figures/manuscript/`

## Data Availability

All source datasets are public Gene Expression Omnibus datasets: `GSE143754`, `GSE208536`, `GSE304572`, `GSE259327`, and `GSE268771`.
