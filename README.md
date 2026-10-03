# Integrated PDAC Biomarker Analysis

This repository contains a public-data reanalysis of pancreatic ductal adenocarcinoma (PDAC) tissue and circulating miRNA studies, alongside the recovered historical project.

**Use the corrected analysis under [reanalysis/](reanalysis/README.md).**
The original scripts and `results/` are preserved for traceability and should
not be used as current evidence for the original manuscript claims.

The reanalysis accounts for repeated spatial ROIs, uses age/sex-adjusted
bulk tissue contrasts, nests plasma feature selection within cross-validation,
and evaluates diagnostic models on the same-assay prediagnostic PLCO cohort.
Every new scientific figure is supplied as both PNG and PDF.

## Current Findings

Spatial tissue-state associations remain, but patient identity is inferred from
deposited clinical profiles and requires confirmation. No bulk PDAC-versus-CP
gene passes FDR < 0.05 after age/sex adjustment in this small cohort. Diagnostic
plasma performance does not transfer convincingly to the full deposited PLCO
sample under the fitted models. These results support an exploratory
reproducibility and transportability study; they do not establish a clinically
validated biomarker panel.

See [methods and results](docs/reanalysis_methods_results.md) and
[claims and limitations](docs/reanalysis_claims_and_limitations.md).

## Repository Layout

```text
scripts/
  01_bulk_tissue_malignant_switch.R
  02_spatial_stage_ppi_hub_analysis.R
  03_ev_mirna_model_transferability.R
  04_curated_multimir_hub_evidence.R
  reanalysis/
    Executable corrected workflow, dependency setup, and output validation.

reanalysis/
  results/       Corrected primary and sensitivity tables.
  figures/       New PNG/PDF scientific figure pairs.
  provenance/    Session information, checksums, and output checks.

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

## Historical Outputs

- `results/tables/discovery/Malignant_Switch_628_Genes.csv`
- `results/tables/spatial/Corrected_80_Gene_Spatial_Progression_Statistics.csv`
- `results/tables/enrichment_ppi/Supplementary_Table_12_Shared_StageAssociated_PPI_Hub_Genes.csv`
- `results/tables/mirna/GSE259327_Validation_DEG.csv`
- `results/tables/mirna/AUC_CV_results.csv`
- `results/tables/mirna/Supplementary_Validated_Evidence_miR107_miR20a5p_Shared12.csv`
- `results/figures/manuscript/`

## Data Availability

All source datasets are public Gene Expression Omnibus datasets: `GSE143754`, `GSE208536`, `GSE304572`, `GSE259327`, and `GSE268771`.
Source studies and public metadata limitations are documented in
[the corrected analysis guide](reanalysis/README.md). Package versions are
recorded in `renv.lock`; run instructions are in
[scripts/reanalysis/README.md](scripts/reanalysis/README.md).
