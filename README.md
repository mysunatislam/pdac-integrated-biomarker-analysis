# Integrated PDAC Biomarker Analysis

This repository contains a public-data reanalysis of pancreatic ductal adenocarcinoma (PDAC) tissue and circulating miRNA studies, alongside the recovered historical project.

**Use the corrected analysis under [reanalysis/](reanalysis/README.md).**
The original scripts and `results/` are preserved for traceability and should
not be used as current evidence for the original manuscript claims.

The reanalysis accounts for repeated spatial ROIs, uses age/sex-adjusted
bulk tissue contrasts, nests plasma feature selection within cross-validation,
and evaluates diagnostic models on the same-assay prediagnostic PLCO cohort.
Every new scientific figure is supplied as both PNG and PDF.

**Revised publication artwork:** [22 main and supplementary figure sets](reanalysis/figures/publication_revision/)
are available as 600-dpi PNG, fully vector PDF and editable SVG. The
[figure revision audit](docs/publication_figure_revision.md) documents the
recovered original Cytoscape session, 56-node STRING graph, 12-gene subgraph,
cached miRNA annotation graph and export checks. The images are rendered with
R/grid; the [editable Cytoscape files](reanalysis/results/figure_revision/networks/)
are supplied separately. These presentation changes do not refit the analyses
or establish a clinical biomarker panel.

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

**Additional public cohorts are now analyzed:** paired tumor/normal tissue,
paired human ADM cultures, and two serum studies. These support parts of the
spatial tissue-state pattern but do not establish CP-specific or prediagnostic
biomarkers. Fixed serum miRNA panels do not show a convincing incremental gain
over binary CA19-9, age and sex. A further tissue cohort was excluded from
inferential validation because chip and disease are perfectly confounded.
See [the external-cohort report](docs/external_validation_methods_results.md)
and [dated specification](docs/external_validation_protocol.md).

The [publication audit](docs/publication_audit_methods_results.md) additionally
cross-checks spatial clinical profiles against the published cases and tests
STAT1's dependence on immune/stromal expression proxies. Paired tumor/normal
STAT1 expression is strongly sensitive to fibroblast-score adjustment; a
tumor-cell-intrinsic effect is not established. The larger E-MTAB-1791
inflammatory-cohort analysis includes 195 PDAC and 59 CP specimens,
with missing individual clinical variables and an unidentified repeated CP
patient explicitly recorded. Of 501 measurable frozen spatial genes, 303 pass
the chip-adjusted family FDR threshold. STAT1's chip-adjusted direction is
positive, but its significance is sensitive to the multiplicity family and
single-CP omissions. This is not proof of malignancy specificity.
The [focused prior-art assessment](docs/publication_prior_art.md) defines the
narrower reproducibility question and the claims that should not be retained.

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
  figures/       Scientific PNG/PDF pairs and revised PNG/PDF/SVG panel sets.
  provenance/    Session information, checksums, and output checks.

data/source/
  Processed/public GEO source matrices and metadata used by the recovered scripts.

data/external/
  Additional public matrices, paired ADM counts and compact platform annotations.

data/publication/
  Pinned composition signatures, source/array manifests and compact inflammatory-cohort matrices.

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

The original source datasets are public Gene Expression Omnibus datasets: `GSE143754`, `GSE208536`, `GSE304572`, `GSE259327`, and `GSE268771`.
The extension uses `GSE15471`, `GSE179248`, `GSE91035`, `GSE59856` and `GSE85589`;
`GSE101462` is retained for its eligibility audit only.
The publication sensitivity extension also uses ArrayExpress/BioStudies
`E-MTAB-1791`; eligible records and analysis status are documented separately.
Source studies and public metadata limitations are documented in
[the corrected analysis guide](reanalysis/README.md). Package versions are
recorded in `renv.lock`; run instructions are in
[scripts/reanalysis/README.md](scripts/reanalysis/README.md).
