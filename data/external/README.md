# Additional Public Inputs

Processed GEO matrices, the human ADM count workbook, and compact platform
annotation tables used by the independent-cohort extension. Source URLs and
MD5 hashes are in `reanalysis/results/external/download_manifest.csv`.

`scripts/reanalysis/12_download_external.R` downloads missing public files and
extracts only platform tables from the SOFT family archives. The large SOFT
archives are download caches and are not committed. The GSE91035 annotated
supplement duplicates the matrix expression values; it is also cached rather
than committed. Its expression match and annotation defects are documented in
the derived probe audit. The downloader restores these caches for reproducible
auditing.

GSE101462 is retained for the sample/chip eligibility audit, not fitted as an
independent cancer-versus-pancreatitis validation cohort. All source data retain
their original provenance and remain subject to applicable source-study reuse
terms.
