# Publication Audit Inputs

These public inputs support scripts 18-23. The dated, explicitly post hoc
[specification](../../docs/publication_sensitivity_protocol.md) defines the
analyses. Completed results and any preparation limitations are recorded in
the [publication audit](../../docs/publication_audit_methods_results.md).

## Original Records And Signatures

| File | Role |
| --- | --- |
| `E-MTAB-1791_metadata.json` | Original BioStudies study/file-list response |
| `E-MTAB-1791.sdrf.txt` | Original 457-record specimen/assay metadata |
| `MCPcounter_genes.txt` | Authors' gene signatures at pinned commit b6eac73e91c246fcff0bb1a5c68a816cd588fc48 |
| `MCPcounter_signature_license.txt` | Source signature license (CC0) |
| `publication_download_manifest.csv` | Input URLs, sizes and MD5 checksums |
| `GSE208536_published_cases.csv` | Clinical attributes transcribed from source-paper Table 1, with DOI attribution |

The spatial case table corroborates the existing unique clinical-profile
grouping. It is not a directly deposited patient-ID map. The source's incomplete
CASE8 nodal-stage label is retained rather than completed by inference.

## Prepared E-MTAB-1791 Data

The compact matrices, annotations and selected-probe map are written only
after all 254 eligible source arrays are checked. The resumable array manifest
and parser-check record are updated during preparation and may be incomplete
until the final output checks pass.

| File | Role |
| --- | --- |
| `E-MTAB-1791_primary_normalized_probes.csv.gz` | All checked probes and 254 eligible specimen columns |
| `E-MTAB-1791_probe_annotation.csv` | Deposited probe identifiers and annotation |
| `E-MTAB-1791_primary_normalized_genes.csv.gz` | One selected probe per unambiguously mapped current gene symbol |
| `E-MTAB-1791_selected_probe_map.csv` | Gene, selected probe, Entrez ID and deposited symbol |
| `E-MTAB-1791_array_manifest.csv` | Every source array's canonical URL, original ZIP, bytes, MD5, probe count and cache checksum |
| `E-MTAB-1791_parser_equivalence.csv` | Exact fast/base-R parser comparison for an identified source array |

The expression values are the authors' deposited quantile-normalized log2
values, not a new raw-array normalization. Entrez identifiers are mapped with
the pinned `org.Hs.eg.db`; ambiguous/unmapped identifiers are excluded. Multiple
probes use the highest overall mean across all eligible specimens, with probe
ID breaking ties. This rule does not use diagnostic labels. No new batch
correction is applied to the matrices.

The first newly parsed array in each preparation run must match the base-R
reference parser exactly. Every array must pass byte-count, finite-value,
probe-uniqueness, ordering and cross-array annotation checks. Original ZIP
retrieval additionally verifies CRC. Source MD5s identify the exact processed
arrays used; they are locally recorded checksums, not publisher-signed hashes.
Reusing a raw download requires a matching integrity record and checksum.
ZIP members are promoted from temporary files only after their checks pass.

The source arrays repeat large annotation columns. They and the resumable RDS
caches are temporary preparation inputs outside the repository. Approximately
2.5 GB of selected compressed ZIP members is required for a fresh preparation,
or 8.4 GB through the uncompressed fallback. The compact matrices and mappings
allow the association analyses to run without downloading those arrays again.
See the [run instructions](../../scripts/reanalysis/README.md).

## Analysis Units And Limits

The exact eligibility table is
[`E-MTAB-1791_eligibility_mapping.csv`](../../reanalysis/results/publication/E-MTAB-1791_eligibility_mapping.csv).
It retains 195 explicitly PDAC specimens and 59 pancreatitis specimens without
an associated-tumor history. Nine tumor-associated pancreatitis specimens are
excluded. Other tissue contexts are not relabeled as primary controls.

The [source article](https://doi.org/10.1002/ijc.31087) reports 59 CP specimens
from 58 patients. The repeated patient and individual age/sex values cannot be
reconstructed from the SDRF. Source/assay labels are specimen identifiers, not
verified unique donor IDs. Chip adjustment and every possible single-CP
omission are sensitivity analyses; neither replaces missing clinical variables
or a verified donor map. The stored data do not establish clinical specificity,
cellular origin or a causal disease trajectory.

## Attribution

- [BioStudies / E-MTAB-1791](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1791)
- [E-MTAB-1791 source publication](https://doi.org/10.1002/ijc.31087)
- [Spatial source publication](https://doi.org/10.3389/fimmu.2022.961457)
- [MCP-counter method](https://doi.org/10.1186/s13059-016-1070-5)
- [Pinned signature/code source](https://github.com/ebecht/MCPcounter/tree/b6eac73e91c246fcff0bb1a5c68a816cd588fc48)

Attribute the original studies when reusing their data. Repository code
licensing does not replace the original data or third-party signature terms.
