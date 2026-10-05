# Publication Sensitivity Specification

Dated 2026-10-04, after the previously reported analyses. These are explicitly
post hoc sensitivity and metadata-audit analyses, not prospectively registered
tests. No previous PLCO predictions or thresholds will be retuned.

## Spatial Case Crosswalk

Compare the four deposited clinical attributes with the eight distinct rows in
Table 1 of the source article (doi:10.3389/fimmu.2022.961457). An exact unique
match corroborates a case label but is not a directly deposited ROI-to-patient
identifier. Do not repair unexpected ROI counts by guessing tissue assignments.
Preserve the existing clinical-profile grouping and its numerical results.

## Composition Proxies

Use the authors' MCP-counter R implementation and gene signatures pinned at
ebecht/MCPcounter commit b6eac73e91c246fcff0bb1a5c68a816cd588fc48.
Compute scores on log2 expression separately within each cohort. Save the marker
coverage and the score of every measurable population. A score is eligible for
the prespecified adjustment only if at least half the published markers and at
least three markers are present and it varies across samples. This coverage rule
is an analyst-selected screening rule, not proof that a truncated signature is
validated for a targeted spatial assay. Scores are not percentages or direct
measurements of cells. They cannot compare different cell populations or scales
across platforms, and they do not resolve epithelial/acinar lineage composition.

The named STAT1 sensitivity uses T-cell and fibroblast scores, separately and
together when eligible. These two populations are fixed before inspecting their
adjusted effects. The primary donor/profile-blocked models remain unchanged.
Spatial analysis retains profile and tissue state; GSE15471 retains donor and
tumor/normal state; GSE143754 retains age, sex and CP/PDAC diagnosis. Diagnose
rank deficiency, residual degrees of freedom and collinearity, and report failed
models. STAT1 is not a marker in either adjustment signature. Adjustment is
conditional and may remove biologically mediated effects; it does not identify
a tumor-cell-intrinsic or causal effect. Report all model variants and pointwise
CIs; label their p-values exploratory rather than selecting an adjusted win.

Amendment on 2026-10-04 after inspecting the initial three-cohort composition
results: add the already analyzed paired human ADM culture cohort GSE179248
using the same T-cell/fibroblast coverage rule and model variants. Retain donor
blocking, outcome-independent count filtering, TMM normalization and voom
precision weights. Recompute voom weights for each design. This addition is a
further post hoc analysis, not a test specified before all preceding results.

## E-MTAB-1791 Eligibility And Analysis

Audit all 457 SDRF records before testing effects. Keep exact specimen, assay,
chip and clinical-history labels. Primary cancer samples must be explicitly
pancreatic ductal adenocarcinoma tissue. Primary inflammatory controls must be
pancreatitis tissue from a pancreas with pancreatitis, with no deposited
associated tumor history. Pancreatitis accompanying other tumors and adjacent
non-tumor specimens cannot be relabeled as independent benign controls.
Source names are specimen labels, not proven unique donor identities. The
source paper reports 59 CP specimens but 58 CP patients; the replicated patient
cannot be identified from the SDRF. Thus 254 specimens are not claimed as 254
independent patients. Do not manufacture age/sex values from group summaries.
The deposited 457 records and the source paper's 452 high-quality records differ;
the primary 195 PDAC and 59 CP counts match the paper's tissue counts. The nine
tumor-associated pancreatitis specimens are excluded without inspecting effects.
The secondary all-pancreatitis comparison is deferred because their inclusion
in the published quality-controlled cohort cannot be established. Assess the
unknown CP duplicate using a conservative omit-one-CP sensitivity envelope;
this does not substitute for a verified donor map.

Use the deposited quantile-normalized log2 column and its probe identifiers.
Resolve unambiguous Entrez IDs to current symbols using org.Hs.eg.db and choose
duplicate probes by the highest outcome-blind overall mean. Verify annotation
and probe order across files. Analyze PDAC versus primary pancreatitis. Fit
limma with chip fixed effects only if the
diagnosis coefficient remains estimable with positive residual degrees of
freedom. Report chip imbalance and both unadjusted and estimable chip-adjusted
models. Complete confounding is an exclusion, not something ComBat repairs.

Correct all measurable genes within each contrast by BH. Separately evaluate
the already frozen 506 spatial genes with BH across their measurable subset.
Report STAT1 regardless of direction or significance. No new gene selection,
pooled cross-platform clinical claim, or causal progression claim is allowed.
If files or identifiers cannot be validated, preserve an eligibility audit and
do not describe the cohort as analyzed.

## Sources

- Spatial article: https://doi.org/10.3389/fimmu.2022.961457
- MCP-counter: https://doi.org/10.1186/s13059-016-1070-5
- Source implementation: https://github.com/ebecht/MCPcounter
- E-MTAB-1791: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1791
