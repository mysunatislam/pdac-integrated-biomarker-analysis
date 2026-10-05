# Focused Prior-Art Assessment

Search date: 2026-10-04. This is a focused overlap assessment, not a systematic
review or a guarantee that a claim is novel. The analysis specification was
written after the earlier results and is explicitly post hoc.

## Closest Published Work

| Primary study | Relevant overlap | Implication for this project |
| --- | --- | --- |
| Yang et al., 2022, [spatial m7G/ADM study](https://doi.org/10.3389/fimmu.2022.961457) | Origin of GSE208536. Already examines ADM-related regulators and PDAC using spatial profiling and immune context. | Spatial ADM/PDAC integration itself is not new. Credit the source and distinguish the profile-aware reanalysis, uncertainty audit and cross-cohort tests. |
| [Transcriptional variations in the wider peritumoral tissue environment of pancreatic cancer](https://doi.org/10.1002/ijc.31087) | Origin of E-MTAB-1791. Distinguishes cancer, pancreatitis, normal and peritumoral tissues. | CP-aware tissue comparisons and peritumoral context are established. Reuse requires exact clinical-history eligibility and acknowledgment of the published patient/specimen difference. |
| Zhou et al., 2020, [qualitative transcriptional signature](https://doi.org/10.3389/fmolb.2020.569842) | A 12-gene-pair tissue classifier was trained using five public cohorts, including E-MTAB-1791, E-MEXP-1121 and GSE91035. | A public-data PDAC/pancreatitis tissue classifier is not a first-use concept. Our frozen-family association checks are not a reproduction or external validation of that signature. |
| Wagle et al., 2023, [shared biochemical processes between pancreatitis and PDAC](https://doi.org/10.31557/APJCP.2023.24.5.1601) | Integrates public pancreatic tissue datasets, including E-MTAB-1791 and GSE15471, with DEG, enrichment and network analyses. | Another overlap/hub/pathway list from these datasets is not a sufficient novelty claim. |
| Perez-Diez et al., 2023, [immune and desmoplastic transcriptional signatures](https://doi.org/10.3390/cancers15112887) | Public-cohort meta-analysis of PDAC transcriptional and microenvironment signatures, including datasets used here. | More datasets or a generic immune/stromal interpretation does not establish novelty. The particular falsification and transportability question needs to be stated. |
| Ebelt et al., 2020, [5-azacytidine and PDAC immunity](https://doi.org/10.3389/fimmu.2020.00538) | Mouse ADM-to-PDAC analyses discuss diminished antiviral/immune expression including Stat1. | Do not claim first observation of STAT1-related immune change across ADM/PDAC. Our human spatial late contrast is not significant and cannot establish a transient peak. |
| [IGF inhibition, fibroblasts and STAT1 signaling](https://doi.org/10.3389/fimmu.2024.1382538) | Experimental pancreatic-tumor work reports STAT1 induction in fibroblasts and macrophage/fibroblast immune signaling. | Bulk STAT1 cannot be assumed to originate in malignant epithelial cells. This does not identify the cell of origin in our datasets. |

## Defensible Study Question

**Which tissue-state associations reproduce after accounting for the available
sampling design and clinical context, and do historical/fixed circulating
miRNA panels transport to a prediagnostic setting or add information to a
clinical comparator?**

The contribution is a transparent, metadata-aware reproducibility assessment:
frozen gene/miRNA families, donor/profile-aware contrasts, assay coverage,
composition sensitivity, explicit exclusion of unresolvable batch confounding,
and retained negative validation results. It is not discovery of a new disease
mechanism or a clinically validated screening test.

## Claims To Remove

- A malignant-specific gene identified by subtracting significant DEG lists.
- A confirmed STAT1 ADM peak or evidence of longitudinal progression.
- A tumor-cell-intrinsic STAT1 effect inferred from bulk or mixed spatial ROIs.
- A novel hub-gene panel based on network degree or database records alone.
- Early-detection performance inferred from advanced/on-treatment diagnostic cases.
- Clinical superiority inferred from a nonsignificant paired AUC comparison.
- Cross-analyte model validation when fixed plasma markers are re-estimated in serum.

## Publication Position

Keep the paired spatial and CP-comparator tissue analyses as one biological
evidence stream. Keep plasma/EV/serum results as a distinct transportability
stream. There are no paired tissue/blood observations that validate a regulatory
bridge between them. A joint manuscript needs an explicit methodological
question; two unrelated biomarker narratives should not be joined by an arrow.

Existing studies make a generic discovery manuscript difficult to justify.
Journal suitability must be assessed after the final results and manuscript
are reviewed. No quartile, impact-factor or acceptance prediction is assigned.
