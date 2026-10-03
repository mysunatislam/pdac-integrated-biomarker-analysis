# Reanalysis scripts

Run from the repository root with R 4.6 and the dependencies in `renv.lock`.
The original scripts outside this directory and `results/` are historical.
The corrected results are under `reanalysis/`.

```r
source("scripts/reanalysis/install_dependencies.R")
source("scripts/reanalysis/run_all.R")
```

For a locked package restore, install `renv` first, then:

```r
dir.create(".R-library", showWarnings = FALSE)
renv::restore(lockfile = "renv.lock", library = ".R-library", prompt = FALSE)
source("scripts/reanalysis/run_all.R")
```

The complete run includes nested cross-validation and spatial sensitivity
models and can take substantial time. It writes tables, paired PNG/PDF
figures, provenance, and output checks. No GitHub credentials are needed to
run the analyses. All input files are already under `data/source/`.

| Script | Purpose |
| --- | --- |
| 00 | Audit GEO sample labels and workbook mappings |
| 01 | Patient-profile-blocked spatial analysis, complete-profile and leave-one-out sensitivity |
| 02 | EV miRNA filtering, differential tests, exploratory candidate ranking |
| 03 | Diagnostic fitting, nested selection, locked PLCO evaluation |
| 04 | Full platform annotation and age/sex-adjusted limma tissue contrasts |
| 05 | Separate descriptive analysis of the secondary EV cohort |
| 06 | Tissue comparison with adjusted and unadjusted evidence kept explicit |
| 07 | GO BP enrichment using jointly measurable, annotated genes as background |
| 08 | Cohort QC and PLCO calibration plots |
| 09 | Session information, checksums, and package lockfile |
| 10 | Structural and numerical output checks |
| 11 | Methods/results report generated from the final tables |

To regenerate plasma figures from saved scores without refitting models:

```powershell
Rscript scripts/reanalysis/03_plasma_plco_validation.R --figures-only
```

Optional export QA requires Python with Pillow and pypdf, plus Poppler:

```powershell
python scripts/reanalysis/verify_figure_exports.py --poppler pdftoppm
```

This writes `reanalysis/provenance/figure_checks.csv` and PDF renders under
the ignored `tmp/figure_qa/` directory. It checks page counts, paired export
names, dimensions, and nonblank pixels; visual inspection is still required.
