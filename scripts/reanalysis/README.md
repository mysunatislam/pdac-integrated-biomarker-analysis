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

## Additional Cohorts

After the original analysis, run:

```r
source("scripts/reanalysis/run_external.R")
```

This runs steps 12-17 without refitting or tuning the original PLCO models.
It downloads missing public inputs/caches, audits sample groups and chip effects,
fits paired tissue and culture models, tests fixed serum marker families and
nested-CV CA19-9 increments, generates four new PNG/PDF figure pairs and an
external-cohort report, and validates the saved tables. Large SOFT family archives
and the duplicated GSE91035 annotated expression supplement are download caches,
not committed data. Compact platform tables and all effective matrices are kept.
Source URLs/hashes are in `reanalysis/results/external/download_manifest.csv`.

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

## Publication Sensitivities

Steps 18-23 audit the published spatial case profiles, test STAT1 sensitivity
to MCP-counter expression proxies, and analyze the larger E-MTAB-1791
PDAC/pancreatitis comparison. They do not refit or retune PLCO predictions.
Read the dated specification and publication audit before interpreting effects.

```powershell
Rscript scripts/reanalysis/run_publication.R
```

When the prepared matrices in `data/publication/` are present, the workflow
uses them directly. Fresh source preparation is explicitly requested with
`--download`. This is a substantial transfer: approximately 2.5 GB of selected
compressed ZIP members, or 8.4 GB through the uncompressed fallback. The array
files repeat annotation, so only compact matrices, mapping tables and source
checksums are committed. Temporary arrays and resumable RDS caches are stored
outside the repository in LocalAppData on Windows or `~/.cache` on Unix;
`PDAC_DATA_CACHE` overrides this location. Final artifacts stay in the repo.

Install the pinned Python download dependencies with Python 3.11 or later in a separate environment:

```powershell
python -m pip install -r scripts/reanalysis/requirements-publication.txt
Rscript scripts/reanalysis/run_publication.R --download "--python=C:/path/to/python.exe"
```

The downloader uses the exact archives named in the source SDRF and checks
ZIP CRCs, decompressed sizes, MD5s and cross-array annotation consistency.
Completed ZIP members are promoted from temporary files only after integrity
checks; reusing a raw file also requires its matching integrity record and MD5.
Offline interrupted-transfer and CRC tests can be run with
`python -m unittest discover -s scripts/reanalysis/tests`.
The first newly parsed source array in a preparation run is compared exactly
with the base-R reference parser; annotations and numeric vectors must agree.
The saved parser-equivalence record identifies the checked file and its hash.
MCP-counter R code and gene signatures are pinned to a specific Git commit.
The 254 observations are **specimens**, not 254 verified independent patients.
The source reports 59 CP specimens from 58 patients and no individual age/sex
fields in the SDRF. Omit-one-CP and composition analyses are sensitivity checks,
not replacements for missing patient identifiers or clinical adjustment.
