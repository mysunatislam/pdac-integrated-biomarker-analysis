# Publication figure revision

All 22 figures used by the combined manuscript have been rebuilt from saved
results: seven main and fifteen supplementary figures. This revision changes
presentation and restores traceable network displays; it does not refit models
or change the conclusions. Earlier figures and the newer experimental figure
scripts remain preserved.

## Figure files

Each panel set is supplied as a 600-dpi PNG, fully vector PDF and editable SVG
under `reanalysis/figures/publication_revision/`. Arial text is generally 9 pt;
the smallest measured text is 8.39 pt at the supplied publication dimensions.
Color, contrast, axis conventions and uncertainty displays are consistent.
Figure 2 combines recovered topology with the corrected candidate evidence.
Figures S13-S15 show the full recovered STRING graph, the cached target
annotation graph and corrected GO enrichment.

The guide principles are legibility at final size and editable vector artwork,
not imitation of a particular journal's house style. See the
[Nature research figure guide](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/)
and [Cytoscape export documentation](https://manual.cytoscape.org/en/latest/Export_Your_Data.html).
The manuscript remains journal neutral; a selected journal may impose different
panel dimensions, type sizes or upload formats.

## Recovered network evidence

The original `Discovery+validation1.cys` session supplies 56 represented genes
and 112 STRING functional-association edges, with saved combined scores from
0.703 upward. All represented genes belong to the historical 80-gene candidate
set; the other 24 are listed in provenance, not invented as observed isolates.
The 12-gene induced subgraph contains 29 saved edges. No new STRING query,
centrality calculation or MCODE reranking was performed.

Gene colors identify the corrected association evidence, not degree or
validation strength. Node sizes and edge widths are constant. STRING evidence
can include coexpression, text mining and database annotations; these edges
are not asserted to be direct physical binding or causal interactions.

The miRNA graph contains 22 cached annotation links among 24 queried pairs.
Dashed edges have no arrowheads. Reporter and protein flags are absent from
these retrieved summaries, and the tissue and blood cohorts were unpaired.
Record counts and PMID strings are not encoded as independently verified PDAC
experiments. The complete pair matrix retains the two unavailable annotations.

The supplied graphics are rendered with R/grid and igraph, not Cytoscape
desktop, which is not installed on this host. The original CYS session, CSV
node/edge tables, GraphML, visual styles and positioned XGMML import files are
under `reanalysis/results/figure_revision/networks/`. XML syntax and table
consistency are checked; opening the new imports in Cytoscape desktop remains
untested. The original session is preserved byte for byte with its hash.

## Verification

`publication_figure_checks.csv` verifies all 22 PNG/PDF/SVG sets: one-page PDFs,
zero raster images in PDFs, font sizes, no out-of-page text, valid SVG XML,
nonblank PNGs and the intended 600-dpi dimensions. Rotated Cairo labels are
measured from their font transform rather than their glyph advance width.
Vector color legends avoid the raster legends produced by default ggplot2.

The renderer records and rechecks MD5 hashes for 31 input files. The verifier
additionally saves canonical CSV checksums that ignore only CRLF/LF
line-ending differences introduced by Git. Binary inputs remain byte-checked;
data fields, row ordering and all other CSV bytes must match. The volcano
panels retain all 23,207 bulk genes and 440 filtered EV miRNAs, compare saved BH
values against their nominal P values, and explicitly show zero FDR discoveries.
All source analysis files were unchanged. The PDF-rendered figure review sheets
were inspected for every figure, with larger inspection of the networks,
heatmap and ROC panels.

The combined manuscript embeds the revised figures in DOCX, PDF, HTML and
Markdown. Its LaTeX source now includes vector PDF figure assets by default;
keep its figure directory when compiling. The built-in LaTeX compiler reports
`Unable to find standard directories for platform`, so TeX compilation is
unverified. The revised manuscript PDFs are independently typeset from the
same structured sources using ReportLab and pypdf, preserving the vector
figure panels. Word export stalled on this host; DOCX-specific pagination is
unverified, and those files remain local rather than in the downloadable
package. Author approval, missing metadata, design limitations and intended-use
validation remain separate from figure quality. Improved artwork is not an
acceptance guarantee.

## Reproduction

From the project root, with its configured R/Python dependencies:

```powershell
python scripts/reanalysis/recover_publication_networks.py
Rscript scripts/reanalysis/27_publication_figure_revision.R
python scripts/reanalysis/28_verify_publication_figures.py
python scripts/reanalysis/29_package_publication_figures.py
```

Recovery uses the user-supplied archive at `F:/Thesis/R_proj.zip` when available
and otherwise uses the supplied `original_validation_session.cys`. Explicit
inputs can be selected with `--archive` or `--session`; no database refresh is
performed. Subsequent rendering reads the local network files and saved results. The script
also accepts figure stems as arguments to regenerate selected exports. Rebuild
the local combined manuscript with its existing document builder after export
verification. The separate 600-dpi and vector panels are the publication assets;
DOCX embeddings are 300 dpi at their displayed size for practical Word export.

For a host without the bundled Poppler path, the verifier also searches `PATH`
or accepts `PDAC_PDFTOPPM`. `--checks-only` performs the export/data checks
without creating PDF-rendered review sheets; it is not a visual-review claim.
