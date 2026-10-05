"""Package verified publication panels, saved networks and their provenance."""
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reanalysis/figures/publication_revision"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    with (ROOT / "reanalysis/provenance/publication_figure_checks.csv").open(newline="", encoding="utf-8-sig") as stream:
        checks = list(csv.DictReader(stream))
    assert len(checks) == 22 and all(row["Checks"] == "PASS" for row in checks)
    stems = ["Figure_" + str(i).zfill(2) for i in range(1, 8)]
    stems += ["Figure_S" + str(i).zfill(2) for i in range(1, 16)]
    figures = [OUT / (stem + suffix) for stem in stems for suffix in (".png", ".pdf", ".svg")]
    assert all(path.is_file() for path in figures)
    support = list((ROOT / "reanalysis/results/figure_revision/networks").glob("*"))
    support += [ROOT / "docs/publication_figure_revision.md"]
    support += [ROOT / "scripts/reanalysis" / name for name in
                ("27_publication_figure_revision.R", "28_verify_publication_figures.py",
                 "29_package_publication_figures.py", "recover_publication_networks.py")]
    support += [ROOT / "reanalysis/provenance" / name for name in
                ("publication_figure_checks.csv", "publication_figure_design_checks.csv", "figure_revision_input_md5.json",
                 "figure_revision_input_canonical_md5.json")]
    files = sorted(figures + [path for path in support if path.is_file()])
    manifest = [{"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)}
                for path in files]
    archive = OUT / "PDAC_publication_figures.zip"
    with ZipFile(archive, "w", ZIP_DEFLATED, compresslevel=9) as zipped:
        for path in files:
            zipped.write(path, path.relative_to(ROOT).as_posix())
        zipped.writestr("figure_package_manifest.json", json.dumps(manifest, indent=2) + "\n")
    with ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        for row in manifest:
            assert hashlib.sha256(zipped.read(row["path"])).hexdigest() == row["sha256"]
    report = {"archive": archive.relative_to(ROOT).as_posix(), "figure_sets": 22,
              "formats": ["600-dpi PNG", "vector PDF", "editable SVG"],
              "members": len(manifest) + 1, "bytes": archive.stat().st_size,
              "sha256": sha(archive), "crc_and_member_hashes": "PASS",
              "Cytoscape_desktop_import_tested": False}
    (ROOT / "reanalysis/provenance/publication_figure_package.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
