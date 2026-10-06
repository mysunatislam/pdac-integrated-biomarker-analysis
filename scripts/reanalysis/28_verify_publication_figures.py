"""Verify publication exports and build internal visual-review sheets."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
import pdfplumber
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "reanalysis/figures/publication_revision"
QA = ROOT / "tmp/publication_figure_qa"
POPPLER = Path(r"C:\Users\Kotha\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe")


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def canonical_digest(path):
    data = path.read_bytes()
    if path.suffix.lower() == ".csv":
        data = data.replace(b"\r\n", b"\n")
    return hashlib.md5(data).hexdigest()


def network_imports():
    net = ROOT / "reanalysis/results/figure_revision/networks"
    palette = {"Broad association support": "#0072B2", "Early spatial association": "#E69F00",
               "Other historical candidate": "#D0D4D8", "Network context": "#F3F4F5", "Cached miRNA query": "#009E73"}
    for stem, prefix, position_file in (("Historical12_publication", "historical12", "historical12_publication_positions.csv"),
                                        ("Recovered_STRING_publication", "ppi", "ppi_publication_positions.csv"),
                                        ("Cached_miRNA_publication", "mirna", "mirna_publication_positions.csv")):
        nodes = read(net / (prefix + "_nodes.csv")); edges = read(net / (prefix + "_edges.csv"))
        node_ids = {row["Gene"]: str(i + 1) for i, row in enumerate(nodes)}
        pos = {r["Gene"]: r for r in read(net / position_file)}
        ns = "http://www.cs.rpi.edu/XGMML"; cy = "http://www.cytoscape.org"
        ET.register_namespace("", ns); ET.register_namespace("cy", cy)
        tag = lambda name: "{" + ns + "}" + name
        g = ET.Element(tag("graph"), {"id": "0", "label": stem, "directed": "0"})
        for r in nodes:
            n = ET.SubElement(g, tag("node"), {"id": node_ids[r["Gene"]], "label": r["Gene"]})
            ET.SubElement(n, tag("att"), {"name": "label", "value": r["Gene"], "type": "string"})
            for key in ("Evidence_class", "Node_type"):
                ET.SubElement(n, tag("att"), {"name": key, "value": r[key], "type": "string"})
            # Match the publication layout in importable graph coordinates; no biological scale is implied.
            graphics = ET.SubElement(n, tag("graphics"), {"x": str(float(pos[r["Gene"]]["x"]) * 4),
                "y": str(-float(pos[r["Gene"]]["y"]) * 4), "w": str(116 if r["Node_type"] == "miRNA" else max(56, len(r["Gene"])*6.8+12)),
                "h": "26", "type": "ROUND_RECTANGLE", "fill": palette[r["Evidence_class"]], "outline": "#202428", "width": "1"})
            for key, val in (("NODE_LABEL", r["Gene"].replace("hsa-", "")), ("NODE_LABEL_FONT_SIZE", "14")):
                ET.SubElement(graphics, tag("att"), {"name": key, "value": val, "type": "string"})
        for i, r in enumerate(edges):
            e = ET.SubElement(g, tag("edge"), {"id": str(len(nodes) + i + 1), "source": node_ids[r["source"]],
                                             "target": node_ids[r["target"]], "label": r["edge_type"]})
            for key in ("edge_type", "combined_score", "experimental_score", "database_score", "textmining_score"):
                if key in r:
                    ET.SubElement(e, tag("att"), {"name": key, "value": r[key], "type": "real" if "score" in key else "string"})
            gr = ET.SubElement(e, tag("graphics"), {"fill": "#A3ABB2", "width": "1"})
            for key, val in (("EDGE_TARGET_ARROW_SHAPE", "NONE"), ("EDGE_SOURCE_ARROW_SHAPE", "NONE"),
                             ("EDGE_LINE_TYPE", "LONG_DASH" if prefix == "mirna" else "SOLID")):
                ET.SubElement(gr, tag("att"), {"name": key, "value": val, "type": "string"})
        ET.ElementTree(g).write(net / (stem + ".xgmml"), encoding="utf-8", xml_declaration=True)
        imported = ET.parse(net / (stem + ".xgmml")).getroot()
        imported_nodes = imported.findall(tag("node"))
        labels = {n.get("id"): n.get("label") for n in imported_nodes}
        assert all(n.get("id").isdigit() and any(a.get("name") == "label" and a.get("value") == n.get("label")
                   for a in n.findall(tag("att"))) for n in imported_nodes)
        imported_edges = [(labels[e.get("source")], labels[e.get("target")]) for e in imported.findall(tag("edge"))]
        assert imported_edges == [(e["source"], e["target"]) for e in edges]
    (net / "README.md").write_text("""# Recovered network assets

The original Cytoscape validation session is preserved as `original_validation_session.cys`.
Open the `*_publication.xgmml` networks in Cytoscape to retain the supplied publication positions.
GraphML and separate CSV node/edge tables are also supplied. Import `cytoscape_styles.json`
as visual styles and apply the expression-evidence or cached-annotation style as appropriate.
XGMML includes labels, evidence classes, coordinates, colors, and saved edge scores where available.
The new import files are syntactically checked; opening them in Cytoscape desktop has not been tested
because Cytoscape is not installed on this host. The delivered PDF/SVG/PNG are R/grid renders.

Edges are undirected functional associations or cached annotations, not causal regulatory arrows.
The 24 historical candidates absent from the saved parent network are listed in provenance,
not invented as isolated nodes. Source hashes and recovery details are in network_provenance.json.
""", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checks-only", action="store_true", help="Skip Poppler rendering and visual-review sheets")
    args = parser.parse_args()
    poppler = os.environ.get("PDAC_PDFTOPPM") or (str(POPPLER) if POPPLER.is_file() else shutil.which("pdftoppm"))
    if not args.checks_only and not poppler:
        parser.error("Install Poppler, set PDAC_PDFTOPPM, or use --checks-only for export checks without rendering")
    QA.mkdir(parents=True, exist_ok=True)
    designs = read(ROOT / "reanalysis/provenance/publication_figure_design_checks.csv")
    assert len(designs) == 22
    inputs = json.loads((ROOT / "reanalysis/provenance/figure_revision_input_md5.json").read_text())
    raw_matches = {path: hashlib.md5((ROOT / path).read_bytes()).hexdigest() == digest for path, digest in inputs.items()}
    canonical_path = ROOT / "reanalysis/provenance/figure_revision_input_canonical_md5.json"
    if all(raw_matches.values()):
        canonical = {path: canonical_digest(ROOT / path) for path in inputs}
        canonical_path.write_text(json.dumps(canonical, indent=2) + "\n")
    else:
        canonical = json.loads(canonical_path.read_text())
        assert set(canonical) == set(inputs)
        assert all(canonical_digest(ROOT / path) == digest for path, digest in canonical.items()), "Changed input data"
    checks = []
    for d in designs:
        stem = d["Figure"]; path = FIG / (stem + ".pdf")
        with pdfplumber.open(path) as pdf:
            assert len(pdf.pages) == 1
            p = pdf.pages[0]
            assert len(p.images) == 0, (stem, "Rasterized PDF")
            # pdfminer reports advance width as 'size' for 90-degree Cairo text.
            # The rotation matrix retains the actual font scale for those labels.
            font_sizes = [c["size"] if c["upright"] else math.hypot(c["matrix"][0], c["matrix"][1]) for c in p.chars]
            assert font_sizes and min(font_sizes) >= 7.9, (stem, "Small text")
            assert not [c for c in p.chars if c["x0"] < -1 or c["x1"] > p.width + 1 or c["top"] < -1 or c["bottom"] > p.height + 1], (stem, "Clipped text")
            assert abs(p.width / 72 * 25.4 - float(d["Width_mm"])) < .36, (stem, "PDF page rounding")
            text = p.extract_text() or ""
            min_font = min(font_sizes)
        with Image.open(FIG / (stem + ".png")) as im:
            expected = (float(d["Width_mm"])/25.4*600, float(d["Height_mm"])/25.4*600)
            assert all(abs(a-b) < 2 for a, b in zip(im.size, expected))
            hist = im.convert("L").histogram(); ink = sum(hist[:245]) / (im.width * im.height)
            assert ink > .007
        ET.parse(FIG / (stem + ".svg"))
        if not args.checks_only:
            subprocess.run([poppler, "-singlefile", "-png", "-r", "160", str(path), str(QA / stem)], check=True,
                           capture_output=True, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        checks.append({"Figure": stem, "PDF_Pages": 1, "PNG_dpi": 600, "Min_PDF_font_pt": round(min_font, 2),
                       "Raster_images_in_PDF": 0, "Editable_SVG": True, "PNG_Nonwhite_Fraction": round(ink, 4), "Checks": "PASS"})
        if stem == "Figure_02":
            assert "29 associations" in text and "-3" in text
        if stem == "Figure_S13":
            assert "56 nodes, 112 associations" in text
        if stem == "Figure_S14":
            assert "22 links" in text
    dest = ROOT / "reanalysis/provenance/publication_figure_checks.csv"
    with dest.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(checks[0])); w.writeheader(); w.writerows(checks)
    for start in ([] if args.checks_only else range(0, len(checks), 6)):
        sheet = Image.new("RGB", (1500, 1600), "white"); draw = ImageDraw.Draw(sheet)
        for j, c in enumerate(checks[start:start+6]):
            x = (j % 2)*750; y = (j // 2)*530
            draw.text((x+12, y+8), c["Figure"], fill="black")
            with Image.open(QA / (c["Figure"] + ".png")) as im:
                im.thumbnail((725, 495)); sheet.paste(im, (x+(750-im.width)//2, y+30+(495-im.height)//2))
        sheet.save(QA / f"review_sheet_{start//6+1}.png")
    network_imports()
    print(f"PASS: {len(checks)} vector PDF / 600-dpi PNG / SVG sets; {len(inputs)} unchanged input files")
    if not all(raw_matches.values()):
        print("CSV CRLF/LF differences accepted using the separately saved canonical checksums; all other bytes checked")


if __name__ == "__main__":
    main()
