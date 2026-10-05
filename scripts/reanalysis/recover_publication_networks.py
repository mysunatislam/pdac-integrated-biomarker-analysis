"""Recover saved Cytoscape topology without refreshing STRING or target queries."""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reanalysis/results/figure_revision/networks"
X = "{http://www.cs.rpi.edu/XGMML}"
CY = "{http://www.cytoscape.org}"
SOURCE = "R_proj/Proj_Panc/Results/Discovery+validation1.cys"


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(name, data):
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)


def export_graphml(stem, nodes, edges, directed=False):
    ns = "http://graphml.graphdrawing.org/xmlns"
    ET.register_namespace("", ns)
    tag = lambda s: "{" + ns + "}" + s
    r = ET.Element(tag("graphml"))
    for key in ("label", "Evidence_class", "Node_type"):
        ET.SubElement(r, tag("key"), {"id": key, "for": "node", "attr.name": key, "attr.type": "string"})
    ET.SubElement(r, tag("key"), {"id": "edge_type", "for": "edge", "attr.name": "edge_type", "attr.type": "string"})
    g = ET.SubElement(r, tag("graph"), {"id": stem, "edgedefault": "directed" if directed else "undirected"})
    for row in nodes:
        n = ET.SubElement(g, tag("node"), {"id": row["Gene"]})
        for key, value in (("label", row["Gene"]), ("Evidence_class", row["Evidence_class"]), ("Node_type", row["Node_type"])):
            ET.SubElement(n, tag("data"), {"key": key}).text = value
    for i, row in enumerate(edges):
        e = ET.SubElement(g, tag("edge"), {"id": f"e{i}", "source": row["source"], "target": row["target"]})
        ET.SubElement(e, tag("data"), {"key": "edge_type"}).text = row["edge_type"]
    ET.ElementTree(r).write(OUT / (stem + ".graphml"), encoding="utf-8", xml_declaration=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--archive", type=Path, help="Original R_proj.zip archive")
    inputs.add_argument("--session", type=Path, help="Previously recovered Cytoscape CYS session")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    archive_path = args.archive or Path("F:/Thesis/R_proj.zip")
    if args.session is not None:
        payload = args.session.read_bytes()
        archive_path = None
        retrieval_source = str(args.session)
    elif args.archive is not None or archive_path.is_file():
        with zipfile.ZipFile(archive_path) as archive:
            payload = archive.read(SOURCE)
        retrieval_source = str(archive_path)
    else:
        session_path = OUT / "original_validation_session.cys"
        payload = session_path.read_bytes()
        archive_path = None
        retrieval_source = str(session_path)
    (OUT / "original_validation_session.cys").write_bytes(payload)
    with zipfile.ZipFile(io.BytesIO(payload)) as session:
        name = next(n for n in session.namelist() if "/networks/" in n and n.endswith(".xgmml"))
        root = ET.fromstring(session.read(name))
        network = next(g for g in root.findall(".//" + X + "graph") if g.get("id") == "155")
        node_map = {n.get("id"): n.get("label") for n in network.findall(X + "node")}
        table_name = next(n for n in session.namelist() if "SHARED_ATTRS" in n and "CyEdge" in n)
        records = list(csv.reader(io.StringIO(session.read(table_name).decode())))
        edge_table = {r[0]: dict(zip(records[1], r)) for r in records[5:] if r}
        raw_xml = session.read(name)
    (OUT / "original_validation_network.xgmml").write_bytes(raw_xml)
    history = rows(ROOT / "reanalysis/results/integrated/historical12_current_evidence.csv")
    h12 = {r["Gene"]: r for r in history}
    candidate80 = set((ROOT / "results/tables/spatial/valid_genes.txt").read_text().split())
    assert set(node_map.values()) <= candidate80
    assert set(h12) <= set(node_map.values())
    supported = {r["Gene"] for r in history if all(float(r[c]) < .05 for c in
                 ("Spatial_omnibus_all_target_BH", "Paired_tumor_normal_BH", "Chip_PDAC_CP_BH"))}
    assert len(supported) == 7
    classes = {g: ("Broad association support" if g in supported else
                   "Early spatial association" if g == "STAT1" else
                   "Other historical candidate" if g in h12 else "Network context") for g in node_map.values()}
    nodes = [{"Gene": g, "Evidence_class": classes[g], "Node_type": "Gene", "Historical12": g in h12}
             for g in sorted(node_map.values())]
    edges = []
    for e in network.findall(X + "edge"):
        t = edge_table[e.get("id")]
        edges.append({"source": node_map[e.get("source")], "target": node_map[e.get("target")],
                      "edge_type": "Saved STRING functional association",
                      "combined_score": float(t["combined_score"]),
                      "experimental_score": float(t["experimentally_determined_interaction"]),
                      "database_score": float(t["database_annotated"]),
                      "textmining_score": float(t["automated_textmining"]), "Original_edge_SUID": e.get("id")})
    assert len(nodes) == 56 and len(edges) == 112
    assert len({tuple(sorted((e["source"], e["target"]))) for e in edges}) == len(edges)
    write_csv("ppi_nodes.csv", nodes)
    write_csv("ppi_edges.csv", edges)
    subset_nodes = [n for n in nodes if n["Gene"] in h12]
    subset_edges = [e for e in edges if e["source"] in h12 and e["target"] in h12]
    assert len(subset_edges) == 29
    write_csv("historical12_nodes.csv", subset_nodes)
    write_csv("historical12_edges.csv", subset_edges)
    export_graphml("Recovered_STRING_56", nodes, edges)
    export_graphml("Historical12_induced", subset_nodes, subset_edges)

    annotations = rows(ROOT / "reanalysis/results/integrated/historical_miRNA_target_evidence.csv")
    assert len(annotations) == 24
    assert all(r["Reporter_Assay_Present"] == r["Protein_Evidence_Present"] == "FALSE" for r in annotations)
    mir_edges = [{"source": r["MicroRNA"], "target": r["Gene"], "edge_type": "Cached target annotation"}
                 for r in annotations if int(r["Validated_Record_Count"]) > 0]
    mir_nodes = subset_nodes + [{"Gene": m, "Evidence_class": "Cached miRNA query", "Node_type": "miRNA", "Historical12": False}
                               for m in sorted({r["MicroRNA"] for r in annotations})]
    assert len(mir_edges) == 22
    write_csv("mirna_nodes.csv", mir_nodes)
    write_csv("mirna_edges.csv", mir_edges)
    export_graphml("Cached_miRNA_annotations", mir_nodes, mir_edges)
    palette = {"Broad association support": "#0072B2", "Early spatial association": "#E69F00",
               "Other historical candidate": "#D0D4D8", "Network context": "#F3F4F5", "Cached miRNA query": "#009E73"}
    styles = []
    for title, annotation in (("PDAC expression evidence", False), ("PDAC cached annotation", True)):
        styles.append({"title": title, "defaults": [
            {"visualProperty": "NETWORK_BACKGROUND_PAINT", "value": "#FFFFFF"},
            {"visualProperty": "NODE_WIDTH", "value": 75}, {"visualProperty": "NODE_HEIGHT", "value": 30},
            {"visualProperty": "NODE_SHAPE", "value": "ROUND_RECTANGLE"},
            {"visualProperty": "NODE_LABEL_FONT_SIZE", "value": 14},
            {"visualProperty": "NODE_BORDER_WIDTH", "value": 1.2},
            {"visualProperty": "EDGE_WIDTH", "value": 1},
            {"visualProperty": "EDGE_STROKE_UNSELECTED_PAINT", "value": "#A3ABB2"},
            {"visualProperty": "EDGE_TARGET_ARROW_SHAPE", "value": "NONE"},
            {"visualProperty": "EDGE_SOURCE_ARROW_SHAPE", "value": "NONE"},
            {"visualProperty": "EDGE_LINE_TYPE", "value": "LONG_DASH" if annotation else "SOLID"}],
            "mappings": [{"mappingType": "passthrough", "mappingColumn": "label", "mappingColumnType": "String", "visualProperty": "NODE_LABEL"},
                         {"mappingType": "discrete", "mappingColumn": "Evidence_class", "mappingColumnType": "String",
                          "visualProperty": "NODE_FILL_COLOR", "map": [{"key": k, "value": v} for k, v in palette.items()]}]})
    (OUT / "cytoscape_styles.json").write_text(json.dumps(styles, indent=2) + "\n")
    provenance = {"archive": str(archive_path) if archive_path is not None else None,
                  "archive_member": SOURCE, "retrieval_source": retrieval_source,
                  "session_sha256": hashlib.sha256(payload).hexdigest(), "network_id": "155",
                  "nodes": len(nodes), "edges": len(edges), "historical12_edges": len(subset_edges),
                  "minimum_saved_combined_score": min(e["combined_score"] for e in edges),
                  "unrepresented_historical80": sorted(candidate80 - set(node_map.values())),
                  "mirna_annotation_edges": len(mir_edges), "fresh_database_query": False,
                  "rendering": "R/grid vector rendering; Cytoscape desktop is not installed on this host",
                  "edge_interpretation": "STRING functional associations, not necessarily direct physical binding; no causal direction"}
    (OUT / "network_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({k: provenance[k] for k in ("nodes", "edges", "historical12_edges", "minimum_saved_combined_score", "mirna_annotation_edges")}, indent=2))


if __name__ == "__main__":
    main()
