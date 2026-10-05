# Recovered network assets

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
