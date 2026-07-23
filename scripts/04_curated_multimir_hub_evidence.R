# Recovered effective analysis script
# Source: R_proj\Proj_Panc\.Rproj.user\A1A16783\sources\session-f90c0a1b\9C7297A0-contents
# Paths were normalized from the original local D:/R_proj layout to project-relative paths.
# The script is preserved from the working project history and may require the listed R packages
# plus public GEO downloads to rerun fully from scratch.

# ==========================================================
# Validated evidence analysis:
# miR-107 and miR-20a-5p versus 12 shared PDAC hub genes
# ==========================================================

# ----------------------------------------------------------
# 0. Install/load packages
# ----------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("multiMiR", quietly = TRUE)) {
  BiocManager::install("multiMiR")
}

library(multiMiR)

# ----------------------------------------------------------
# 1. Input settings
# ----------------------------------------------------------

results_dir <- "results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 12 shared stage-associated PPI hub genes
shared12_genes <- c(
  "ACTB",
  "ANLN",
  "CASP3",
  "CCNA2",
  "CCNB1",
  "CENPF",
  "FOXM1",
  "ITGA2",
  "MKI67",
  "PCNA",
  "STAT1",
  "TOP2A"
)

# Final-panel miRNAs that were associated with all 12 genes
test_mirnas <- c(
  "hsa-miR-107",
  "hsa-miR-20a-5p"
)

# ----------------------------------------------------------
# 2. Helper functions
# ----------------------------------------------------------

# For case-insensitive miRNA matching
mirna_key <- function(x) {
  
  x <- trimws(as.character(x))
  x <- sub("^hsa-", "", x, ignore.case = TRUE)
  tolower(x)
}

# Safely obtain a column as character
get_char_column <- function(df, column_name) {
  
  if (column_name %in% names(df)) {
    return(as.character(df[[column_name]]))
  }
  
  return(rep("", nrow(df)))
}

# Correctly extract every PMID, even if multiple IDs are
# stored together in one field
extract_pmids <- function(x) {
  
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  text_all <- paste(x, collapse = "; ")
  
  pmids <- unlist(
    regmatches(
      text_all,
      gregexpr("[0-9]+", text_all)
    )
  )
  
  pmids <- pmids[nchar(pmids) > 0]
  
  sort(unique(pmids))
}

# ----------------------------------------------------------
# 3. Query experimentally validated miRNA-target records
# ----------------------------------------------------------

mm_query <- tryCatch(
  {
    get_multimir(
      org = "hsa",
      mirna = test_mirnas,
      target = shared12_genes,
      table = "validated",
      summary = FALSE
    )
  },
  error = function(e) {
    stop(
      "multiMiR query failed. Check internet connection, ",
      "package installation, and miRNA/gene names.\n\n",
      e$message
    )
  }
)

validated_raw <- as.data.frame(mm_query@data)

if (nrow(validated_raw) == 0) {
  stop(
    "No validated interactions were returned for these ",
    "miRNA-gene pairs."
  )
}

# Keep only validated records if the type field exists
if ("type" %in% names(validated_raw)) {
  validated_raw <- validated_raw[
    validated_raw$type == "validated",
    ,
    drop = FALSE
  ]
}

# ----------------------------------------------------------
# 4. Standardize columns
# ----------------------------------------------------------

validated_raw$Raw_miRNA <- get_char_column(
  validated_raw,
  "mature_mirna_id"
)

validated_raw$Raw_Gene <- get_char_column(
  validated_raw,
  "target_symbol"
)

validated_raw$Database <- get_char_column(
  validated_raw,
  "database"
)

validated_raw$Experiment <- get_char_column(
  validated_raw,
  "experiment"
)

validated_raw$Support_Type <- get_char_column(
  validated_raw,
  "support_type"
)

validated_raw$PubMed_ID_Raw <- get_char_column(
  validated_raw,
  "pubmed_id"
)

validated_raw$miRNA_Key <- mirna_key(
  validated_raw$Raw_miRNA
)

validated_raw$Gene <- toupper(
  trimws(validated_raw$Raw_Gene)
)

# Retain only exact requested miRNA-gene pairs
validated_raw <- validated_raw[
  validated_raw$miRNA_Key %in% mirna_key(test_mirnas) &
    validated_raw$Gene %in% shared12_genes,
  ,
  drop = FALSE
]

if (nrow(validated_raw) == 0) {
  stop(
    "The query returned records, but none matched the exact ",
    "two-miRNA / twelve-gene comparison after filtering."
  )
}

# Restore preferred miRNA labels
validated_raw$MicroRNA <- ifelse(
  validated_raw$miRNA_Key == mirna_key("hsa-miR-107"),
  "hsa-miR-107",
  "hsa-miR-20a-5p"
)

# ----------------------------------------------------------
# 5. Define evidence categories from experiment descriptions
# ----------------------------------------------------------

evidence_text <- paste(
  validated_raw$Experiment,
  validated_raw$Support_Type,
  sep = " | "
)

validated_raw$Reporter_Assay <- grepl(
  "luciferase|reporter[ -]?assay",
  evidence_text,
  ignore.case = TRUE
)

validated_raw$Protein_Evidence <- grepl(
  "western[ -]?blot|immunoblot|protein expression|elisa",
  evidence_text,
  ignore.case = TRUE
)

validated_raw$Transcript_Evidence <- grepl(
  "q-?rt-?pcr|rt-?qpcr|real[ -]?time pcr|microarray|rna-?seq",
  evidence_text,
  ignore.case = TRUE
)

validated_raw$CLIP_or_HighThroughput <- grepl(
  "h?its-?clip|par-?clip|i?clip|clip-?seq|clash|qclash",
  evidence_text,
  ignore.case = TRUE
)

# Save raw retrieved records
write.csv(
  validated_raw,
  file.path(
    results_dir,
    "miR107_miR20a5p_Shared12_Validated_Raw_Records.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------
# 6. Build complete 24-pair matrix
# ----------------------------------------------------------

pair_grid <- expand.grid(
  MicroRNA = test_mirnas,
  Gene = shared12_genes,
  stringsAsFactors = FALSE
)

pair_key <- paste(
  validated_raw$MicroRNA,
  validated_raw$Gene,
  sep = "___"
)

pair_summary_list <- lapply(
  split(validated_raw, pair_key),
  function(x) {
    
    pmids <- extract_pmids(x$PubMed_ID_Raw)
    
    experiment_values <- sort(unique(
      x$Experiment[x$Experiment != ""]
    ))
    
    database_values <- sort(unique(
      x$Database[x$Database != ""]
    ))
    
    support_values <- sort(unique(
      x$Support_Type[x$Support_Type != ""]
    ))
    
    data.frame(
      MicroRNA = unique(x$MicroRNA)[1],
      Gene = unique(x$Gene)[1],
      
      Validated_Record_Count = nrow(x),
      
      Unique_PubMed_Count = length(pmids),
      
      Supporting_Database_Count = length(
        database_values
      ),
      
      Supporting_Databases = ifelse(
        length(database_values) == 0,
        "None",
        paste(database_values, collapse = "; ")
      ),
      
      Reporter_Assay_Present = any(
        x$Reporter_Assay
      ),
      
      Protein_Evidence_Present = any(
        x$Protein_Evidence
      ),
      
      Transcript_Evidence_Present = any(
        x$Transcript_Evidence
      ),
      
      CLIP_or_HighThroughput_Present = any(
        x$CLIP_or_HighThroughput
      ),
      
      Experiments = ifelse(
        length(experiment_values) == 0,
        "None",
        paste(experiment_values, collapse = " | ")
      ),
      
      Support_Types = ifelse(
        length(support_values) == 0,
        "None",
        paste(support_values, collapse = " | ")
      ),
      
      PubMed_IDs = ifelse(
        length(pmids) == 0,
        "None",
        paste(pmids, collapse = "; ")
      ),
      
      stringsAsFactors = FALSE
    )
  }
)

pair_summary <- do.call(
  rbind,
  pair_summary_list
)

# Merge with all 24 possible pairs
evidence_matrix <- merge(
  pair_grid,
  pair_summary,
  by = c("MicroRNA", "Gene"),
  all.x = TRUE
)

# Fill absent interactions
numeric_cols <- c(
  "Validated_Record_Count",
  "Unique_PubMed_Count",
  "Supporting_Database_Count"
)

for (col in numeric_cols) {
  evidence_matrix[[col]][
    is.na(evidence_matrix[[col]])
  ] <- 0
}

logical_cols <- c(
  "Reporter_Assay_Present",
  "Protein_Evidence_Present",
  "Transcript_Evidence_Present",
  "CLIP_or_HighThroughput_Present"
)

for (col in logical_cols) {
  evidence_matrix[[col]][
    is.na(evidence_matrix[[col]])
  ] <- FALSE
}

text_cols <- c(
  "Supporting_Databases",
  "Experiments",
  "Support_Types",
  "PubMed_IDs"
)

for (col in text_cols) {
  evidence_matrix[[col]][
    is.na(evidence_matrix[[col]])
  ] <- "None"
}

# Evidence classification
evidence_matrix$Evidence_Tier <- ifelse(
  evidence_matrix$Reporter_Assay_Present,
  "Tier A: Reporter-assay evidence",
  ifelse(
    evidence_matrix$Protein_Evidence_Present,
    "Tier B: Protein-level evidence",
    ifelse(
      evidence_matrix$CLIP_or_HighThroughput_Present,
      "Tier C: CLIP/high-throughput evidence",
      ifelse(
        evidence_matrix$Validated_Record_Count > 0,
        "Tier D: Other curated validated evidence",
        "Tier E: No validated record"
      )
    )
  )
)

# Sort strongest pairs first
evidence_matrix <- evidence_matrix[
  order(
    evidence_matrix$MicroRNA,
    -as.integer(evidence_matrix$Reporter_Assay_Present),
    -as.integer(evidence_matrix$Protein_Evidence_Present),
    -as.integer(evidence_matrix$CLIP_or_HighThroughput_Present),
    -evidence_matrix$Unique_PubMed_Count,
    -evidence_matrix$Supporting_Database_Count,
    -evidence_matrix$Validated_Record_Count,
    evidence_matrix$Gene
  ),
]

# ----------------------------------------------------------
# 7. Gene-level ranking across both miRNAs
# ----------------------------------------------------------

gene_evidence_summary <- do.call(
  rbind,
  lapply(shared12_genes, function(gene_name) {
    
    x <- evidence_matrix[
      evidence_matrix$Gene == gene_name,
      ,
      drop = FALSE
    ]
    
    x_raw <- validated_raw[
      validated_raw$Gene == gene_name,
      ,
      drop = FALSE
    ]
    
    total_pmids <- extract_pmids(
      x_raw$PubMed_ID_Raw
    )
    
    total_databases <- sort(unique(
      x_raw$Database[x_raw$Database != ""]
    ))
    
    data.frame(
      Gene = gene_name,
      
      miR107_Validated_Records =
        x$Validated_Record_Count[
          x$MicroRNA == "hsa-miR-107"
        ],
      
      miR20a5p_Validated_Records =
        x$Validated_Record_Count[
          x$MicroRNA == "hsa-miR-20a-5p"
        ],
      
      miR107_Unique_PMIDs =
        x$Unique_PubMed_Count[
          x$MicroRNA == "hsa-miR-107"
        ],
      
      miR20a5p_Unique_PMIDs =
        x$Unique_PubMed_Count[
          x$MicroRNA == "hsa-miR-20a-5p"
        ],
      
      miR107_Reporter_Assay =
        x$Reporter_Assay_Present[
          x$MicroRNA == "hsa-miR-107"
        ],
      
      miR20a5p_Reporter_Assay =
        x$Reporter_Assay_Present[
          x$MicroRNA == "hsa-miR-20a-5p"
        ],
      
      miR107_Protein_Evidence =
        x$Protein_Evidence_Present[
          x$MicroRNA == "hsa-miR-107"
        ],
      
      miR20a5p_Protein_Evidence =
        x$Protein_Evidence_Present[
          x$MicroRNA == "hsa-miR-20a-5p"
        ],
      
      miR107_CLIP_or_HighThroughput =
        x$CLIP_or_HighThroughput_Present[
          x$MicroRNA == "hsa-miR-107"
        ],
      
      miR20a5p_CLIP_or_HighThroughput =
        x$CLIP_or_HighThroughput_Present[
          x$MicroRNA == "hsa-miR-20a-5p"
        ],
      
      Total_Validated_Records = nrow(x_raw),
      
      Total_Unique_PMIDs = length(total_pmids),
      
      Total_Supporting_Databases = length(
        total_databases
      ),
      
      Supporting_Databases = ifelse(
        length(total_databases) == 0,
        "None",
        paste(total_databases, collapse = "; ")
      ),
      
      Direct_Reporter_Pairs = sum(
        x$Reporter_Assay_Present
      ),
      
      Protein_Supported_Pairs = sum(
        x$Protein_Evidence_Present
      ),
      
      CLIP_or_HighThroughput_Pairs = sum(
        x$CLIP_or_HighThroughput_Present
      ),
      
      stringsAsFactors = FALSE
    )
  })
)

gene_evidence_summary <- gene_evidence_summary[
  order(
    -gene_evidence_summary$Direct_Reporter_Pairs,
    -gene_evidence_summary$Protein_Supported_Pairs,
    -gene_evidence_summary$CLIP_or_HighThroughput_Pairs,
    -gene_evidence_summary$Total_Unique_PMIDs,
    -gene_evidence_summary$Total_Supporting_Databases,
    -gene_evidence_summary$Total_Validated_Records,
    gene_evidence_summary$Gene
  ),
]

# ----------------------------------------------------------
# 8. miRNA-level summary
# ----------------------------------------------------------

mirna_evidence_summary <- do.call(
  rbind,
  lapply(test_mirnas, function(mirna_name) {
    
    x <- evidence_matrix[
      evidence_matrix$MicroRNA == mirna_name,
      ,
      drop = FALSE
    ]
    
    data.frame(
      MicroRNA = mirna_name,
      Shared_Genes_With_Validated_Record = sum(
        x$Validated_Record_Count > 0
      ),
      Shared_Genes_With_Reporter_Assay = sum(
        x$Reporter_Assay_Present
      ),
      Shared_Genes_With_Protein_Evidence = sum(
        x$Protein_Evidence_Present
      ),
      Shared_Genes_With_CLIP_or_HighThroughput = sum(
        x$CLIP_or_HighThroughput_Present
      ),
      stringsAsFactors = FALSE
    )
  })
)

# ----------------------------------------------------------
# 9. Save final outputs
# ----------------------------------------------------------

write.csv(
  evidence_matrix,
  file.path(
    results_dir,
    "Supplementary_Validated_Evidence_miR107_miR20a5p_Shared12.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_evidence_summary,
  file.path(
    results_dir,
    "Shared12_Gene_Ranking_Validated_miR107_miR20a5p.csv"
  ),
  row.names = FALSE
)

write.csv(
  mirna_evidence_summary,
  file.path(
    results_dir,
    "Shared12_miRNA_Validated_Evidence_Summary.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------
# 10. Print key results
# ----------------------------------------------------------

cat("\n============================================\n")
cat("PAIR-LEVEL VALIDATED EVIDENCE MATRIX\n")
cat("============================================\n")
print(evidence_matrix)

cat("\n============================================\n")
cat("GENE-LEVEL EVIDENCE RANKING\n")
cat("============================================\n")
print(gene_evidence_summary)

cat("\n============================================\n")
cat("miRNA-LEVEL EVIDENCE SUMMARY\n")
cat("============================================\n")
print(mirna_evidence_summary)

cat("\nFiles saved in:\n", results_dir, "\n")

# ==========================================================
# Create Cytoscape files for curated miRNA-hub-gene network
# ==========================================================

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

library(igraph)

results_dir <- "results"

# Final six-panel miRNAs with validated interaction evidence
test_mirnas <- c(
  "hsa-miR-107",
  "hsa-miR-20a-5p"
)

# Keep only pairs with curated validated records
cyto_edges <- evidence_matrix[
  evidence_matrix$Validated_Record_Count > 0,
  c(
    "MicroRNA",
    "Gene",
    "Validated_Record_Count",
    "Unique_PubMed_Count",
    "Supporting_Database_Count",
    "Supporting_Databases",
    "CLIP_or_HighThroughput_Present",
    "Evidence_Tier"
  )
]

# Cytoscape source-target columns
cyto_edges$source <- cyto_edges$MicroRNA
cyto_edges$target <- cyto_edges$Gene

cyto_edges$Interaction <- "Curated high-throughput miRNA-target association"

cyto_edges$Edge_Width_Value <- cyto_edges$Validated_Record_Count

cyto_edges$Edge_Label <- paste0(
  cyto_edges$MicroRNA,
  " â†’ ",
  cyto_edges$Gene,
  " | Records: ",
  cyto_edges$Validated_Record_Count
)

# Put source and target first
cyto_edges <- cyto_edges[
  ,
  c(
    "source",
    "target",
    "Interaction",
    "MicroRNA",
    "Gene",
    "Validated_Record_Count",
    "Unique_PubMed_Count",
    "Supporting_Database_Count",
    "Supporting_Databases",
    "CLIP_or_HighThroughput_Present",
    "Evidence_Tier",
    "Edge_Width_Value",
    "Edge_Label"
  )
]

# ----------------------------------------------------------
# Create node table
# ----------------------------------------------------------

gene_nodes <- data.frame(
  name = gene_evidence_summary$Gene,
  Label = gene_evidence_summary$Gene,
  Node_Type = "Shared stage-associated PPI hub gene",
  Total_Validated_Records =
    gene_evidence_summary$Total_Validated_Records,
  Total_Associated_PMIDs =
    gene_evidence_summary$Total_Unique_PMIDs,
  miR107_Records =
    gene_evidence_summary$miR107_Validated_Records,
  miR20a5p_Records =
    gene_evidence_summary$miR20a5p_Validated_Records,
  stringsAsFactors = FALSE
)

mirna_nodes <- data.frame(
  name = test_mirnas,
  Label = test_mirnas,
  Node_Type = "Final six-miRNA diagnostic-panel member",
  Total_Validated_Records = sapply(
    test_mirnas,
    function(x) {
      sum(
        cyto_edges$Validated_Record_Count[
          cyto_edges$MicroRNA == x
        ]
      )
    }
  ),
  Total_Associated_PMIDs = sapply(
    test_mirnas,
    function(x) {
      sum(
        cyto_edges$Unique_PubMed_Count[
          cyto_edges$MicroRNA == x
        ]
      )
    }
  ),
  miR107_Records = NA,
  miR20a5p_Records = NA,
  stringsAsFactors = FALSE
)

cyto_nodes <- rbind(mirna_nodes, gene_nodes)

# Optional annotation for Cytoscape tooltips
cyto_nodes$Evidence_Note <- "Other shared hub gene"

cyto_nodes$Evidence_Note[
  cyto_nodes$name == "ITGA2"
] <- "Broadest curated record support across both miRNAs"

cyto_nodes$Evidence_Note[
  cyto_nodes$name == "CCNB1"
] <- "Highest associated PMID count across both miRNAs"

# ----------------------------------------------------------
# Export CSV files and GraphML network
# ----------------------------------------------------------

write.csv(
  cyto_edges,
  file.path(
    results_dir,
    "Cytoscape_Curated_miRNA_HubGene_Edges.csv"
  ),
  row.names = FALSE
)

write.csv(
  cyto_nodes,
  file.path(
    results_dir,
    "Cytoscape_Curated_miRNA_HubGene_Nodes.csv"
  ),
  row.names = FALSE
)

g_cyto <- graph_from_data_frame(
  d = cyto_edges,
  directed = TRUE,
  vertices = cyto_nodes
)

write_graph(
  g_cyto,
  file.path(
    results_dir,
    "Supplementary_Curated_miRNA_HubGene_Network.graphml"
  ),
  format = "graphml"
)

cat("Nodes:", vcount(g_cyto), "\n")
cat("Edges:", ecount(g_cyto), "\n")
cat(
  "\nSaved GraphML:\n",
  file.path(
    results_dir,
    "Supplementary_Curated_miRNA_HubGene_Network.graphml"
  ),
  "\n"
)

