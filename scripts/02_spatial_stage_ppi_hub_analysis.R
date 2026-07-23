# Recovered effective analysis script
# Source: R_proj\Proj_Panc\.Rproj.user\A1A16783\sources\session-f90c0a1b\398C0BD5-contents
# Paths were normalized from the original local D:/R_proj layout to project-relative paths.
# The script is preserved from the working project history and may require the listed R packages
# plus public GEO downloads to rerun fully from scratch.

casp3_objects <- find_gene_objects("CASP3")
stat1_objects <- find_gene_objects("STAT1")

casp3_objects
stat1_objects


graphml_files <- list.files(
  ".",
  pattern = "\\.graphml$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

library(igraph)

graph_path <- "./mirnet1.graphml"

g <- read_graph(graph_path, format = "graphml")

vertex_df <- as_data_frame(g, what = "vertices")
edge_df   <- as_data_frame(g, what = "edges")

names(vertex_df)
head(vertex_df)

names(edge_df)
head(edge_df)

core_genes <- c("CASP3", "STAT1")

# Find every GraphML node whose information contains CASP3 or STAT1
core_nodes <- vertex_df$name[
  apply(vertex_df, 1, function(x) {
    any(toupper(trimws(as.character(x))) %in% core_genes)
  })
]

core_nodes

core_edges <- edge_df[
  edge_df$from %in% core_nodes |
    edge_df$to %in% core_nodes,
  ,
  drop = FALSE
]

core_edges

node_lookup <- vertex_df[, c("name", setdiff(names(vertex_df), "name"))]

core_edges_labeled <- merge(
  core_edges,
  node_lookup,
  by.x = "from",
  by.y = "name",
  all.x = TRUE
)

core_edges_labeled

names(vertex_df)
head(vertex_df)
core_nodes
head(core_edges)

graphml_files

hub_genes <- c(
  "FOXM1", "PCNA", "TOP2A", "CCNB1", "MDM2", "STAT1",
  "GSK3B", "CCNA2", "AURKA", "ENO1", "CASP3", "MKI67",
  "ITGA2", "ACTB", "ANLN", "CENPF", "GPI", "IDH1"
)

setdiff(hub_genes, vertex_df$name)

grep("STAT", vertex_df$name, value = TRUE, ignore.case = TRUE)

get_targeting_mirs <- function(gene, edge_df) {
  
  gene <- toupper(gene)
  
  connected <- unique(c(
    edge_df$from[toupper(edge_df$to) == gene],
    edge_df$to[toupper(edge_df$from) == gene]
  ))
  
  connected[
    grepl("^(hsa-)?(miR|let)-", connected, ignore.case = TRUE)
  ]
}

clean_mir <- function(x) {
  x <- trimws(x)
  tolower(sub("^hsa-", "", x, ignore.case = TRUE))
}

# All miRNAs connected to CASP3 in miRNet
casp3_mirs <- get_targeting_mirs("CASP3", edge_df)

# Your 17 EV-associated candidates
ev17_list <- c(
  "miR-15a-5p", "miR-20a-5p", "miR-21-5p", "miR-98-5p",
  "miR-107", "miR-34a-5p", "miR-182-5p", "miR-424-5p",
  "miR-184", "miR-630", "miR-216b-5p", "miR-509-3-5p",
  "miR-1908-5p", "miR-4756-3p", "miR-1306-5p",
  "miR-215-3p", "miR-134-3p"
)

# Final diagnostic panel
panel6_list <- c(
  "miR-184", "miR-107", "miR-216b-5p",
  "miR-20a-5p", "miR-134-3p", "miR-215-3p"
)

# CASP3-connected miRNAs that are EV candidates
casp3_ev17 <- ev17_list[
  clean_mir(ev17_list) %in% clean_mir(casp3_mirs)
]

# CASP3-connected miRNAs that are in the final 6-miRNA panel
casp3_panel6 <- panel6_list[
  clean_mir(panel6_list) %in% clean_mir(casp3_mirs)
]

casp3_ev17
casp3_panel6


casp3_focus <- data.frame(
  miRNA = casp3_ev17,
  target_gene = "CASP3",
  in_final_6_panel = casp3_ev17 %in% casp3_panel6,
  stringsAsFactors = FALSE
)

casp3_focus$category <- ifelse(
  casp3_focus$in_final_6_panel,
  "Final 6-miRNA panel",
  "EV candidate only"
)

casp3_focus

write.csv(
  casp3_focus,
  "results/CASP3_EV_miRNA_subnetwork.csv",
  row.names = FALSE
)

library(igraph)

casp3_edges <- data.frame(
  from = casp3_focus$miRNA,
  to = "CASP3"
)

g_casp3 <- graph_from_data_frame(casp3_edges, directed = TRUE)

V(g_casp3)$NodeType <- ifelse(
  V(g_casp3)$name == "CASP3",
  "Core progression-hub gene",
  ifelse(
    V(g_casp3)$name %in% casp3_panel6,
    "Final 6-miRNA panel",
    "EV candidate miRNA"
  )
)

write_graph(
  g_casp3,
  "results/CASP3_EV_miRNA_subnetwork.graphml",
  format = "graphml"
)


library(igraph)

graph_path <- "results/mirnet1for2genes.graphml"

g2 <- read_graph(graph_path, format = "graphml")

vertex2 <- as_data_frame(g2, what = "vertices")
edge2   <- as_data_frame(g2, what = "edges")

head(vertex2)
head(edge2)

core_genes <- c("CASP3", "STAT1")

is_mir_from <- grepl(
  "^(hsa-)?(miR|let)-",
  edge2$from,
  ignore.case = TRUE
)

is_mir_to <- grepl(
  "^(hsa-)?(miR|let)-",
  
  
  # Check original edge-column names
  names(edge2)
  head(edge2)
  
  # Detect the source and target columns automatically
  if (all(c("from", "to") %in% names(edge2))) {
    source_col <- "from"
    target_col <- "to"
  } else if (all(c("source", "target") %in% names(edge2))) {
    source_col <- "source"
    target_col <- "target"
  } else {
    stop("Could not find source/target columns. Run names(edge2) and check the output.")
  }
  
  # Create a clean edge table
  edges_clean <- data.frame(
    From = as.character(edge2[[source_col]]),
    To   = as.character(edge2[[target_col]]),
    stringsAsFactors = FALSE
  )
  
  core_genes <- c("CASP3", "STAT1")
  
  is_mir_from <- grepl(
    "^(hsa-)?(miR|let)-",
    edges_clean$From,
    ignore.case = TRUE
  )
  
  is_mir_to <- grepl(
    "^(hsa-)?(miR|let)-",
    edges_clean$To,
    ignore.case = TRUE
  )
  
  is_gene_from <- toupper(edges_clean$From) %in% core_genes
  is_gene_to   <- toupper(edges_clean$To) %in% core_genes
  
  # Keep only miRNA â†” CASP3/STAT1 interactions
  core_mir_edges <- edges_clean[
    (is_mir_from & is_gene_to) |
      (is_mir_to & is_gene_from),
    ,
    drop = FALSE
  ]
  
  # Create the two needed columns
  core_mir_edges$miRNA <- ifelse(
    grepl("^(hsa-)?(miR|let)-", core_mir_edges$From, ignore.case = TRUE),
    core_mir_edges$From,
    core_mir_edges$To
  )
  
  core_mir_edges$Target_Gene <- ifelse(
    toupper(core_mir_edges$From) %in% core_genes,
    core_mir_edges$From,
    core_mir_edges$To
  )
  
  # Keep only required columns
  core_mir_edges <- unique(
    core_mir_edges[, c("miRNA", "Target_Gene"), drop = FALSE]
  )
  
  core_mir_edges
  edge2$to,
  ignore.case = TRUE
)

is_core_from <- toupper(edge2$from) %in% core_genes
is_core_to   <- toupper(edge2$to) %in% core_genes

core_mir_edges <- edge2[
  (is_mir_from & is_core_to) |
    (is_mir_to & is_core_from),
  ,
  drop = FALSE
]

core_mir_edges$miRNA <- ifelse(
  grepl("^(hsa-)?(miR|let)-", core_mir_edges$from, ignore.case = TRUE),
  core_mir_edges$from,
  core_mir_edges$to
)

core_mir_edges$Target_Gene <- ifelse(
  toupper(core_mir_edges$from) %in% core_genes,
  core_mir_edges$from,
  core_mir_edges$to
)

core_mir_edges <- core_mir_edges[
  , c("miRNA", "Target_Gene", "direction")
]

core_mir_edges <- unique(core_mir_edges)

core_mir_edges

ev17_list <- c(
  "miR-15a-5p", "miR-20a-5p", "miR-21-5p", "miR-98-5p",
  "miR-107", "miR-34a-5p", "miR-182-5p", "miR-424-5p",
  "miR-184", "miR-630", "miR-216b-5p", "miR-509-3-5p",
  "miR-1908-5p", "miR-4756-3p", "miR-1306-5p",
  "miR-215-3p", "miR-134-3p"
)

panel6_list <- c(
  "miR-184", "miR-107", "miR-216b-5p",
  "miR-20a-5p", "miR-134-3p", "miR-215-3p"
)

clean_mir <- function(x) {
  tolower(sub("^hsa-", "", trimws(x), ignore.case = TRUE))
}

# Keep only EV-supported miRNAs
core_ev_edges <- core_mir_edges[
  clean_mir(core_mir_edges$miRNA) %in% clean_mir(ev17_list),
  ,
  drop = FALSE
]

# Mark whether each miRNA belongs to your final six-miRNA panel
core_ev_edges$Final_6_Panel <- clean_mir(core_ev_edges$miRNA) %in%
  clean_mir(panel6_list)

core_ev_edges$Category <- ifelse(
  core_ev_edges$Final_6_Panel,
  "Final 6-miRNA panel",
  "EV candidate only"
)

core_ev_edges <- core_ev_edges[
  order(core_ev_edges$Target_Gene, core_ev_edges$Category, core_ev_edges$miRNA),
]

core_ev_edges

table(core_ev_edges$Target_Gene)

unique(core_ev_edges$miRNA)

subset(core_ev_edges, Final_6_Panel == TRUE)


library(igraph)

core_final_edges <- subset(
  core_ev_edges,
  Final_6_Panel == TRUE
)

core_final_edges$miRNA <- sub(
  "^hsa-",
  "",
  core_final_edges$miRNA,
  ignore.case = TRUE
)

final_network_edges <- data.frame(
  from = core_final_edges$miRNA,
  to = core_final_edges$Target_Gene,
  stringsAsFactors = FALSE
)

g_final <- graph_from_data_frame(
  final_network_edges,
  directed = TRUE
)

V(g_final)$NodeType <- ifelse(
  V(g_final)$name %in% c("CASP3", "STAT1"),
  "Core progression-hub gene",
  "Final diagnostic miRNA"
)

write_graph(
  g_final,
  "results/Core_CASP3_STAT1_miRNA_network.graphml",
  format = "graphml"
)


write.csv(
  core_final_edges,
  "results/Core_CASP3_STAT1_Final6_Interactions.csv",
  row.names = FALSE
)


# Core genes
core_genes <- c("CASP3", "STAT1")

# Find objects containing BOTH CASP3 and STAT1
locate_core_expression <- function() {
  
  out <- data.frame(
    Object = character(),
    Class = character(),
    Dimensions = character(),
    Gene_location = character(),
    stringsAsFactors = FALSE
  )
  
  for (nm in ls(envir = .GlobalEnv)) {
    
    x <- get(nm, envir = .GlobalEnv, inherits = FALSE)
    
    if (!(is.data.frame(x) || is.matrix(x))) next
    
    rn <- rownames(x)
    cn <- colnames(x)
    
    if (is.null(rn)) rn <- character(0)
    if (is.null(cn)) cn <- character(0)
    
    if (all(core_genes %in% rn)) {
      out <- rbind(
        out,
        data.frame(
          Object = nm,
          Class = class(x)[1],
          Dimensions = paste(dim(x), collapse = " x "),
          Gene_location = "rownames",
          stringsAsFactors = FALSE
        )
      )
    }
    
    if (all(core_genes %in% cn)) {
      out <- rbind(
        out,
        data.frame(
          Object = nm,
          Class = class(x)[1],
          Dimensions = paste(dim(x), collapse = " x "),
          Gene_location = "columns",
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  out
}

locate_core_expression()

locate_progression_tables <- function() {
  
  out <- data.frame(
    Object = character(),
    Columns = character(),
    stringsAsFactors = FALSE
  )
  
  for (nm in ls(envir = .GlobalEnv)) {
    
    x <- get(nm, envir = .GlobalEnv, inherits = FALSE)
    
    if (!is.data.frame(x)) next
    
    cols <- names(x)
    
    has_anova <- any(grepl("anova", cols, ignore.case = TRUE))
    has_spearman <- any(grepl("spearman|rho|trend", cols, ignore.case = TRUE))
    
    if (has_anova && has_spearman) {
      out <- rbind(
        out,
        data.frame(
          Object = nm,
          Columns = paste(cols, collapse = ", "),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  out
}

locate_progression_tables()



sig13 <- read.csv(
  "./sig_genes_13.csv.xls",
  check.names = FALSE
)

hub_genes <- c(
  "FOXM1", "PCNA", "TOP2A", "CCNB1", "MDM2", "STAT1",
  "GSK3B", "CCNA2", "AURKA", "ENO1", "CASP3", "MKI67",
  "ITGA2", "ACTB", "ANLN", "CENPF", "GPI", "IDH1"
)

core_genes <- intersect(sig13$Gene, hub_genes)

core_genes


core_genes <- c("CASP3", "PCNA")

core_edges_2 <- edge_df[
  toupper(edge_df$from) %in% core_genes |
    toupper(edge_df$to) %in% core_genes,
  ,
  drop = FALSE
]

is_mir_from <- grepl("^(hsa-)?(miR|let)-", core_edges_2$from, ignore.case = TRUE)
is_mir_to   <- grepl("^(hsa-)?(miR|let)-", core_edges_2$to, ignore.case = TRUE)

core_mir_edges_2 <- core_edges_2[
  is_mir_from | is_mir_to,
  ,
  drop = FALSE
]

core_mir_edges_2$miRNA <- ifelse(
  is_mir_from[is_mir_from | is_mir_to],
  core_mir_edges_2$from,
  core_mir_edges_2$to
)

core_mir_edges_2$Target_Gene <- ifelse(
  toupper(core_mir_edges_2$from) %in% core_genes,
  core_mir_edges_2$from,
  core_mir_edges_2$to
)

core_mir_edges_2 <- unique(
  core_mir_edges_2[, c("miRNA", "Target_Gene")]
)

core_mir_edges_2

clean_mir <- function(x) {
  tolower(sub("^hsa-", "", trimws(x), ignore.case = TRUE))
}

core_ev_edges_2 <- core_mir_edges_2[
  clean_mir(core_mir_edges_2$miRNA) %in% clean_mir(ev17_list),
  ,
  drop = FALSE
]

core_ev_edges_2$Final_6_Panel <- clean_mir(core_ev_edges_2$miRNA) %in%
  clean_mir(panel6_list)

core_ev_edges_2

core_ev_edges_2
subset(core_ev_edges_2, Final_6_Panel == TRUE)


library(igraph)

# Keep only interactions involving the final six-miRNA panel
core_final <- subset(core_ev_edges_2, Final_6_Panel == TRUE)

# Remove hsa- prefix for cleaner labels
core_final$miRNA <- sub(
  "^hsa-",
  "",
  core_final$miRNA,
  ignore.case = TRUE
)

core_final <- unique(
  core_final[, c("miRNA", "Target_Gene")]
)

core_final



# Save table
write.csv(
  core_final,
  "results/Core_CASP3_PCNA_Final6_Interactions.csv",
  row.names = FALSE
)

# Create Cytoscape network edges
network_edges <- data.frame(
  from = core_final$miRNA,
  to = core_final$Target_Gene,
  stringsAsFactors = FALSE
)

# Define node types
node_info <- data.frame(
  name = unique(c(network_edges$from, network_edges$to)),
  stringsAsFactors = FALSE
)

node_info$NodeType <- ifelse(
  node_info$name %in% c("CASP3", "PCNA"),
  "Core progression-hub gene",
  "Final diagnostic miRNA"
)

# Build and export GraphML
g_core <- graph_from_data_frame(
  network_edges,
  directed = TRUE,
  vertices = node_info
)

write_graph(
  g_core,
  "results/Core_CASP3_PCNA_miRNA_network.graphml",
  format = "graphml"
)


all_files <- list.files(
  ".",
  recursive = TRUE,
  full.names = TRUE
)

spatial_files <- all_files[
  grepl(
    "208536|spatial|expression|expr|metadata|pheno|processed|normal|adm",
    basename(all_files),
    ignore.case = TRUE
  )
]

spatial_files


library(readxl)

spatial_path <- "data/source/GSE208536_Processed_data.xlsx"

# See all sheets in the Excel file
spatial_sheets <- excel_sheets(spatial_path)
spatial_sheets

# Preview every sheet
for (sheet_name in spatial_sheets) {
  
  cat("\n==============================\n")
  cat("SHEET:", sheet_name, "\n")
  cat("==============================\n")
  
  preview <- read_excel(
    spatial_path,
    sheet = sheet_name,
    n_max = 5
  )
  
  print(names(preview))
  print(preview)
}


spatial_raw <- read_excel(
  spatial_path,
  sheet = 1
)

dim(spatial_raw)

names(spatial_raw)

head(spatial_raw)


grep(
  "condition|group|stage|class|sample|adm|normal|pdac",
  names(spatial_raw),
  ignore.case = TRUE,
  value = TRUE
)
spatial_sheets
dim(spatial_raw)
names(spatial_raw)
grep(
  "condition|group|stage|class|sample|adm|normal|pdac",
  names(spatial_raw),
  ignore.case = TRUE,
  value = TRUE
)



# --------------------------------------------------
# STEP 1: Extract sample metadata from GEO series matrix
# --------------------------------------------------

series_path <- "data/source/GSE208536_series_matrix.txt.gz"

geo_lines <- readLines(gzfile(series_path), warn = FALSE)

# Extract values from one GEO metadata line
extract_values <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  values <- parts[-1]
  gsub('^"|"$', "", values)
}

# Extract a named metadata field safely
get_geo_field <- function(field, n_samples) {
  
  hit <- grep(
    paste0("^!", field, "\\t"),
    geo_lines,
    value = TRUE
  )
  
  if (length(hit) == 0) {
    return(rep(NA_character_, n_samples))
  }
  
  extract_values(hit[1])
}

# Number of GEO samples
gsm_ids <- get_geo_field(
  "Sample_geo_accession",
  n_samples = 48
)

n_samples <- length(gsm_ids)

# Basic sample information
sample_meta <- data.frame(
  GSM = gsm_ids,
  Title = get_geo_field("Sample_title", n_samples),
  Source = get_geo_field("Sample_source_name_ch1", n_samples),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Add all sample-characteristics fields
characteristic_lines <- grep(
  "^!Sample_characteristics_ch1\\t",
  geo_lines,
  value = TRUE
)

for (i in seq_along(characteristic_lines)) {
  
  values <- extract_values(characteristic_lines[i])
  
  label <- sub(":.*$", "", values[1])
  label <- make.names(tolower(label))
  
  column_name <- paste0("Characteristic_", i, "_", label)
  
  sample_meta[[column_name]] <- values
}

# Combine all metadata into one searchable column
sample_meta$All_Metadata <- apply(
  sample_meta[, setdiff(names(sample_meta), "GSM"), drop = FALSE],
  1,
  function(x) paste(x, collapse = " | ")
)

# Automatically identify Normal, ADM, and PDAC if they appear in metadata
infer_condition <- function(x) {
  
  x <- toupper(x)
  
  if (grepl("\\bADM\\b|ACINAR.*DUCTAL", x)) {
    return("ADM")
  }
  
  if (grepl("\\bPDAC\\b|PANCREATIC.*ADENOCARCINOMA", x)) {
    return("PDAC")
  }
  
  if (grepl("\\bNORMAL\\b|ADJACENT.*NORMAL", x)) {
    return("Normal")
  }
  
  return(NA_character_)
}

sample_meta$Condition_auto <- vapply(
  sample_meta$All_Metadata,
  infer_condition,
  character(1)
)

# Save metadata for safety
write.csv(
  sample_meta,
  "results/GSE208536_sample_metadata.csv",
  row.names = FALSE
)

# Check extracted metadata
table(sample_meta$Condition_auto, useNA = "ifany")

sample_meta[, c("GSM", "Title", "Source", "Condition_auto", "All_Metadata")]





# --------------------------------------------------
# STEP 2: Match GEO metadata with processed Excel columns
# --------------------------------------------------

make_tma_key <- function(x) {
  
  x <- toupper(x)
  
  match_result <- regexec(
    "TMA\\s*([0-9]+)\\D+0*([0-9]+)",
    x,
    perl = TRUE
  )
  
  hit <- regmatches(x, match_result)[[1]]
  
  if (length(hit) == 3) {
    return(sprintf("TMA%s_%03d", hit[2], as.integer(hit[3])))
  }
  
  NA_character_
}

# Expression sample names from the Excel file
expr_columns <- names(spatial_raw)[-1]

expr_map <- data.frame(
  Expression_Column = expr_columns,
  TMA_Key = vapply(expr_columns, make_tma_key, character(1)),
  stringsAsFactors = FALSE
)

# Extract TMA identity from GEO metadata text
sample_meta$TMA_Key <- vapply(
  sample_meta$All_Metadata,
  make_tma_key,
  character(1)
)

# Match expression columns to metadata
mapping <- merge(
  expr_map,
  sample_meta,
  by = "TMA_Key",
  all.x = TRUE,
  sort = FALSE
)

# Restore original Excel-column order
mapping <- mapping[
  match(expr_map$TMA_Key, mapping$TMA_Key),
]

write.csv(
  mapping,
  "results/GSE208536_expression_metadata_mapping.csv",
  row.names = FALSE
)

mapping[, c(
  "Expression_Column",
  "TMA_Key",
  "GSM",
  "Title",
  "Condition_auto"
)]

table(mapping$Condition_auto, useNA = "ifany")


table(sample_meta$Condition_auto, useNA = "ifany")

table(mapping$Condition_auto, useNA = "ifany")





# --------------------------------------------------
# Rebuild TMA ID mapping with a more flexible parser
# --------------------------------------------------

make_tma_key_v2 <- function(x) {
  
  x <- toupper(as.character(x))
  
  # Accept formats such as:
  # TMA1 | 001, TMA 1-001, TMA_1_001, TMA1_1
  hit <- regexec(
    "TMA[^0-9]*([0-9]+)[^0-9]+0*([0-9]+)",
    x,
    perl = TRUE
  )
  
  parts <- regmatches(x, hit)[[1]]
  
  if (length(parts) == 3) {
    return(sprintf("TMA%s_%03d", parts[2], as.integer(parts[3])))
  }
  
  NA_character_
}

# Create keys from the 48 Excel expression columns
expr_columns <- names(spatial_raw)[-1]

expr_map2 <- data.frame(
  Expression_Column = expr_columns,
  TMA_Key = vapply(expr_columns, make_tma_key_v2, character(1)),
  stringsAsFactors = FALSE
)

# Create keys from GEO sample metadata
sample_meta$TMA_Key_v2 <- vapply(
  sample_meta$All_Metadata,
  make_tma_key_v2,
  character(1)
)

# Check whether TMA IDs were found
table(is.na(expr_map2$TMA_Key))
table(is.na(sample_meta$TMA_Key_v2))

# Match Excel expression columns with GEO metadata
match_index <- match(expr_map2$TMA_Key, sample_meta$TMA_Key_v2)

mapping2 <- cbind(
  expr_map2,
  sample_meta[match_index, c(
    "GSM", "Title", "Source",
    "Condition_auto", "All_Metadata"
  )]
)

# Check matching result
mapping2[, c(
  "Expression_Column",
  "TMA_Key",
  "GSM",
  "Title",
  "Condition_auto"
)]

table(mapping2$Condition_auto, useNA = "ifany")

# Save the verified mapping
write.csv(
  mapping2,
  "results/GSE208536_verified_expression_metadata_mapping.csv",
  row.names = FALSE
)


# Create the TMA-ID parser first
make_tma_key_v2 <- function(x) {
  
  x <- toupper(as.character(x))
  
  hit <- regexec(
    "TMA[^0-9]*([0-9]+)[^0-9]+0*([0-9]+)",
    x,
    perl = TRUE
  )
  
  parts <- regmatches(x, hit)[[1]]
  
  if (length(parts) == 3) {
    return(sprintf("TMA%s_%03d", parts[2], as.integer(parts[3])))
  }
  
  return(NA_character_)
}

# Create matching IDs for the Excel expression columns
expr_columns <- names(spatial_raw)[-1]

expr_map2 <- data.frame(
  Expression_Column = expr_columns,
  TMA_Key = vapply(expr_columns, make_tma_key_v2, character(1)),
  stringsAsFactors = FALSE
)

# Create matching IDs from GEO metadata
sample_meta$TMA_Key_v2 <- vapply(
  sample_meta$All_Metadata,
  make_tma_key_v2,
  character(1)
)

# Check whether IDs were extracted
table(is.na(expr_map2$TMA_Key))
table(is.na(sample_meta$TMA_Key_v2))

# Match each Excel expression column to its GEO condition
match_index <- match(expr_map2$TMA_Key, sample_meta$TMA_Key_v2)

mapping2 <- cbind(
  expr_map2,
  sample_meta[match_index, c(
    "GSM", "Title", "Source",
    "Condition_auto", "All_Metadata"
  )]
)

# Final check
table(mapping2$Condition_auto, useNA = "ifany")

# View first few rows
head(mapping2[, c(
  "Expression_Column", "TMA_Key",
  "GSM", "Title", "Condition_auto"
)], 12)



table(is.na(expr_map2$TMA_Key))
table(is.na(sample_meta$TMA_Key_v2))
table(mapping2$Condition_auto, useNA = "ifany")


options(width = 250)

meta_check <- data.frame(
  GEO_Order = seq_len(nrow(sample_meta)),
  GSM = sample_meta$GSM,
  Title = sample_meta$Title,
  Source = sample_meta$Source,
  Condition = sample_meta$Condition_auto,
  stringsAsFactors = FALSE
)

print(meta_check, row.names = FALSE)

# Show every metadata field and its unique values
char_cols <- grep("^Characteristic_", names(sample_meta), value = TRUE)

for (col in char_cols) {
  cat("\n---------------------------------\n")
  cat(col, "\n")
  print(unique(sample_meta[[col]]))
}

# Check whether the original GEO file contains any TMA labels at all
tma_lines <- grep("TMA", geo_lines, ignore.case = TRUE, value = TRUE)

length(tma_lines)
head(tma_lines, 20)


head(meta_check, 20)
length(tma_lines)
head(tma_lines, 20)




# --------------------------------------------------
# Correctly map the Excel TMA columns to GEO conditions
# --------------------------------------------------

# Extract tab-separated GEO values
extract_geo_values <- function(line) {
  x <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", x)
}

# Get the TMA identity for each GEO sample
description_line <- grep(
  "^!Sample_description\\t",
  geo_lines,
  value = TRUE
)

sample_meta$TMA_Description <- extract_geo_values(description_line[1])

# Make one consistent TMA key, e.g. TMA1_001
make_tma_key <- function(x) {
  
  x <- toupper(as.character(x))
  x <- gsub("\\s+", "", x)
  
  hit <- regexec(
    "TMA([0-9]+)[^0-9]+0*([0-9]+)",
    x,
    perl = TRUE
  )
  
  parts <- regmatches(x, hit)[[1]]
  
  if (length(parts) == 3) {
    return(sprintf("TMA%s_%03d", parts[2], as.integer(parts[3])))
  }
  
  NA_character_
}

# Keys from the processed Excel expression matrix
expr_columns <- names(spatial_raw)[-1]

expr_map <- data.frame(
  Expression_Column = expr_columns,
  TMA_Key = vapply(expr_columns, make_tma_key, character(1)),
  stringsAsFactors = FALSE
)

# Keys from GEO sample descriptions
sample_meta$TMA_Key <- vapply(
  sample_meta$TMA_Description,
  make_tma_key,
  character(1)
)

# Match each Excel column to its GEO sample and condition
match_index <- match(expr_map$TMA_Key, sample_meta$TMA_Key)

mapping2 <- data.frame(
  Expression_Column = expr_map$Expression_Column,
  TMA_Key = expr_map$TMA_Key,
  GSM = sample_meta$GSM[match_index],
  Title = sample_meta$Title[match_index],
  Condition = sample_meta$Condition_auto[match_index],
  stringsAsFactors = FALSE
)

# Confirm that all 48 samples were mapped
table(mapping2$Condition, useNA = "ifany")

head(mapping2, 12)

# Save verified mapping
write.csv(
  mapping2,
  "results/GSE208536_verified_expression_metadata_mapping.csv",
  row.names = FALSE
)

# --------------------------------------------------
# Create CASP3 and PCNA spatial expression dataset
# --------------------------------------------------

core_genes <- c("CASP3", "PCNA")

# Convert Excel expression values to matrix
spatial_expr <- as.matrix(spatial_raw[, -1])
storage.mode(spatial_expr) <- "numeric"

rownames(spatial_expr) <- spatial_raw$TargetName
colnames(spatial_expr) <- names(spatial_raw)[-1]

# Check genes are present
setdiff(core_genes, rownames(spatial_expr))

# Extract expression values in the verified sample order
core_expression <- t(
  spatial_expr[
    core_genes,
    mapping2$Expression_Column,
    drop = FALSE
  ]
)

core_spatial_df <- data.frame(
  Sample = mapping2$Expression_Column,
  Condition = factor(
    mapping2$Condition,
    levels = c("Normal", "ADM", "PDAC")
  ),
  Stage_Number = as.numeric(
    factor(
      mapping2$Condition,
      levels = c("Normal", "ADM", "PDAC")
    )
  ),
  CASP3 = core_expression[, "CASP3"],
  PCNA = core_expression[, "PCNA"],
  row.names = NULL
)

table(core_spatial_df$Condition)

head(core_spatial_df)


table(mapping2$Condition, useNA = "ifany")
head(core_spatial_df)



# --------------------------------------------------
# Rebuild Excel-column to GEO-condition mapping
# --------------------------------------------------

# Convert any TMA label to one standard form: TMA1_001
make_tma_key_safe <- function(x) {
  
  x <- toupper(trimws(as.character(x)))
  
  # Extract all digit groups, e.g. TMA1 | 001 -> c("1", "001")
  nums <- regmatches(x, gregexpr("[0-9]+", x, perl = TRUE))[[1]]
  
  if (length(nums) < 2) {
    return(NA_character_)
  }
  
  paste0(
    "TMA", as.integer(nums[1]),
    "_", sprintf("%03d", as.integer(nums[2]))
  )
}

# Re-extract the sample-description row from the GEO series matrix
description_line <- grep(
  "^!Sample_description\\t",
  geo_lines,
  value = TRUE
)

if (length(description_line) != 1) {
  stop("Could not find exactly one !Sample_description line.")
}

description_values <- strsplit(
  description_line[1],
  "\t",
  fixed = TRUE
)[[1]][-1]

description_values <- gsub('^"|"$', "", description_values)

if (length(description_values) != nrow(sample_meta)) {
  stop("Number of TMA descriptions does not match the number of GEO samples.")
}

# Add the true TMA descriptions to metadata
sample_meta$TMA_Description <- description_values

# Keys from Excel expression-column names
expr_columns <- names(spatial_raw)[-1]

expr_map3 <- data.frame(
  Expression_Column = expr_columns,
  TMA_Key = vapply(expr_columns, make_tma_key_safe, character(1)),
  stringsAsFactors = FALSE
)

# Keys from GEO sample-description metadata
meta_map3 <- data.frame(
  GSM = sample_meta$GSM,
  Title = sample_meta$Title,
  TMA_Description = sample_meta$TMA_Description,
  TMA_Key = vapply(
    sample_meta$TMA_Description,
    make_tma_key_safe,
    character(1)
  ),
  Condition = as.character(sample_meta$Condition_auto),
  stringsAsFactors = FALSE
)

# Diagnostic checks
cat("\nExcel key examples:\n")
print(head(expr_map3, 10))

cat("\nGEO key examples:\n")
print(head(meta_map3, 10))

cat("\nNumber of common TMA keys:\n")
print(length(intersect(expr_map3$TMA_Key, meta_map3$TMA_Key)))

cat("\nDuplicate Excel keys:\n")
print(anyDuplicated(expr_map3$TMA_Key))

cat("\nDuplicate GEO keys:\n")
print(anyDuplicated(meta_map3$TMA_Key))

# Match Excel columns to GEO sample conditions
match_index <- match(expr_map3$TMA_Key, meta_map3$TMA_Key)

mapping2 <- data.frame(
  Expression_Column = expr_map3$Expression_Column,
  TMA_Key = expr_map3$TMA_Key,
  GSM = meta_map3$GSM[match_index],
  Title = meta_map3$Title[match_index],
  Condition = meta_map3$Condition[match_index],
  stringsAsFactors = FALSE
)

mapping2$Condition <- factor(
  mapping2$Condition,
  levels = c("Normal", "ADM", "PDAC")
)

cat("\nFinal condition count:\n")
print(table(mapping2$Condition, useNA = "ifany"))

# Save verified mapping
write.csv(
  mapping2,
  "results/GSE208536_verified_expression_metadata_mapping.csv",
  row.names = FALSE
)




# --------------------------------------------------
# Correct TMA mapping: select the Sample_description
# row that contains the TMA IDs
# --------------------------------------------------

# Find all Sample_description rows
description_lines <- grep(
  "^!Sample_description\\t",
  geo_lines,
  value = TRUE
)

# Keep only the row containing TMA labels
tma_description_line <- description_lines[
  grepl("TMA", description_lines, ignore.case = TRUE)
]

# Check that exactly one TMA-description row was found
length(tma_description_line)

# Extract the 48 TMA labels
description_values <- strsplit(
  tma_description_line[1],
  "\t",
  fixed = TRUE
)[[1]][-1]

description_values <- gsub('^"|"$', "", description_values)

# Convert all TMA formats into TMA1_001 format
make_tma_key_safe <- function(x) {
  
  x <- toupper(trimws(as.character(x)))
  
  nums <- regmatches(
    x,
    gregexpr("[0-9]+", x, perl = TRUE)
  )[[1]]
  
  if (length(nums) < 2) {
    return(NA_character_)
  }
  
  paste0(
    "TMA", as.integer(nums[1]),
    "_", sprintf("%03d", as.integer(nums[2]))
  )
}

# Add GEO TMA IDs to sample metadata
sample_meta$TMA_Description <- description_values

# Keys from Excel expression columns
expr_columns <- names(spatial_raw)[-1]

expr_map3 <- data.frame(
  Expression_Column = expr_columns,
  TMA_Key = vapply(expr_columns, make_tma_key_safe, character(1)),
  stringsAsFactors = FALSE
)

# Keys from GEO metadata
meta_map3 <- data.frame(
  GSM = sample_meta$GSM,
  Title = sample_meta$Title,
  TMA_Key = vapply(
    sample_meta$TMA_Description,
    make_tma_key_safe,
    character(1)
  ),
  Condition = as.character(sample_meta$Condition_auto),
  stringsAsFactors = FALSE
)

# Match Excel columns to Normal / ADM / PDAC conditions
match_index <- match(expr_map3$TMA_Key, meta_map3$TMA_Key)

mapping2 <- data.frame(
  Expression_Column = expr_map3$Expression_Column,
  TMA_Key = expr_map3$TMA_Key,
  GSM = meta_map3$GSM[match_index],
  Title = meta_map3$Title[match_index],
  Condition = meta_map3$Condition[match_index],
  stringsAsFactors = FALSE
)

mapping2$Condition <- factor(
  mapping2$Condition,
  levels = c("Normal", "ADM", "PDAC")
)

# Verify mapping
length(intersect(expr_map3$TMA_Key, meta_map3$TMA_Key))

table(mapping2$Condition, useNA = "ifany")

head(mapping2, 12)



# ==================================================
# CASP3 and PCNA: verified spatial progression analysis
# ==================================================

core_genes <- c("CASP3", "PCNA")

# Safety checks
stopifnot(all(core_genes %in% spatial_raw$TargetName))
stopifnot(!any(is.na(mapping2$Condition)))
stopifnot(all(table(mapping2$Condition) == 16))

# Expression matrix: genes in rows, 48 ROIs in columns
spatial_expr <- as.matrix(spatial_raw[, -1])
storage.mode(spatial_expr) <- "numeric"

rownames(spatial_expr) <- spatial_raw$TargetName
colnames(spatial_expr) <- names(spatial_raw)[-1]

# Put groups in the same order as the Excel expression columns
condition <- factor(
  mapping2$Condition,
  levels = c("Normal", "ADM", "PDAC")
)

stage_number <- unname(
  c(Normal = 1, ADM = 2, PDAC = 3)[as.character(condition)]
)

# Create focused CASP3-PCNA expression dataset
core_spatial_df <- data.frame(
  Sample = mapping2$Expression_Column,
  Condition = condition,
  Stage_Number = stage_number,
  CASP3 = as.numeric(spatial_expr["CASP3", mapping2$Expression_Column]),
  PCNA  = as.numeric(spatial_expr["PCNA",  mapping2$Expression_Column]),
  check.names = FALSE
)

write.csv(
  core_spatial_df,
  "results/Core_CASP3_PCNA_spatial_expression.csv",
  row.names = FALSE
)

print(table(core_spatial_df$Condition))
head(core_spatial_df)


# ==================================================
# ANOVA + Spearman tests for every measured gene
# ==================================================

run_progression_stats <- function(gene_name) {
  
  gene_df <- data.frame(
    Expression = as.numeric(
      spatial_expr[gene_name, mapping2$Expression_Column]
    ),
    Condition = condition,
    Stage_Number = stage_number
  )
  
  fit <- aov(Expression ~ Condition, data = gene_df)
  anova_p <- summary(fit)[[1]]["Condition", "Pr(>F)"]
  
  spearman_test <- suppressWarnings(
    cor.test(
      gene_df$Expression,
      gene_df$Stage_Number,
      method = "spearman",
      exact = FALSE
    )
  )
  
  group_means <- tapply(
    gene_df$Expression,
    gene_df$Condition,
    mean,
    na.rm = TRUE
  )
  
  group_medians <- tapply(
    gene_df$Expression,
    gene_df$Condition,
    median,
    na.rm = TRUE
  )
  
  data.frame(
    Gene = gene_name,
    Mean_Normal = group_means["Normal"],
    Mean_ADM = group_means["ADM"],
    Mean_PDAC = group_means["PDAC"],
    Median_Normal = group_medians["Normal"],
    Median_ADM = group_medians["ADM"],
    Median_PDAC = group_medians["PDAC"],
    ANOVA_p = anova_p,
    Spearman_rho = unname(spearman_test$estimate),
    Spearman_p = spearman_test$p.value,
    Trend = ifelse(
      unname(spearman_test$estimate) > 0,
      "Increasing",
      "Decreasing"
    ),
    stringsAsFactors = FALSE
  )
}

# Run tests across all genes in the GeoMx expression panel
all_spatial_stats <- do.call(
  rbind,
  lapply(rownames(spatial_expr), run_progression_stats)
)

# Conservative FDR correction across the complete measured panel
all_spatial_stats$ANOVA_BH_FDR_all_genes <- p.adjust(
  all_spatial_stats$ANOVA_p,
  method = "BH"
)

all_spatial_stats$Spearman_BH_FDR_all_genes <- p.adjust(
  all_spatial_stats$Spearman_p,
  method = "BH"
)

# Extract CASP3 and PCNA results
core_stats <- all_spatial_stats[
  all_spatial_stats$Gene %in% core_genes,
]

core_stats <- core_stats[
  match(core_genes, core_stats$Gene),
]

print(core_stats)

write.csv(
  all_spatial_stats,
  "results/GSE208536_all_gene_progression_statistics.csv",
  row.names = FALSE
)

write.csv(
  core_stats,
  "results/Core_CASP3_PCNA_progression_statistics.csv",
  row.names = FALSE
)





# ==================================================
# Tukey post-hoc comparisons
# Normal vs ADM, ADM vs PDAC, Normal vs PDAC
# ==================================================

run_tukey_posthoc <- function(gene_name) {
  
  gene_df <- data.frame(
    Expression = as.numeric(
      spatial_expr[gene_name, mapping2$Expression_Column]
    ),
    Condition = condition
  )
  
  fit <- aov(Expression ~ Condition, data = gene_df)
  
  tukey <- TukeyHSD(fit, "Condition")$Condition
  
  data.frame(
    Gene = gene_name,
    Comparison = rownames(tukey),
    Mean_Difference = tukey[, "diff"],
    Lower_CI = tukey[, "lwr"],
    Upper_CI = tukey[, "upr"],
    Tukey_Adjusted_p = tukey[, "p adj"],
    row.names = NULL,
    check.names = FALSE
  )
}

core_posthoc <- do.call(
  rbind,
  lapply(core_genes, run_tukey_posthoc)
)

print(core_posthoc)

write.csv(
  core_posthoc,
  "results/Core_CASP3_PCNA_posthoc_Tukey.csv",
  row.names = FALSE
)


core_stats
core_posthoc


library(ggplot2)

# Convert CASP3 and PCNA expression to long format
plot_df <- rbind(
  data.frame(
    Sample = core_spatial_df$Sample,
    Condition = core_spatial_df$Condition,
    Gene = "CASP3",
    Expression = core_spatial_df$CASP3
  ),
  data.frame(
    Sample = core_spatial_df$Sample,
    Condition = core_spatial_df$Condition,
    Gene = "PCNA",
    Expression = core_spatial_df$PCNA
  )
)

plot_df$Condition <- factor(
  plot_df$Condition,
  levels = c("Normal", "ADM", "PDAC")
)

# Format p-values for the plot
format_p <- function(p) {
  ifelse(
    p < 0.001,
    "<0.001",
    sprintf("%.3f", p)
  )
}

# Create figure annotation for each gene
annotation_df <- do.call(
  rbind,
  lapply(core_stats$Gene, function(gene_name) {
    
    gene_values <- plot_df$Expression[plot_df$Gene == gene_name]
    y_range <- diff(range(gene_values, na.rm = TRUE))
    
    if (y_range == 0) y_range <- 1
    
    stat_row <- core_stats[core_stats$Gene == gene_name, ]
    
    trend_text <- if (
      stat_row$Spearman_BH_FDR_all_genes < 0.05
    ) {
      stat_row$Trend
    } else {
      "No significant monotonic trend"
    }
    
    data.frame(
      Gene = gene_name,
      x = 2,
      y = max(gene_values, na.rm = TRUE) + 0.30 * y_range,
      label = paste0(
        "ANOVA p = ", format_p(stat_row$ANOVA_p),
        "; FDR = ", format_p(stat_row$ANOVA_BH_FDR_all_genes),
        "\nSpearman rho = ", sprintf("%.3f", stat_row$Spearman_rho),
        "; p = ", format_p(stat_row$Spearman_p),
        "; FDR = ", format_p(stat_row$Spearman_BH_FDR_all_genes),
        "\nTrend: ", trend_text
      ),
      stringsAsFactors = FALSE
    )
  })
)

# Build the figure
p_core_expression <- ggplot(
  plot_df,
  aes(x = Condition, y = Expression, fill = Condition)
) +
  geom_boxplot(
    width = 0.62,
    alpha = 0.75,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.80,
    color = "black"
  ) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    lineheight = 0.95
  ) +
  facet_wrap(
    ~Gene,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#4DAF4A",
      "ADM" = "#FFB000",
      "PDAC" = "#E64B35"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.40))
  ) +
  labs(
    title = "Spatial expression of core candidate genes across pancreatic progression",
    subtitle = "Normal â†’ ADM â†’ PDAC; FDR values are adjusted across all measured GeoMx genes",
    x = NULL,
    y = "Normalized expression"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(8, 15, 8, 8)
  )

print(p_core_expression)

# Save publication-quality outputs
ggsave(
  "results/Figure_Core_CASP3_PCNA_Spatial_Expression.png",
  p_core_expression,
  width = 11,
  height = 6,
  dpi = 600
)

ggsave(
  "results/Figure_Core_CASP3_PCNA_Spatial_Expression.pdf",
  p_core_expression,
  width = 11,
  height = 6
)

if (exists("validated")) {
  print(class(validated))
  print(dim(validated))
  print(names(validated))
  print(head(validated))
} else {
  cat("Object 'validated' was not found in the workspace.\n")
}


# ==================================================
# Find the 80 spatially validated gene list
# ==================================================

extract_gene_vector <- function(x) {
  
  # Simple character vector
  if (is.character(x) && is.null(dim(x))) {
    return(unique(toupper(trimws(x))))
  }
  
  # Data frame: look for a gene-symbol-like column
  if (is.data.frame(x)) {
    
    gene_col <- grep(
      "^gene$|gene_symbol|symbol|targetname|target_name",
      names(x),
      ignore.case = TRUE,
      value = TRUE
    )
    
    if (length(gene_col) >= 1) {
      return(unique(toupper(trimws(as.character(x[[gene_col[1]]])))))
    }
  }
  
  # Matrix: use row names if available
  if (is.matrix(x) && !is.null(rownames(x))) {
    return(unique(toupper(trimws(rownames(x)))))
  }
  
  return(NULL)
}

inventory <- data.frame(
  Object = character(),
  Total_Items = integer(),
  GeoMx_Genes_Matched = integer(),
  Contains_CASP3 = logical(),
  Contains_PCNA = logical(),
  stringsAsFactors = FALSE
)

for (nm in ls(envir = .GlobalEnv)) {
  
  x <- get(nm, envir = .GlobalEnv)
  genes <- extract_gene_vector(x)
  
  if (is.null(genes)) next
  
  matched <- intersect(genes, toupper(rownames(spatial_expr)))
  
  if (length(matched) >= 20) {
    inventory <- rbind(
      inventory,
      data.frame(
        Object = nm,
        Total_Items = length(genes),
        GeoMx_Genes_Matched = length(matched),
        Contains_CASP3 = "CASP3" %in% matched,
        Contains_PCNA = "PCNA" %in% matched,
        stringsAsFactors = FALSE
      )
    )
  }
}

inventory <- inventory[
  order(inventory$GeoMx_Genes_Matched),
]

print(inventory)

write.csv(
  inventory,
  "results/Gene_List_Object_Inventory.csv",
  row.names = FALSE
)


# ==================================================
# Find possible 80-gene / 628-gene objects in workspace
# ==================================================

possible_objects <- grep(
  "overlap|candidate|cand|switch|common|sig|gene|spatial",
  ls(),
  ignore.case = TRUE,
  value = TRUE
)

possible_objects

for (nm in possible_objects) {
  
  x <- get(nm)
  
  cat("\n====================================\n")
  cat("OBJECT:", nm, "\n")
  cat("CLASS:", class(x)[1], "\n")
  
  if (is.data.frame(x) || is.matrix(x)) {
    cat("DIMENSION:", paste(dim(x), collapse = " x "), "\n")
    cat("COLUMNS:", paste(head(colnames(x), 10), collapse = ", "), "\n")
    print(head(x, 5))
    
  } else if (is.character(x)) {
    cat("LENGTH:", length(x), "\n")
    print(head(x, 15))
    
  } else {
    str(x, max.level = 1)
  }
}

# ==================================================
# Find saved files related to the 80/628 gene sets
# ==================================================

all_files <- list.files(
  ".",
  recursive = TRUE,
  full.names = TRUE
)

gene_list_files <- all_files[
  grepl(
    "80|628|switch|overlap|candidate|validated|spatial|progress|gene|deg",
    basename(all_files),
    ignore.case = TRUE
  ) &
    grepl(
      "\\.(csv|txt|xlsx|rds|rda|RData)$",
      all_files,
      ignore.case = TRUE
    )
]

gene_list_files


possible_objects
gene_list_files





# ==================================================
# Corrected 80-gene progression analysis
# Use ONLY valid_genes.txt
# ==================================================

# Read the 80 spatially validated genes
spatial80_genes <- readLines(
  "./valid_genes.txt",
  warn = FALSE
)

spatial80_genes <- toupper(trimws(spatial80_genes))
spatial80_genes <- spatial80_genes[nzchar(spatial80_genes)]

# Remove possible header, if present
spatial80_genes <- spatial80_genes[
  !spatial80_genes %in% c(
    "GENE", "GENES", "SYMBOL",
    "GENE_SYMBOL", "TARGETNAME"
  )
]

spatial80_genes <- unique(spatial80_genes)

# Check list size and matching with GeoMx matrix
length(spatial80_genes)

missing_from_geomx <- setdiff(
  spatial80_genes,
  toupper(rownames(spatial_expr))
)

missing_from_geomx

# Keep only genes present in the GeoMx data
spatial80_genes <- intersect(
  spatial80_genes,
  toupper(rownames(spatial_expr))
)

length(spatial80_genes)

# Extract already recalculated statistics for only these 80 genes
stats80 <- all_spatial_stats[
  all_spatial_stats$Gene %in% spatial80_genes,
]

# Correct FDR within the predefined 80-gene candidate set
stats80$ANOVA_BH_FDR_80 <- p.adjust(
  stats80$ANOVA_p,
  method = "BH"
)

stats80$Spearman_BH_FDR_80 <- p.adjust(
  stats80$Spearman_p,
  method = "BH"
)

# Stage-associated genes:
# significant overall difference among Normal, ADM, PDAC
stats80$Stage_Associated <- stats80$ANOVA_BH_FDR_80 < 0.05

# Significant monotonic trend across Normal -> ADM -> PDAC
stats80$Monotonic_Trend <- stats80$Spearman_BH_FDR_80 < 0.05

# Sort by overall stage-association significance
stats80 <- stats80[order(stats80$ANOVA_BH_FDR_80), ]

# Corrected stage-associated gene list
stage_genes_corrected <- stats80$Gene[stats80$Stage_Associated]

# Your 18 hub genes
hub_genes <- c(
  "FOXM1", "PCNA", "TOP2A", "CCNB1", "MDM2", "STAT1",
  "GSK3B", "CCNA2", "AURKA", "ENO1", "CASP3", "MKI67",
  "ITGA2", "ACTB", "ANLN", "CENPF", "GPI", "IDH1"
)

# Corrected core overlap
core_genes_corrected <- intersect(
  stage_genes_corrected,
  hub_genes
)

# Display key outputs
cat("Number of spatially validated genes used:", length(spatial80_genes), "\n\n")

cat("Number of corrected stage-associated genes:",
    length(stage_genes_corrected), "\n\n")

cat("Corrected core progression-associated hub genes:\n")
print(core_genes_corrected)

# Show full statistics for hub genes within the 80-gene set
hub_stats80 <- stats80[
  stats80$Gene %in% hub_genes,
  c(
    "Gene",
    "Mean_Normal", "Mean_ADM", "Mean_PDAC",
    "ANOVA_p", "ANOVA_BH_FDR_80",
    "Spearman_rho", "Spearman_p", "Spearman_BH_FDR_80",
    "Trend", "Stage_Associated", "Monotonic_Trend"
  )
]

print(hub_stats80)

# Save final corrected results
write.csv(
  stats80,
  "results/Corrected_80_Gene_Spatial_Progression_Statistics.csv",
  row.names = FALSE
)

write.csv(
  hub_stats80,
  "results/Corrected_Hub_Genes_Within_80_Spatial_Genes.csv",
  row.names = FALSE
)

write.table(
  stage_genes_corrected,
  "results/Corrected_Stage_Associated_Genes.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  core_genes_corrected,
  "results/Corrected_Core_Progression_Hub_Genes.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)



length(spatial80_genes)
length(stage_genes_corrected)
core_genes_corrected
hub_stats80



# ==================================================
# Corrected hub-overlap and progression classification
# ==================================================

core_hub_stage_genes <- intersect(
  stage_genes_corrected,
  hub_genes
)

core_progression_hubs <- stats80[
  stats80$Gene %in% core_hub_stage_genes &
    stats80$Stage_Associated &
    stats80$Monotonic_Trend,
  ,
  drop = FALSE
]

stage_only_hubs <- stats80[
  stats80$Gene %in% core_hub_stage_genes &
    stats80$Stage_Associated &
    !stats80$Monotonic_Trend,
  ,
  drop = FALSE
]

core_hub_summary <- stats80[
  stats80$Gene %in% core_hub_stage_genes,
  c(
    "Gene",
    "Mean_Normal", "Mean_ADM", "Mean_PDAC",
    "ANOVA_p", "ANOVA_BH_FDR_80",
    "Spearman_rho", "Spearman_p", "Spearman_BH_FDR_80",
    "Trend", "Stage_Associated", "Monotonic_Trend"
  )
]

core_hub_summary$Classification <- ifelse(
  core_hub_summary$Monotonic_Trend,
  "Core progression-associated hub gene",
  "Stage-associated hub gene"
)

core_hub_summary <- core_hub_summary[
  order(
    core_hub_summary$Classification,
    core_hub_summary$ANOVA_BH_FDR_80
  ),
]

print(core_hub_summary)

cat("\nCore progression-associated hubs:\n")
print(core_progression_hubs$Gene)

cat("\nStage-associated hubs without significant monotonic trend:\n")
print(stage_only_hubs$Gene)

write.csv(
  core_hub_summary,
  "results/Corrected_Core_Hub_Gene_Summary.csv",
  row.names = FALSE
)



# ==================================================
# Focused post-hoc analysis: CASP3 and STAT1
# ==================================================

focused_genes <- c("CASP3", "STAT1")

# Correct FDR statistics within the 80-gene candidate set
focused_stats <- stats80[
  match(focused_genes, stats80$Gene),
  c(
    "Gene",
    "Mean_Normal", "Mean_ADM", "Mean_PDAC",
    "ANOVA_p", "ANOVA_BH_FDR_80",
    "Spearman_rho", "Spearman_p", "Spearman_BH_FDR_80",
    "Trend", "Stage_Associated", "Monotonic_Trend"
  )
]

print(focused_stats)

# Tukey post-hoc comparisons
run_focused_tukey <- function(gene_name) {
  
  gene_df <- data.frame(
    Expression = as.numeric(
      spatial_expr[gene_name, mapping2$Expression_Column]
    ),
    Condition = factor(
      mapping2$Condition,
      levels = c("Normal", "ADM", "PDAC")
    )
  )
  
  fit <- aov(Expression ~ Condition, data = gene_df)
  tukey <- TukeyHSD(fit, "Condition")$Condition
  
  data.frame(
    Gene = gene_name,
    Comparison = rownames(tukey),
    Mean_Difference = tukey[, "diff"],
    Lower_CI = tukey[, "lwr"],
    Upper_CI = tukey[, "upr"],
    Tukey_Adjusted_p = tukey[, "p adj"],
    row.names = NULL
  )
}

focused_posthoc <- do.call(
  rbind,
  lapply(focused_genes, run_focused_tukey)
)

print(focused_posthoc)

write.csv(
  focused_stats,
  "results/Focused_CASP3_STAT1_Statistics.csv",
  row.names = FALSE
)

write.csv(
  focused_posthoc,
  "results/Focused_CASP3_STAT1_Tukey_Posthoc.csv",
  row.names = FALSE
)




# ==================================================
# Focused statistics: CASP3 and STAT1
# ==================================================

focused_genes <- c("CASP3", "STAT1")

# Extract their corrected 80-gene FDR statistics
focused_stats <- stats80[
  match(focused_genes, stats80$Gene),
  c(
    "Gene",
    "Mean_Normal", "Mean_ADM", "Mean_PDAC",
    "ANOVA_p", "ANOVA_BH_FDR_80",
    "Spearman_rho", "Spearman_p", "Spearman_BH_FDR_80",
    "Trend", "Stage_Associated", "Monotonic_Trend"
  )
]

print(focused_stats)

# Tukey post-hoc comparisons
run_focused_tukey <- function(gene_name) {
  
  gene_df <- data.frame(
    Expression = as.numeric(
      spatial_expr[gene_name, mapping2$Expression_Column]
    ),
    Condition = factor(
      mapping2$Condition,
      levels = c("Normal", "ADM", "PDAC")
    )
  )
  
  fit <- aov(Expression ~ Condition, data = gene_df)
  tukey <- TukeyHSD(fit, "Condition")$Condition
  
  data.frame(
    Gene = gene_name,
    Comparison = rownames(tukey),
    Mean_Difference = tukey[, "diff"],
    Lower_CI = tukey[, "lwr"],
    Upper_CI = tukey[, "upr"],
    Tukey_Adjusted_p = tukey[, "p adj"],
    row.names = NULL
  )
}

focused_posthoc <- do.call(
  rbind,
  lapply(focused_genes, run_focused_tukey)
)

print(focused_posthoc)

write.csv(
  focused_stats,
  "results/Focused_CASP3_STAT1_Statistics.csv",
  row.names = FALSE
)

write.csv(
  focused_posthoc,
  "results/Focused_CASP3_STAT1_Tukey_Posthoc.csv",
  row.names = FALSE
)



focused_stats
focused_posthoc



# ==================================================
# Figure: CASP3 and STAT1 across Normal -> ADM -> PDAC
# ==================================================

library(ggplot2)

focused_genes <- c("CASP3", "STAT1")

# Long-format data for plotting
plot_df <- rbind(
  data.frame(
    Sample = core_spatial_df$Sample,
    Condition = factor(
      mapping2$Condition,
      levels = c("Normal", "ADM", "PDAC")
    ),
    Gene = "CASP3",
    Expression = as.numeric(
      spatial_expr["CASP3", mapping2$Expression_Column]
    )
  ),
  data.frame(
    Sample = core_spatial_df$Sample,
    Condition = factor(
      mapping2$Condition,
      levels = c("Normal", "ADM", "PDAC")
    ),
    Gene = "STAT1",
    Expression = as.numeric(
      spatial_expr["STAT1", mapping2$Expression_Column]
    )
  )
)

# Format p-values
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

# Annotation text for each panel
annotation_df <- data.frame(
  Gene = c("CASP3", "STAT1"),
  x = c("ADM", "ADM"),
  y = c(
    max(plot_df$Expression[plot_df$Gene == "CASP3"]) * 1.30,
    max(plot_df$Expression[plot_df$Gene == "STAT1"]) * 1.28
  ),
  label = c(
    paste0(
      "ANOVA p = ", format_p(focused_stats$ANOVA_p[focused_stats$Gene == "CASP3"]),
      "; BH-FDR = ", format_p(focused_stats$ANOVA_BH_FDR_80[focused_stats$Gene == "CASP3"]),
      "\nSpearman rho = ", sprintf("%.3f", focused_stats$Spearman_rho[focused_stats$Gene == "CASP3"]),
      "; BH-FDR = ", format_p(focused_stats$Spearman_BH_FDR_80[focused_stats$Gene == "CASP3"]),
      "\nPattern: PDAC-associated increase"
    ),
    paste0(
      "ANOVA p = ", format_p(focused_stats$ANOVA_p[focused_stats$Gene == "STAT1"]),
      "; BH-FDR = ", format_p(focused_stats$ANOVA_BH_FDR_80[focused_stats$Gene == "STAT1"]),
      "\nSpearman rho = ", sprintf("%.3f", focused_stats$Spearman_rho[focused_stats$Gene == "STAT1"]),
      "; BH-FDR = ", format_p(focused_stats$Spearman_BH_FDR_80[focused_stats$Gene == "STAT1"]),
      "\nPattern: ADM-enriched, non-monotonic"
    )
  ),
  stringsAsFactors = FALSE
)

# Create plot
p_casp3_stat1 <- ggplot(
  plot_df,
  aes(x = Condition, y = Expression, fill = Condition)
) +
  geom_boxplot(
    width = 0.62,
    alpha = 0.75,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.8,
    color = "black"
  ) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    lineheight = 0.95
  ) +
  facet_wrap(
    ~Gene,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#4DAF4A",
      "ADM" = "#FFB000",
      "PDAC" = "#E64B35"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.42))
  ) +
  labs(
    title = "Spatial expression of selected core hub genes across pancreatic progression",
    subtitle = "Normal â†’ ADM â†’ PDAC; FDR adjusted across 80 spatially validated genes",
    x = NULL,
    y = "Normalized GeoMx expression"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(8, 15, 8, 8)
  )

print(p_casp3_stat1)

# Save figure
ggsave(
  "results/Figure_CASP3_STAT1_Spatial_Expression.png",
  p_casp3_stat1,
  width = 11,
  height = 6,
  dpi = 600
)

ggsave(
  "results/Figure_CASP3_STAT1_Spatial_Expression.pdf",
  p_casp3_stat1,
  width = 11,
  height = 6
)

# ==================================================
# 12 shared stage-associated PPI hub genes:
# table + overlap figure + spatial heatmap
# ==================================================

library(ggplot2)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Check needed objects
needed_objects <- c("core_hub_summary", "spatial_expr", "mapping2")

missing_objects <- needed_objects[
  !vapply(needed_objects, exists, logical(1), envir = .GlobalEnv)
]

if (length(missing_objects) > 0) {
  stop(
    paste(
      "Missing object(s):",
      paste(missing_objects, collapse = ", ")
    )
  )
}

# --------------------------------------------------
# 1. Create final 12-gene shared-result table
# --------------------------------------------------

shared12_table <- core_hub_summary

shared12_table <- shared12_table[
  order(shared12_table$ANOVA_BH_FDR_80),
]

rownames(shared12_table) <- NULL

# Rename this carefully:
# Significant Spearman FDR = ordinal stage association,
# not necessarily perfectly monotonic mean expression.
shared12_table$Ordinal_Association_FDR_Significant <- 
  shared12_table$Spearman_BH_FDR_80 < 0.05

mean_columns <- c("Mean_Normal", "Mean_ADM", "Mean_PDAC")

shared12_table$Peak_Expression_Group <- c(
  "Normal", "ADM", "PDAC"
)[
  max.col(
    as.matrix(shared12_table[, mean_columns]),
    ties.method = "first"
  )
]

shared12_table$Interpretation <- ifelse(
  shared12_table$Ordinal_Association_FDR_Significant &
    shared12_table$Peak_Expression_Group == "PDAC",
  
  "PDAC-enriched with significant positive ordinal stage association",
  
  ifelse(
    shared12_table$Ordinal_Association_FDR_Significant &
      shared12_table$Peak_Expression_Group == "ADM",
    
    "ADM-enriched with significant positive ordinal stage association",
    
    "Stage-associated without FDR-significant ordinal stage association"
  )
)

print(shared12_table)

write.csv(
  shared12_table,
  file.path(
    results_dir,
    "Table_12_Shared_Stage_Associated_PPI_Hub_Genes.csv"
  ),
  row.names = FALSE
)

# --------------------------------------------------
# 2. Create overlap summary figure
# 59 spatial stage-associated genes âˆ© 18 hubs = 12 genes
# --------------------------------------------------

wrap_gene_names <- function(x, genes_per_line = 4) {
  
  groups <- split(
    x,
    ceiling(seq_along(x) / genes_per_line)
  )
  
  paste(
    vapply(
      groups,
      function(z) paste(z, collapse = ", "),
      character(1)
    ),
    collapse = "\n"
  )
}

shared_gene_text <- wrap_gene_names(
  shared12_table$Gene,
  genes_per_line = 4
)

p_overlap <- ggplot() +
  
  annotate(
    "label",
    x = 1,
    y = 1,
    label = "59 genes\nSpatial stage-associated genes",
    size = 5,
    fontface = "bold",
    fill = "#D9EAF7",
    label.size = 0.6
  ) +
  
  annotate(
    "label",
    x = 3,
    y = 1,
    label = "18 genes\nPPI-derived hub genes",
    size = 5,
    fontface = "bold",
    fill = "#FCE4D6",
    label.size = 0.6
  ) +
  
  annotate(
    "text",
    x = 2,
    y = 1,
    label = "shared\nintersection",
    size = 4.5,
    fontface = "bold"
  ) +
  
  annotate(
    "label",
    x = 2,
    y = -0.25,
    label = paste0(
      "12 shared stage-associated hub genes\n\n",
      shared_gene_text
    ),
    size = 4,
    fontface = "bold",
    fill = "#E2F0D9",
    label.size = 0.7
  ) +
  
  coord_cartesian(
    xlim = c(0.2, 3.8),
    ylim = c(-1.1, 1.55),
    clip = "off"
  ) +
  
  labs(
    title = "Identification of shared stage-associated PPI hub genes"
  ) +
  
  theme_void(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

print(p_overlap)

ggsave(
  file.path(
    results_dir,
    "Figure_Overlap_59StageAssociated_18PPIHubs_12Shared.png"
  ),
  p_overlap,
  width = 10,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    results_dir,
    "Figure_Overlap_59StageAssociated_18PPIHubs_12Shared.pdf"
  ),
  p_overlap,
  width = 10,
  height = 7
)

# --------------------------------------------------
# 3. Create heatmap of the 12 shared genes
# Expression is row-scaled to Z-scores
# --------------------------------------------------

shared12_genes <- shared12_table$Gene

# Order samples: Normal -> ADM -> PDAC
sample_order <- order(
  factor(
    mapping2$Condition,
    levels = c("Normal", "ADM", "PDAC")
  )
)

sample_info <- mapping2[
  sample_order,
  c("Expression_Column", "Condition")
]

expr12 <- spatial_expr[
  shared12_genes,
  sample_info$Expression_Column,
  drop = FALSE
]

# Row-wise Z-score scaling
expr12_z <- t(scale(t(expr12)))

# Convert matrix to long format for ggplot
heat_df <- data.frame(
  Gene = rep(rownames(expr12_z), times = ncol(expr12_z)),
  Sample = rep(colnames(expr12_z), each = nrow(expr12_z)),
  Z_score = as.vector(expr12_z),
  stringsAsFactors = FALSE
)

sample_condition <- setNames(
  as.character(sample_info$Condition),
  sample_info$Expression_Column
)

heat_df$Condition <- factor(
  sample_condition[heat_df$Sample],
  levels = c("Normal", "ADM", "PDAC")
)

heat_df$Sample <- factor(
  heat_df$Sample,
  levels = colnames(expr12_z)
)

heat_df$Gene <- factor(
  heat_df$Gene,
  levels = rev(shared12_genes)
)

p_heatmap <- ggplot(
  heat_df,
  aes(
    x = Sample,
    y = Gene,
    fill = Z_score
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.15
  ) +
  
  facet_grid(
    . ~ Condition,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Row Z-score"
  ) +
  
  labs(
    title = "Spatial expression heatmap of 12 shared stage-associated PPI hub genes",
    subtitle = "Expression is standardized within each gene across Normal, ADM, and PDAC regions",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    panel.spacing.x = grid::unit(0.4, "lines")
  )

print(p_heatmap)

ggsave(
  file.path(
    results_dir,
    "Figure_Heatmap_12_Shared_StageAssociated_PPIHub_Genes.png"
  ),
  p_heatmap,
  width = 13,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    results_dir,
    "Figure_Heatmap_12_Shared_StageAssociated_PPIHub_Genes.pdf"
  ),
  p_heatmap,
  width = 13,
  height = 7
)

# ==================================================
# Clean Venn diagram: no cut labels
# ==================================================

library(VennDiagram)
library(grid)

results_dir <- "results"

spatial_stage_genes <- unique(stage_genes_corrected)
ppi_hub_genes <- unique(hub_genes)
shared_genes <- sort(intersect(spatial_stage_genes, ppi_hub_genes))

# Create diagram without automatic category labels
venn_plot <- draw.pairwise.venn(
  area1 = length(spatial_stage_genes),   # 59
  area2 = length(ppi_hub_genes),         # 18
  cross.area = length(shared_genes),     # 12
  
  category = c("", ""),
  
  euler.d = TRUE,
  scaled = FALSE,
  
  fill = c("#A6CEE3", "#FDBF6F"),
  alpha = c(0.60, 0.60),
  
  lwd = 2,
  lty = "solid",
  
  cex = 1.5,
  fontface = "bold",
  
  ind = FALSE
)

draw_clean_venn <- function() {
  
  grid.newpage()
  grid.draw(venn_plot)
  
  # Category labels: move upward
  grid.text(
    "Spatial stage-associated genes\n(n = 59)",
    x = 0.30,
    y = 0.85,
    gp = gpar(fontsize = 13, fontface = "bold", col = "#1F4E79")
  )
  
  grid.text(
    "PPI-derived hub genes\n(n = 18)",
    x = 0.70,
    y = 0.85,
    gp = gpar(fontsize = 13, fontface = "bold", col = "#8A4B00")
  )
  
  
  grid.text(
    paste0(
      "Shared genes (n = 12):\n",
      paste(shared_genes, collapse = ", ")
    ),
    x = 0.5,
    y = 0.08,
    gp = gpar(fontsize = 11)
  )
}

# Preview in RStudio
draw_clean_venn()

# Save PNG
png(
  file.path(
    results_dir,
    "Figure_Venn_59Spatial_18Hubs_12Shared_Clean.png"
  ),
  width = 3000,
  height = 2300,
  res = 300
)

draw_clean_venn()
dev.off()

# Save PDF
pdf(
  file.path(
    results_dir,
    "Figure_Venn_59Spatial_18Hubs_12Shared_Clean.pdf"
  ),
  width = 10,
  height = 7.5
)

draw_clean_venn()
dev.off()


# ==================================================
# Supplementary table:
# 12 shared stage-associated PPI hub genes
# ==================================================

results_dir <- "results"

# --------------------------------------------------
# 1. Start with the corrected 12-gene summary table
# --------------------------------------------------

supp_table_12 <- core_hub_summary[, c(
  "Gene",
  "Mean_Normal",
  "Mean_ADM",
  "Mean_PDAC",
  "ANOVA_p",
  "ANOVA_BH_FDR_80",
  "Spearman_rho",
  "Spearman_p",
  "Spearman_BH_FDR_80",
  "Stage_Associated",
  "Monotonic_Trend"
)]

# Better manuscript-friendly classification
supp_table_12$Classification <- ifelse(
  supp_table_12$Monotonic_Trend,
  "Core progression-associated hub gene",
  "Stage-associated hub gene"
)

# Identify expression pattern from group means
supp_table_12$Peak_Expression_Group <- c(
  "Normal", "ADM", "PDAC"
)[
  max.col(
    as.matrix(
      supp_table_12[, c("Mean_Normal", "Mean_ADM", "Mean_PDAC")]
    ),
    ties.method = "first"
  )
]

supp_table_12$Expression_Pattern <- ifelse(
  supp_table_12$Peak_Expression_Group == "PDAC" &
    supp_table_12$Monotonic_Trend,
  "PDAC-enriched with significant positive ordinal association",
  
  ifelse(
    supp_table_12$Peak_Expression_Group == "ADM" &
      supp_table_12$Monotonic_Trend,
    "ADM-enriched with significant positive ordinal association",
    
    "Stage-associated without significant ordinal association"
  )
)

# --------------------------------------------------
# 2. Function to load PPI rank files safely
# --------------------------------------------------

read_ppi_rank <- function(path, rank_label) {
  
  if (!file.exists(path)) {
    message("File not found: ", path)
    return(NULL)
  }
  
  df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names_lower <- tolower(names(df))
  
  # Find a likely gene/name column
  gene_col <- grep(
    "gene|symbol|node|name",
    names_lower,
    value = TRUE
  )
  
  if (length(gene_col) == 0) {
    gene_col <- names(df)[1]
  } else {
    gene_col <- names(df)[match(gene_col[1], names_lower)]
  }
  
  # Find rank column
  rank_col <- grep(
    "^rank$|rank",
    names_lower,
    value = TRUE
  )
  
  if (length(rank_col) > 0) {
    rank_col <- names(df)[match(rank_col[1], names_lower)]
    
    output <- data.frame(
      Gene = toupper(trimws(as.character(df[[gene_col]]))),
      Rank = as.numeric(df[[rank_col]]),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    # If rank is not present, search for score/degree/MCC column
    score_col <- grep(
      "score|mcc|degree",
      names_lower,
      value = TRUE
    )
    
    if (length(score_col) == 0) {
      message("No rank or score column found in: ", path)
      return(NULL)
    }
    
    score_col <- names(df)[match(score_col[1], names_lower)]
    
    output <- data.frame(
      Gene = toupper(trimws(as.character(df[[gene_col]]))),
      Score = as.numeric(df[[score_col]]),
      stringsAsFactors = FALSE
    )
    
    output <- output[order(-output$Score), ]
    output$Rank <- seq_len(nrow(output))
    output$Score <- NULL
  }
  
  names(output)[names(output) == "Rank"] <- rank_label
  
  output <- output[!duplicated(output$Gene), ]
  
  return(output)
}

# --------------------------------------------------
# 3. Load MCC and degree rankings, if available
# --------------------------------------------------

mcc_rank <- read_ppi_rank(
  "results/MCC_top20_80genes_rank.csv",
  "MCC_PPI_Rank"
)

degree_rank <- read_ppi_rank(
  "results/Degree_top20_80genes.csv",
  "Degree_PPI_Rank"
)

# Add MCC rank
if (!is.null(mcc_rank)) {
  supp_table_12 <- merge(
    supp_table_12,
    mcc_rank,
    by = "Gene",
    all.x = TRUE,
    sort = FALSE
  )
}

# Add degree rank
if (!is.null(degree_rank)) {
  supp_table_12 <- merge(
    supp_table_12,
    degree_rank,
    by = "Gene",
    all.x = TRUE,
    sort = FALSE
  )
}

# --------------------------------------------------
# 4. Keep original 12-gene order and round values
# --------------------------------------------------

supp_table_12 <- supp_table_12[
  match(core_hub_summary$Gene, supp_table_12$Gene),
]

rownames(supp_table_12) <- NULL

supp_table_12$Mean_Normal <- round(supp_table_12$Mean_Normal, 3)
supp_table_12$Mean_ADM <- round(supp_table_12$Mean_ADM, 3)
supp_table_12$Mean_PDAC <- round(supp_table_12$Mean_PDAC, 3)

supp_table_12$ANOVA_p <- signif(supp_table_12$ANOVA_p, 4)
supp_table_12$ANOVA_BH_FDR_80 <- signif(
  supp_table_12$ANOVA_BH_FDR_80,
  4
)

supp_table_12$Spearman_rho <- round(
  supp_table_12$Spearman_rho,
  3
)

supp_table_12$Spearman_p <- signif(
  supp_table_12$Spearman_p,
  4
)

supp_table_12$Spearman_BH_FDR_80 <- signif(
  supp_table_12$Spearman_BH_FDR_80,
  4
)

# Rename columns for paper table
names(supp_table_12)[names(supp_table_12) == "Mean_Normal"] <- 
  "Mean_Expression_Normal"

names(supp_table_12)[names(supp_table_12) == "Mean_ADM"] <- 
  "Mean_Expression_ADM"

names(supp_table_12)[names(supp_table_12) == "Mean_PDAC"] <- 
  "Mean_Expression_PDAC"

names(supp_table_12)[names(supp_table_12) == "ANOVA_BH_FDR_80"] <- 
  "ANOVA_BH_FDR"

names(supp_table_12)[names(supp_table_12) == "Spearman_BH_FDR_80"] <- 
  "Spearman_BH_FDR"

# Show final table
print(supp_table_12)

# Save final supplementary table
write.csv(
  supp_table_12,
  file.path(
    results_dir,
    "Supplementary_Table_12_Shared_StageAssociated_PPI_Hub_Genes.csv"
  ),
  row.names = FALSE
)

# Optional Excel-friendly tab-separated version
write.table(
  supp_table_12,
  file.path(
    results_dir,
    "Supplementary_Table_12_Shared_StageAssociated_PPI_Hub_Genes.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ==================================================
# Robust reader for MCC / Degree ranking files
# ==================================================

read_ppi_rank_safe <- function(path, rank_label) {
  
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  
  if (length(lines) == 0) {
    stop("The file is empty: ", path)
  }
  
  # Detect the separator
  test_lines <- lines[seq_len(min(10, length(lines)))]
  
  comma_n <- sum(nchar(gsub("[^,]", "", test_lines)))
  tab_n   <- sum(nchar(gsub("[^\t]", "", test_lines)))
  semi_n  <- sum(nchar(gsub("[^;]", "", test_lines)))
  
  sep <- c("," = ",", "\t" = "\t", ";" = ";")[
    which.max(c(comma_n, tab_n, semi_n))
  ]
  
  # Read safely without assuming correct header formatting
  x <- read.table(
    text = lines,
    header = FALSE,
    sep = sep,
    fill = TRUE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Remove fully empty columns
  keep_cols <- sapply(x, function(z) {
    any(!is.na(z) & trimws(as.character(z)) != "")
  })
  
  x <- x[, keep_cols, drop = FALSE]
  
  # Check whether the first row is a header
  first_row <- tolower(trimws(as.character(unlist(x[1, ]))))
  
  header_words <- c(
    "rank", "gene", "symbol", "node",
    "name", "mcc", "degree", "score"
  )
  
  first_row_is_header <- any(
    vapply(
      header_words,
      function(word) any(grepl(word, first_row)),
      logical(1)
    )
  )
  
  if (first_row_is_header) {
    
    colnames(x) <- make.unique(
      trimws(as.character(unlist(x[1, ])))
    )
    
    x <- x[-1, , drop = FALSE]
    
  } else {
    
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  
  # Convert blank strings into NA
  x[] <- lapply(x, function(z) {
    z <- trimws(as.character(z))
    z[z == ""] <- NA
    z
  })
  
  name_low <- tolower(names(x))
  
  # Find gene column
  gene_col <- grep(
    "gene|symbol|node|name|target",
    name_low,
    value = TRUE
  )
  
  # If no obvious gene column exists, use first non-numeric column
  if (length(gene_col) == 0) {
    
    numeric_fraction <- sapply(x, function(z) {
      mean(!is.na(suppressWarnings(as.numeric(z))))
    })
    
    gene_col <- names(x)[which.min(numeric_fraction)]
    
  } else {
    gene_col <- gene_col[1]
  }
  
  # Find rank column
  rank_col <- grep(
    "^rank$|rank",
    name_low,
    value = TRUE
  )
  
  # Find score column
  score_col <- grep(
    "mcc|degree|score",
    name_low,
    value = TRUE
  )
  
  output <- data.frame(
    Gene = toupper(trimws(x[[gene_col]])),
    stringsAsFactors = FALSE
  )
  
  output <- output[
    !is.na(output$Gene) &
      output$Gene != "" &
      output$Gene != "GENE",
    ,
    drop = FALSE
  ]
  
  if (length(rank_col) > 0) {
    
    output$Rank <- suppressWarnings(
      as.numeric(x[[rank_col[1]]])
    )
    
  } else if (length(score_col) > 0) {
    
    output$Score <- suppressWarnings(
      as.numeric(x[[score_col[1]]])
    )
    
    output <- output[order(-output$Score), , drop = FALSE]
    output$Rank <- seq_len(nrow(output))
    output$Score <- NULL
    
  } else {
    
    output$Rank <- seq_len(nrow(output))
  }
  
  output <- output[!duplicated(output$Gene), , drop = FALSE]
  
  names(output)[names(output) == "Rank"] <- rank_label
  
  return(output)
}

# Read MCC ranking file
mcc_rank <- read_ppi_rank_safe(
  "results/MCC_top20_80genes_rank.csv",
  "MCC_PPI_Rank"
)

print(mcc_rank)

# Read Degree ranking file
degree_rank <- read_ppi_rank_safe(
  "results/Degree_top20_80genes.csv",
  "Degree_PPI_Rank"
)

print(degree_rank)


# ==================================================
# Read Cytoscape Degree and MCC ranking files correctly
# ==================================================

read_cytoscape_rank <- function(path, rank_name, score_name) {
  
  x <- read.csv(
    path,
    skip = 1,              # Skip descriptive first line
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Expected columns: Rank, Name, Score
  out <- data.frame(
    Gene = toupper(trimws(x$Name)),
    Rank = as.numeric(x$Rank),
    Score = as.numeric(x$Score),
    stringsAsFactors = FALSE
  )
  
  names(out)[names(out) == "Rank"] <- rank_name
  names(out)[names(out) == "Score"] <- score_name
  
  out
}

# Read both ranking files
mcc_rank <- read_cytoscape_rank(
  "results/MCC_top20_80genes_rank.csv",
  "MCC_PPI_Rank",
  "MCC_Score"
)

degree_rank <- read_cytoscape_rank(
  "results/Degree_top20_80genes.csv",
  "Degree_PPI_Rank",
  "Degree_Score"
)

print(mcc_rank)
print(degree_rank)

# ==================================================
# Final supplementary table for 12 shared genes
# ==================================================

supp_table_12 <- core_hub_summary[, c(
  "Gene",
  "Mean_Normal",
  "Mean_ADM",
  "Mean_PDAC",
  "ANOVA_p",
  "ANOVA_BH_FDR_80",
  "Spearman_rho",
  "Spearman_p",
  "Spearman_BH_FDR_80",
  "Stage_Associated",
  "Monotonic_Trend"
)]

# Classification
supp_table_12$Classification <- ifelse(
  supp_table_12$Monotonic_Trend,
  "Core progression-associated hub gene",
  "Stage-associated hub gene"
)

# Add MCC values
mcc_match <- match(supp_table_12$Gene, mcc_rank$Gene)

supp_table_12$MCC_PPI_Rank <- mcc_rank$MCC_PPI_Rank[mcc_match]
supp_table_12$MCC_Score <- mcc_rank$MCC_Score[mcc_match]

# Add Degree values
degree_match <- match(supp_table_12$Gene, degree_rank$Gene)

supp_table_12$Degree_PPI_Rank <- degree_rank$Degree_PPI_Rank[degree_match]
supp_table_12$Degree_Score <- degree_rank$Degree_Score[degree_match]

# Round values for paper readability
supp_table_12$Mean_Normal <- round(supp_table_12$Mean_Normal, 2)
supp_table_12$Mean_ADM <- round(supp_table_12$Mean_ADM, 2)
supp_table_12$Mean_PDAC <- round(supp_table_12$Mean_PDAC, 2)

supp_table_12$ANOVA_p <- signif(supp_table_12$ANOVA_p, 4)
supp_table_12$ANOVA_BH_FDR_80 <- signif(
  supp_table_12$ANOVA_BH_FDR_80, 4
)

supp_table_12$Spearman_rho <- round(
  supp_table_12$Spearman_rho, 3
)

supp_table_12$Spearman_p <- signif(
  supp_table_12$Spearman_p, 4
)

supp_table_12$Spearman_BH_FDR_80 <- signif(
  supp_table_12$Spearman_BH_FDR_80, 4
)

# Rename for final paper table
names(supp_table_12)[names(supp_table_12) == "Mean_Normal"] <- "Mean_Normal"
names(supp_table_12)[names(supp_table_12) == "Mean_ADM"] <- "Mean_ADM"
names(supp_table_12)[names(supp_table_12) == "Mean_PDAC"] <- "Mean_PDAC"
names(supp_table_12)[names(supp_table_12) == "ANOVA_BH_FDR_80"] <- "ANOVA_BH_FDR"
names(supp_table_12)[names(supp_table_12) == "Spearman_BH_FDR_80"] <- "Spearman_BH_FDR"

print(supp_table_12)

write.csv(
  supp_table_12,
  "results/Supplementary_Table_12_Shared_StageAssociated_PPI_Hub_Genes.csv",
  row.names = FALSE
)


# ==================================================
# Analyze pancreas-specific miRNA network
# for all 12 shared stage-associated PPI hub genes
# ==================================================

library(igraph)

graph_path <- "results/mirnet12genes.graphml"

g12 <- read_graph(
  graph_path,
  format = "graphml"
)

vertex12 <- as_data_frame(g12, what = "vertices")
edge12 <- as_data_frame(g12, what = "edges")

names(vertex12)
names(edge12)

head(vertex12)
head(edge12)


# ==================================================
# Filter 12-gene pancreas-specific miRNet network
# by 17 EV candidates and final 6-miRNA panel
# ==================================================

library(igraph)

# The 12 shared stage-associated PPI hub genes
shared12_genes <- unique(toupper(core_hub_summary$Gene))

# Your 17 EV-associated miRNA candidates
ev17_list <- c(
  "miR-15a-5p",
  "miR-20a-5p",
  "miR-21-5p",
  "miR-98-5p",
  "miR-107",
  "miR-34a-5p",
  "miR-182-5p",
  "miR-424-5p",
  "miR-184",
  "miR-630",
  "miR-216b-5p",
  "miR-509-3-5p",
  "miR-1908-5p",
  "miR-4756-3p",
  "miR-1306-5p",
  "miR-215-3p",
  "miR-134-3p"
)

# Final diagnostic panel
panel6_list <- c(
  "miR-184",
  "miR-107",
  "miR-216b-5p",
  "miR-20a-5p",
  "miR-134-3p",
  "miR-215-3p"
)

# Standardize miRNA names for matching
clean_mir <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^hsa-", "", x, ignore.case = TRUE)
  tolower(x)
}

# Identify miRNA nodes
is_mirna <- function(x) {
  grepl("^(hsa-)?(miR|let)-", x, ignore.case = TRUE)
}

# Keep only gene-miRNA associations involving the 12 shared genes
gene_mir_edges12 <- edge12[
  (
    toupper(edge12$from) %in% shared12_genes &
      is_mirna(edge12$to)
  ) |
    (
      toupper(edge12$to) %in% shared12_genes &
        is_mirna(edge12$from)
    ),
  ,
  drop = FALSE
]

# Create a clean Gene-miRNA table
gene_mir_edges12$Gene <- ifelse(
  toupper(gene_mir_edges12$from) %in% shared12_genes,
  toupper(gene_mir_edges12$from),
  toupper(gene_mir_edges12$to)
)

gene_mir_edges12$miRNA <- ifelse(
  is_mirna(gene_mir_edges12$from),
  gene_mir_edges12$from,
  gene_mir_edges12$to
)

gene_mir_edges12 <- unique(
  gene_mir_edges12[, c("Gene", "miRNA")]
)

# Filter using the 17 EV-miRNA candidates
shared12_ev_edges <- gene_mir_edges12[
  clean_mir(gene_mir_edges12$miRNA) %in% clean_mir(ev17_list),
  ,
  drop = FALSE
]

# Mark whether each EV-miRNA belongs to the final six-miRNA panel
shared12_ev_edges$Final_6_Panel <- clean_mir(
  shared12_ev_edges$miRNA
) %in% clean_mir(panel6_list)

shared12_ev_edges$Category <- ifelse(
  shared12_ev_edges$Final_6_Panel,
  "Final six-miRNA panel",
  "EV candidate only"
)

# Keep only final-six interactions
shared12_final6_edges <- subset(
  shared12_ev_edges,
  Final_6_Panel == TRUE
)

# Order tables
shared12_ev_edges <- shared12_ev_edges[
  order(
    shared12_ev_edges$Gene,
    shared12_ev_edges$Category,
    shared12_ev_edges$miRNA
  ),
]

shared12_final6_edges <- shared12_final6_edges[
  order(
    shared12_final6_edges$Gene,
    shared12_final6_edges$miRNA
  ),
]

# --------------------------------------------------
# Gene-wise summary: objective ranking
# --------------------------------------------------

gene_mir_summary <- data.frame(
  Gene = shared12_genes,
  EV_Candidate_miRNA_Count = integer(length(shared12_genes)),
  Final6_miRNA_Count = integer(length(shared12_genes)),
  EV_Candidate_miRNAs = character(length(shared12_genes)),
  Final6_miRNAs = character(length(shared12_genes)),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(gene_mir_summary))) {
  
  this_gene <- gene_mir_summary$Gene[i]
  
  ev_mirs <- sort(unique(
    shared12_ev_edges$miRNA[
      shared12_ev_edges$Gene == this_gene
    ]
  ))
  
  final_mirs <- sort(unique(
    shared12_final6_edges$miRNA[
      shared12_final6_edges$Gene == this_gene
    ]
  ))
  
  gene_mir_summary$EV_Candidate_miRNA_Count[i] <- length(ev_mirs)
  gene_mir_summary$Final6_miRNA_Count[i] <- length(final_mirs)
  
  gene_mir_summary$EV_Candidate_miRNAs[i] <- ifelse(
    length(ev_mirs) == 0,
    "None",
    paste(ev_mirs, collapse = ", ")
  )
  
  gene_mir_summary$Final6_miRNAs[i] <- ifelse(
    length(final_mirs) == 0,
    "None",
    paste(final_mirs, collapse = ", ")
  )
}

gene_mir_summary <- gene_mir_summary[
  order(
    -gene_mir_summary$Final6_miRNA_Count,
    -gene_mir_summary$EV_Candidate_miRNA_Count,
    gene_mir_summary$Gene
  ),
]

# --------------------------------------------------
# miRNA-wise summary
# --------------------------------------------------

mirna_summary <- aggregate(
  Gene ~ miRNA + Category,
  data = shared12_ev_edges,
  FUN = function(x) paste(sort(unique(x)), collapse = ", ")
)

mirna_count <- aggregate(
  Gene ~ miRNA + Category,
  data = shared12_ev_edges,
  FUN = function(x) length(unique(x))
)

names(mirna_count)[names(mirna_count) == "Gene"] <- "Shared_Genes_Connected"

mirna_summary <- merge(
  mirna_summary,
  mirna_count,
  by = c("miRNA", "Category"),
  all.x = TRUE
)

names(mirna_summary)[names(mirna_summary) == "Gene"] <- "Connected_Genes"

mirna_summary <- mirna_summary[
  order(
    -mirna_summary$Shared_Genes_Connected,
    mirna_summary$miRNA
  ),
]

# --------------------------------------------------
# Save results
# --------------------------------------------------

results_dir <- "results"

write.csv(
  shared12_ev_edges,
  file.path(
    results_dir,
    "Shared12_Genes_17EVmiRNA_Associations.csv"
  ),
  row.names = FALSE
)

write.csv(
  shared12_final6_edges,
  file.path(
    results_dir,
    "Shared12_Genes_Final6miRNA_Associations.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_mir_summary,
  file.path(
    results_dir,
    "Shared12_Gene_miRNA_Association_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  mirna_summary,
  file.path(
    results_dir,
    "Shared12_miRNA_Association_Summary.csv"
  ),
  row.names = FALSE
)

# --------------------------------------------------
# Create Cytoscape networks
# No arrows: these are miRNet-derived associations
# --------------------------------------------------

make_network <- function(edge_table, output_file) {
  
  if (nrow(edge_table) == 0) {
    warning("No edges found for: ", output_file)
    return(NULL)
  }
  
  network_edges <- data.frame(
    from = edge_table$miRNA,
    to = edge_table$Gene,
    stringsAsFactors = FALSE
  )
  
  node_info <- data.frame(
    name = unique(c(network_edges$from, network_edges$to)),
    stringsAsFactors = FALSE
  )
  
  node_info$NodeType <- ifelse(
    node_info$name %in% shared12_genes,
    "Shared stage-associated PPI hub gene",
    "miRNA"
  )
  
  g_network <- graph_from_data_frame(
    network_edges,
    directed = FALSE,
    vertices = node_info
  )
  
  write_graph(
    g_network,
    file.path(results_dir, output_file),
    format = "graphml"
  )
  
  return(g_network)
}

g_shared12_ev <- make_network(
  shared12_ev_edges,
  "Shared12_Genes_17EVmiRNA_Network.graphml"
)

g_shared12_final6 <- make_network(
  shared12_final6_edges,
  "Shared12_Genes_Final6miRNA_Network.graphml"
)

# Main outputs to inspect
cat("\nAll 12-gene / EV-miRNA associations:\n")
print(shared12_ev_edges)

cat("\nFinal six-miRNA associations:\n")
print(shared12_final6_edges)

cat("\nGene-level ranking:\n")
print(gene_mir_summary)

# ==========================================================
# Figure 3: Spatial validation and stage-association analysis
# ==========================================================

# Install packages once, only if needed
required_packages <- c("ggplot2", "patchwork")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(patchwork)
library(grid)

# ----------------------------------------------------------
# 1. File path
# ----------------------------------------------------------

results_dir <- "results"

stats_file <- file.path(
  results_dir,
  "Corrected_80_Gene_Spatial_Progression_Statistics.csv"
)

if (!file.exists(stats_file)) {
  stop(
    "File not found: ", stats_file,
    "\nCheck your Results folder path."
  )
}

stats_df <- read.csv(
  stats_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ----------------------------------------------------------
# 2. Check required columns
# ----------------------------------------------------------

required_columns <- c(
  "Gene",
  "ANOVA_BH_FDR_80",
  "Spearman_BH_FDR_80",
  "Stage_Associated"
)

missing_columns <- setdiff(required_columns, colnames(stats_df))

if (length(missing_columns) > 0) {
  stop(
    "These columns are missing from the CSV file: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Convert logical columns safely
stats_df$Stage_Associated <- as.logical(stats_df$Stage_Associated)

# Core progression-associated genes:
# significant in both ANOVA and ordinal Spearman analyses
stats_df$Core_Progression <- (
  stats_df$Stage_Associated == TRUE &
    stats_df$Spearman_BH_FDR_80 < 0.05
)

# ----------------------------------------------------------
# 3. Calculate counts directly from your corrected results
# ----------------------------------------------------------

total_genes <- nrow(stats_df)

stage_associated_n <- sum(
  stats_df$Stage_Associated == TRUE,
  na.rm = TRUE
)

not_stage_associated_n <- sum(
  stats_df$Stage_Associated == FALSE,
  na.rm = TRUE
)

core_genes_n <- sum(
  stats_df$Core_Progression == TRUE,
  na.rm = TRUE
)

stage_only_n <- stage_associated_n - core_genes_n

cat("\nCorrected spatial-analysis counts:\n")
cat("Total spatially validated genes:", total_genes, "\n")
cat("Stage-associated genes:", stage_associated_n, "\n")
cat("Core progression-associated genes:", core_genes_n, "\n")
cat("Stage-associated only:", stage_only_n, "\n")
cat("Not stage-associated:", not_stage_associated_n, "\n")

# Expected values:
# Total = 80
# Stage-associated = 59
# Core = 9
# Stage-associated only = 50
# Not stage-associated = 21

# ----------------------------------------------------------
# 4. Shared journal-style theme
# ----------------------------------------------------------

figure_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    plot.margin = margin(
      t = 16,
      r = 12,
      b = 16,
      l = 12
    )
  )

# ----------------------------------------------------------
# Panel A: Spatial cohort composition
# ----------------------------------------------------------

cohort_df <- data.frame(
  Group = c(
    "Adjacent normal\n(n = 16)",
    "ADM\n(n = 16)",
    "PDAC\n(n = 16)"
  ),
  y_position = c(3, 2, 1),
  Fill = c("#5B9BD5", "#F4B183", "#C55A5A"),
  TextColor = c("white", "black", "white")
)

panel_A <- ggplot(cohort_df, aes(x = 1, y = y_position)) +
  geom_label(
    aes(
      label = Group,
      fill = Fill,
      color = TextColor
    ),
    label.size = 0,
    size = 5.2,
    fontface = "bold",
    label.r = unit(0.18, "lines"),
    label.padding = unit(0.65, "lines"),
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 1,
    y = 4.05,
    label = "GSE208536 GeoMx cohort",
    fontface = "bold",
    size = 5.2
  ) +
  annotate(
    "text",
    x = 1,
    y = 0.15,
    label = paste0(
      "48 spatial regions\n",
      "verified by matching GeoMx segment identifiers\n",
      "with GEO sample metadata"
    ),
    size = 3.6,
    color = "grey30"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(
    xlim = c(0.35, 1.65),
    ylim = c(0, 4.35),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(16, 8, 14, 8)
  )

# ----------------------------------------------------------
# Panel B: Clean spatial-validation workflow
# ----------------------------------------------------------
# ----------------------------------------------------------
# Panel B: Clean spatial-validation workflow
# ----------------------------------------------------------

panel_B <- ggplot() +
  
  # Box 1
  annotate(
    "label",
    x = 1.25,
    y = 2.85,
    label = "628\nmalignant-switch genes\nfrom GSE143754",
    fill = "#FFF2CC",
    color = "#C98500",
    size = 4.0,
    fontface = "bold",
    label.r = unit(0.15, "lines"),
    label.padding = unit(0.55, "lines")
  ) +
  
  # Box 2
  annotate(
    "label",
    x = 4.25,
    y = 2.85,
    label = "GeoMx expression\nplatform coverage",
    fill = "#EAF2FB",
    color = "#2F75C1",
    size = 4.0,
    fontface = "bold",
    label.r = unit(0.15, "lines"),
    label.padding = unit(0.55, "lines")
  ) +
  
  # Box 3
  annotate(
    "label",
    x = 2.75,
    y = 0.95,
    label = "80 spatially validated\ngenes",
    fill = "#EAF2FB",
    color = "#2F75C1",
    size = 4.8,
    fontface = "bold",
    label.r = unit(0.15, "lines"),
    label.padding = unit(0.70, "lines")
  ) +
  
  # Small straight arrow: box 1 -> box 2
  annotate(
    "segment",
    x = 2.10,
    xend = 3.38,
    y = 2.85,
    yend = 2.85,
    linewidth = 0.7,
    arrow = arrow(
      length = unit(0.14, "cm"),
      type = "closed"
    )
  ) +
  
  # Smaller curved arrow: box 2 -> box 3
  annotate(
    "curve",
    x = 4.05,
    y = 2.20,
    xend = 3.15,
    yend = 1.45,
    curvature = 0.20,
    linewidth = 0.7,
    arrow = arrow(
      length = unit(0.14, "cm"),
      type = "closed"
    )
  ) +
  
  # Title
  annotate(
    "text",
    x = 2.75,
    y = 4.05,
    label = "Spatial validation",
    fontface = "bold",
    size = 5.2
  ) +
  
  # Bottom text
  annotate(
    "text",
    x = 2.75,
    y = 0.10,
    label = "Genes were retained only when represented\non the GeoMx expression platform.",
    size = 3.6,
    color = "grey30"
  ) +
  
  coord_cartesian(
    xlim = c(0.1, 5.4),
    ylim = c(0, 4.35),
    clip = "off"
  ) +
  
  theme_void() +
  theme(
    plot.margin = margin(16, 12, 14, 12)
  )
# ----------------------------------------------------------
# Panel C: ANOVA stage-association result
# ----------------------------------------------------------

anova_summary <- data.frame(
  Status = factor(
    c(
      "Stage-associated\n(ANOVA FDR < 0.05)",
      "Not stage-associated"
    ),
    levels = c(
      "Stage-associated\n(ANOVA FDR < 0.05)",
      "Not stage-associated"
    )
  ),
  Count = c(
    stage_associated_n,
    not_stage_associated_n
  )
)

panel_C <- ggplot(
  anova_summary,
  aes(x = Status, y = Count, fill = Status)
) +
  geom_col(
    width = 0.8,
    color = "grey25",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = Count),
    vjust = -0.45,
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Stage-associated\n(ANOVA FDR < 0.05)" = "#70AD47",
      "Not stage-associated" = "#BFBFBF"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 80, 10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "ANOVA stage-association results",
    x = NULL,
    y = "Number of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 10.5,
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 12
    ),
    plot.margin = margin(16, 8, 40, 8)
  ) +
  annotate(
    "text",
    x = 1.5,
    y = -13,
    label = "Benjaminiâ€“Hochberg correction across\n80 spatially validated genes.",
    size = 3.5,
    color = "grey30"
  ) +
  coord_cartesian(clip = "off")
# ----------------------------------------------------------
# Panel D: Ordinal association classification
# ----------------------------------------------------------

ordinal_summary <- data.frame(
  Status = factor(
    c(
      "Core progression-associated\n(ANOVA + Spearman FDR)",
      "Stage-associated only\n(ANOVA FDR only)",
      "Not stage-associated"
    ),
    levels = c(
      "Core progression-associated\n(ANOVA + Spearman FDR)",
      "Stage-associated only\n(ANOVA FDR only)",
      "Not stage-associated"
    )
  ),
  Count = c(
    9,
    50,
    21
  )
)

panel_D <- ggplot(
  ordinal_summary,
  aes(x = Status, y = Count, fill = Status)
) +
  
  geom_col(
    width = 0.78,
    color = "grey25",
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(label = Count),
    vjust = -0.45,
    size = 6,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Core progression-associated\n(ANOVA + Spearman FDR)" = "#1F8E9E",
      "Stage-associated only\n(ANOVA FDR only)" = "#9DC3E6",
      "Not stage-associated" = "#D0D0D0"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 64),
    breaks = seq(0, 60, 10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    title = "Ordinal association classification",
    x = NULL,
    y = "Number of genes"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.x = element_text(
      size = 9.2,
      color = "black",
      lineheight = 0.95
    ),
    
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    plot.margin = margin(
      t = 16,
      r = 10,
      b = 18,
      l = 10
    )
  ) +
  
  annotate(
    "text",
    x = 2,
    y = 61.5,
    label = "9 genes were significant in both ANOVA FDR\nand ordinal Spearman FDR.",
    size = 3.6,
    color = "grey30"
  )

# ----------------------------------------------------------
# 5. Combine Panels A-D
# ----------------------------------------------------------

figure_3 <- (
  panel_A +
    panel_B +
    panel_C +
    panel_D +
    plot_layout(widths = c(1.0, 1.15, 1.18, 1.25))
) +
  plot_annotation(
    tag_levels = "A"
  )

print(figure_3)
  

# ----------------------------------------------------------
# 6. Save final figures
# ----------------------------------------------------------

ggsave(
  filename = file.path(
    results_dir,
    "Figure_3_Spatial_Validation_and_Stage_Association.png"
  ),
  plot = figure_3,
  width = 18,
  height = 7.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    results_dir,
    "Figure_3_Spatial_Validation_and_Stage_Association.pdf"
  ),
  plot = figure_3,
  width = 18,
  height = 7.5,
  bg = "white"
)

print(figure_3)

