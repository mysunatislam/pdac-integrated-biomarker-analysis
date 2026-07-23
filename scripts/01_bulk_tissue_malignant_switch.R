# Recovered effective analysis script
# Source: R_proj\Proj_Panc\Results\phase1.R
# Paths were normalized from the original local D:/R_proj layout to project-relative paths.
# The script is preserved from the working project history and may require the listed R packages
# plus public GEO downloads to rerun fully from scratch.

rm(list = ls())
gc()
library(oligo)
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hta20transcriptcluster.db)
library(ggplot2)
library(pheatmap)

celFiles <- list.files(
  "./GSE143754_RAW",
  pattern = "CEL.gz$",
  full.names = TRUE
)

length(celFiles)

rawData <- read.celfiles(celFiles)
rawData


normData <- rma(rawData)

expr <- exprs(normData)

dim(expr)

saveRDS(
  normData,
  "./normData_HTA2.rds"
)

gse <- getGEO(
  "GSE143754",
  GSEMatrix = TRUE
)

pheno <- pData(gse[[1]])

sample_ids <- sub(
  "_.*",
  "",
  colnames(expr)
)

head(sample_ids)

pheno <- pheno[sample_ids, ]

all(sample_ids == rownames(pheno))

table(
  pheno$`disease state:ch1`
)

keep <- pheno$`disease state:ch1` %in%
  c(
    "Adjacent Normal",
    "Tumor"
  )

expr2 <- expr[, keep]

group <- factor(
  pheno$`disease state:ch1`[keep],
  levels = c(
    "Adjacent Normal",
    "Tumor"
  )
)

table(group)


dim(expr2)
length(group)

design <- model.matrix(
  ~0 + group
)

colnames(design) <- c(
  "Normal",
  "Tumor"
)

design

fit <- lmFit(
  expr2,
  design
)

cont.matrix <- makeContrasts(
  Tumor - Normal,
  levels = design
)

fit2 <- contrasts.fit(
  fit,
  cont.matrix
)

fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

head(deg)

sig <- subset(
  deg,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

nrow(sig)

table(
  ifelse(
    sig$logFC > 0,
    "Up",
    "Down"
  )
)

probe_ids <- rownames(sig)

annot <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = probe_ids,
  columns = c(
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "PROBEID"
)

sig$PROBEID <- rownames(sig)

sig_annot <- merge(
  sig,
  annot,
  by = "PROBEID",
  all.x = TRUE
)

sig_genes <- sig_annot[
  !is.na(sig_annot$SYMBOL),
]

sig_genes <- sig_genes[
  order(sig_genes$adj.P.Val),
]

sig_genes <- sig_genes[
  !duplicated(sig_genes$SYMBOL),
]

write.csv(
  sig_genes,
  "./GSE143754_DEGs.csv",
  row.names = FALSE
)

volc <- deg

volc$Status <- "NS"

volc$Status[
  volc$adj.P.Val < 0.05 &
    volc$logFC > 1
] <- "Up"

volc$Status[
  volc$adj.P.Val < 0.05 &
    volc$logFC < -1
] <- "Down"

p <- ggplot(
  volc,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Status
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      Up = "red",
      Down = "blue",
      NS = "grey"
    )
  ) +
  theme_bw() +
  labs(
    title = "GSE143754",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P"
  )

p

ggsave(
  "./Volcano_GSE143754.png",
  p,
  width = 8,
  height = 6,
  dpi = 300
)

top50 <- sig_genes$PROBEID[1:50]

mat <- expr2[
  top50,
]

mat <- t(
  scale(
    t(mat)
  )
)

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(mat)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = FALSE
)

png(
  "./Heatmap_top50.png",
  width = 1200,
  height = 1000
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = FALSE
)

dev.off()

dim(expr2)
table(group)

all(sample_ids == rownames(pheno))
save.image(
  "./Panc_Workspace.RData"
)


group_all <- pheno$`disease state:ch1`

keep <- group_all %in% c(
  "Adjacent Normal",
  "Chronic Pancreatitis"
)

table(keep)

expr_cp <- expr[, keep]

group_cp <- factor(
  group_all[keep],
  levels = c(
    "Adjacent Normal",
    "Chronic Pancreatitis"
  )
)

table(group_cp)

dim(expr_cp)
length(group_cp)

library(limma)

design_cp <- model.matrix(~0 + group_cp)

colnames(design_cp) <- c(
  "Normal",
  "CP"
)

design_cp


fit_cp <- lmFit(
  expr_cp,
  design_cp
)

cont_cp <- makeContrasts(
  CP - Normal,
  levels = design_cp
)

fit_cp2 <- contrasts.fit(
  fit_cp,
  cont_cp
)

fit_cp2 <- eBayes(
  fit_cp2
)

deg_cp <- topTable(
  fit_cp2,
  number = Inf,
  adjust.method = "BH"
)

head(deg_cp)

sig_cp <- subset(
  deg_cp,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

nrow(sig_cp)

table(
  ifelse(
    sig_cp$logFC > 0,
    "Up",
    "Down"
  )
)

probe_ids <- rownames(sig_cp)

annot_cp <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = probe_ids,
  columns = c(
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "PROBEID"
)

sig_cp$PROBEID <- rownames(sig_cp)

sig_cp_annot <- merge(
  sig_cp,
  annot_cp,
  by = "PROBEID",
  all.x = TRUE
)

sig_cp_genes <- sig_cp_annot[
  !is.na(sig_cp_annot$SYMBOL),
]

sig_cp_genes <- sig_cp_genes[
  order(sig_cp_genes$adj.P.Val),
]

sig_cp_genes <- sig_cp_genes[
  !duplicated(sig_cp_genes$SYMBOL),
]

write.csv(
  sig_cp_genes,
  "CP_vs_Normal_DEGs.csv",
  row.names = FALSE
)


library(ggplot2)

volc_cp <- deg_cp

volc_cp$Status <- "NS"

volc_cp$Status[
  volc_cp$adj.P.Val < 0.05 &
    volc_cp$logFC > 1
] <- "Up"

volc_cp$Status[
  volc_cp$adj.P.Val < 0.05 &
    volc_cp$logFC < -1
] <- "Down"

ggplot(
  volc_cp,
  aes(
    logFC,
    -log10(adj.P.Val),
    color = Status
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      Down = "blue",
      NS = "grey",
      Up = "red"
    )
  ) +
  theme_bw() +
  labs(
    title = "Chronic Pancreatitis vs Normal",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

ggsave(
  "Volcano_CP_vs_Normal.png",
  width = 8,
  height = 6,
  dpi = 300
)

library(pheatmap)

top50_cp <- sig_cp_genes$PROBEID[1:50]

mat_cp <- expr_cp[top50_cp, ]

mat_cp <- t(scale(t(mat_cp)))

annotation_col_cp <- data.frame(
  Group = group_cp
)

rownames(annotation_col_cp) <- colnames(mat_cp)

png(
  "Heatmap_CP_vs_Normal.png",
  width = 1400,
  height = 1200,
  res = 150
)

pheatmap(
  mat_cp,
  annotation_col = annotation_col_cp,
  show_rownames = FALSE,
  fontsize_col = 10,
  main = "CP vs Normal"
)

dev.off()

save(
  deg_cp,
  sig_cp,
  sig_cp_genes,
  file = "./CP_vs_Normal.RData"
)

group_all <- pheno$`disease state:ch1`

keep <- group_all %in% c(
  "Tumor",
  "Chronic Pancreatitis"
)

expr_tcp <- expr[, keep]

group_tcp <- factor(
  group_all[keep],
  levels = c(
    "Chronic Pancreatitis",
    "Tumor"
  )
)

table(group_tcp)

dim(expr_tcp)
library(limma)

design <- model.matrix(~0 + group_tcp)

colnames(design) <- make.names(levels(group_tcp))

design
fit <- lmFit(expr_tcp, design)

cont.matrix <- makeContrasts(
  Tumor - Chronic.Pancreatitis,
  levels = design
)

cont.matrix

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

deg_tcp <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

head(deg_tcp)


sig_tcp <- subset(
  deg_tcp,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

nrow(sig_tcp)

table(
  ifelse(sig_tcp$logFC > 0,
         "Up",
         "Down")
)

library(AnnotationDbi)
library(hta20transcriptcluster.db)

probe_ids <- rownames(sig_tcp)

annot <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = probe_ids,
  columns = c("SYMBOL","GENENAME"),
  keytype = "PROBEID"
)

sig_tcp$PROBEID <- rownames(sig_tcp)

sig_tcp_annot <- merge(
  sig_tcp,
  annot,
  by = "PROBEID",
  all.x = TRUE
)

sig_tcp_genes <- sig_tcp_annot[
  !is.na(sig_tcp_annot$SYMBOL),
]

sig_tcp_genes <- sig_tcp_genes[
  order(sig_tcp_genes$adj.P.Val),
]

sig_tcp_genes <- sig_tcp_genes[
  !duplicated(sig_tcp_genes$SYMBOL),
]

head(
  sig_tcp_genes[
    ,
    c("SYMBOL","logFC","adj.P.Val")
  ],
  30
)

write.csv(
  sig_tcp_genes,
  "./Tumor_vs_CP_DEGs.csv",
  row.names = FALSE
)

library(ggplot2)

volc <- deg_tcp

volc$Status <- "NS"

volc$Status[
  volc$adj.P.Val < 0.1 &
    volc$logFC > 1
] <- "Up"

volc$Status[
  volc$adj.P.Val < 0.1 &
    volc$logFC < -1
] <- "Down"

ggplot(volc,
       aes(logFC,
           -log10(adj.P.Val),
           color=Status)) +
  geom_point(alpha=0.7,size=1.2) +
  
  geom_vline(xintercept=c(-1,1),
             linetype="dashed") +
  
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed") +
  
  scale_color_manual(
    values=c(
      Down="blue",
      NS="grey80",
      Up="red"
    )
  ) +
  theme_bw(base_size=14)

ggsave(
  "./Volcano_Tumor_vs_CP.png",
  width = 8,
  height = 6,
  dpi = 300
)

library(pheatmap)

top50_tcp <- sig_tcp_genes$PROBEID[1:50]

mat_tcp <- expr_tcp[top50_tcp, ]

mat_tcp <- t(scale(t(mat_tcp)))

annotation_col <- data.frame(
  Group = group_tcp
)

rownames(annotation_col) <- colnames(mat_tcp)

png(
  "./Heatmap_Tumor_vs_CP.png",
  width = 1400,
  height = 1200,
  res = 150
)

pheatmap(
  mat_tcp,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Tumor vs Chronic Pancreatitis"
)

dev.off()

tumor_genes <- unique(sig_genes$SYMBOL)
cp_genes <- unique(sig_cp_genes$SYMBOL)

common_genes <- intersect(
  tumor_genes,
  cp_genes
)

tumor_specific <- setdiff(
  tumor_genes,
  cp_genes
)

length(common_genes)
length(tumor_specific)


write.csv(
  data.frame(SYMBOL = common_genes),
  "./Common_Tumor_CP_Genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(SYMBOL = tumor_specific),
  "./Tumor_Specific_Genes.csv",
  row.names = FALSE
)


tumor_specific_table <- sig_genes[
  sig_genes$SYMBOL %in% tumor_specific,
]

tumor_specific_table <- tumor_specific_table[
  order(tumor_specific_table$adj.P.Val),
]

head(
  tumor_specific_table[
    ,
    c("SYMBOL","logFC","adj.P.Val")
  ],
  50
)
install.packages("VennDiagram")

library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    Tumor = tumor_genes,
    Chronic_Pancreatitis = cp_genes
  ),
  filename = NULL,
  fill = c("red","blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

grid::grid.draw(venn.plot)

png(
  "./Venn_Tumor_CP.png",
  width = 1200,
  height = 1000,
  res = 150
)

grid::grid.draw(venn.plot)

dev.off()

install.packages("ggVennDiagram")
library(ggVennDiagram)

gene_tn <- unique(sig_genes$SYMBOL)
gene_tcp <- unique(sig_tcp_genes$SYMBOL)

ggVennDiagram(
  list(
    Tumor_vs_Normal = gene_tn,
    Tumor_vs_CP = gene_tcp
  )
)
ggsave(
  "./Venn_Tumor_Normal_CP.png",
  width = 7,
  height = 6,
  dpi = 300
)
install.packages("ggVennDiagram")
library(ggVennDiagram)
gene_tn  <- unique(na.omit(sig_genes$SYMBOL))
gene_cp  <- unique(na.omit(sig_cp_genes$SYMBOL))
gene_tcp <- unique(na.omit(sig_tcp_genes$SYMBOL))
length(gene_tn)
length(gene_cp)
length(gene_tcp)
library(ggVennDiagram)
library(ggVennDiagram)

p <- ggVennDiagram(
  list(
    TN = gene_tn,
    CPN = gene_cp,
    TCP = gene_tcp
  ),
  label_alpha = 0
) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(50,100,50,100)
  )

p
ggsave(
  "./Venn_3Groups.png",
  plot = p,
  width = 12,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)
install.packages("VennDiagram")
library(VennDiagram)

png(
  "./Venn_Flower.png",
  width = 2000,
  height = 1800,
  res = 300
)

draw.triple.venn(
  area1 = length(gene_tn),
  area2 = length(gene_cp),
  area3 = length(gene_tcp),
  
  n12 = length(intersect(gene_tn,gene_cp)),
  n23 = length(intersect(gene_cp,gene_tcp)),
  n13 = length(intersect(gene_tn,gene_tcp)),
  
  n123 = length(
    Reduce(
      intersect,
      list(gene_tn,gene_cp,gene_tcp)
    )
  ),
  
  category = c(
    "Tumor vs Normal",
    "CP vs Normal",
    "Tumor vs CP"
  ),
  
  fill = c(
    "red",
    "skyblue",
    "green"
  ),
  
  alpha = 0.5,
  
  cex = 2,
  cat.cex = 1.5,
  
  scaled = FALSE
)

dev.off()
ls()
gene_tn <- unique(sig_genes$SYMBOL)

gene_cp <- unique(sig_cp_genes$SYMBOL)

gene_tcp <- unique(sig_tcp_genes$SYMBOL)

common_all <- Reduce(
  intersect,
  list(
    gene_tn,
    gene_cp,
    gene_tcp
  )
)

length(common_all)

common_all

write.csv(
  data.frame(Gene = common_all),
  "./Common_3way_Genes.csv",
  row.names = FALSE
)

tn_only <- setdiff(
  gene_tn,
  union(gene_cp, gene_tcp)
)

cp_only <- setdiff(
  gene_cp,
  union(gene_tn, gene_tcp)
)

tcp_only <- setdiff(
  gene_tcp,
  union(gene_tn, gene_cp)
)

length(tn_only)
length(cp_only)
length(tcp_only)

length(gene_tn)
length(gene_cp)
length(gene_tcp)

library(pheatmap)

common_all
# 16 genes

# Get probe IDs corresponding to common genes
common_probes <- unique(
  c(
    sig_genes$PROBEID[sig_genes$SYMBOL %in% common_all],
    sig_cp_genes$PROBEID[sig_cp_genes$SYMBOL %in% common_all],
    sig_tcp_genes$PROBEID[sig_tcp_genes$SYMBOL %in% common_all]
  )
)

common_probes
length(common_probes)

# Expression matrix for all 26 samples
mat_common <- expr[common_probes, ]

# Gene symbols as row names
gene_names <- annot$SYMBOL[
  match(common_probes, annot$PROBEID)
]

rownames(mat_common) <- gene_names

# Remove duplicate symbols if present
mat_common <- mat_common[
  !duplicated(rownames(mat_common)),
]

# Z-score normalization
mat_common <- t(scale(t(mat_common)))

# Sample groups
group26 <- pheno$`disease state:ch1`

annotation_col <- data.frame(
  Group = group26
)

rownames(annotation_col) <- colnames(mat_common)

# Save heatmap
png(
  "./Heatmap_16_Common_Genes.png",
  width = 1800,
  height = 1400,
  res = 200
)

pheatmap(
  mat_common,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 12,
  fontsize_col = 8,
  main = "16 Common DEGs Across Normal, CP and Tumor"
)

dev.off()


group26 <- pheno$`disease state:ch1`

ord <- order(group26)

mat_common2 <- mat_common[, ord]

annotation_col2 <- data.frame(
  Group = group26[ord]
)

rownames(annotation_col2) <- colnames(mat_common2)

png(
  "./Heatmap_16_Common_Genes_Ordered.png",
  width = 2200,
  height = 1600,
  res = 250
)

pheatmap(
  mat_common2,
  annotation_col = annotation_col2,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 14,
  fontsize_col = 8,
  main = "Common DEGs Across Normal, CP and Tumor"
)

dev.off()

common_all <- c(
  "GNMT","TMED11P","CTRL","SYCN","SERPINI2",
  "KLK1","CTRC","CELA2B","AOX1","RBPJL",
  "CELA2A","PAK3","PAIP2B","ERP27",
  "CELA3B","CLPS"
)

common_annot <- annot[
  annot$SYMBOL %in% common_all,
]

head(common_annot)

table(common_annot$SYMBOL)

common_annot <- common_annot[
  !duplicated(common_annot$SYMBOL),
]

common_annot

common_expr <- expr[
  common_annot$PROBEID,
]

dim(common_expr)
rownames(common_expr) <- common_annot$SYMBOL
group_all <- pheno$`disease state:ch1`

table(group_all)

annotation_col <- data.frame(
  Group = group_all
)

rownames(annotation_col) <- colnames(common_expr)

common_expr_scaled <- t(
  scale(
    t(common_expr)
  )
)

library(pheatmap)

pheatmap(
  common_expr_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  main = "16 Common Genes Across Normal, CP and Tumor"
)

png(
  "./Heatmap_16_Common_Genes.png",
  width = 1800,
  height = 1400,
  res = 200
)

pheatmap(
  common_expr_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  main = "16 Common Genes Across Normal, CP and Tumor"
)

dev.off()


library(ggplot2)

# Expression matrix
expr_all <- expr

# Sample groups
group_all <- pheno$`disease state:ch1`

# PCA
pca <- prcomp(
  t(expr_all),
  scale. = TRUE
)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = group_all
)

ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = Group
  )
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCA of GSE143754 Samples"
  )
ggsave(
  "./PCA_AllSamples.png",
  width = 8,
  height = 6,
  dpi = 300
)

sample_dist <- dist(
  t(expr)
)

hc <- hclust(sample_dist)

plot(
  hc,
  main = "Hierarchical Clustering",
  xlab = "",
  sub = ""
)

png(
  "./Hierarchical_Clustering.png",
  width = 1400,
  height = 1000,
  res = 150
)

plot(
  hc,
  main = "Hierarchical Clustering",
  xlab = "",
  sub = ""
)

dev.off()

boxplot(
  expr,
  outline = FALSE,
  las = 2,
  main = "Normalized Expression"
)

png(
  "./Boxplot_Normalized.png",
  width = 1600,
  height = 1000,
  res = 150
)

boxplot(
  expr,
  outline = FALSE,
  las = 2,
  main = "Normalized Expression"
)

dev.off()

common_all

common_rows <- sig_tcp_genes[
  sig_tcp_genes$SYMBOL %in% common_all,
]

common_rows <- common_rows[
  !duplicated(common_rows$SYMBOL),
]

common_rows

mat16 <- expr[
  common_rows$PROBEID,
]

rownames(mat16) <- common_rows$SYMBOL

mat16 <- t(scale(t(mat16)))

ann <- data.frame(
  Group = pheno$`disease state:ch1`
)

rownames(ann) <- colnames(mat16)

library(pheatmap)

pheatmap(
  mat16,
  annotation_col = ann,
  fontsize_row = 10,
  main = "16 Common Genes"
)

png(
  "./Common16_Heatmap.png",
  width = 1800,
  height = 1400,
  res = 150
)

pheatmap(
  mat16,
  annotation_col = ann,
  fontsize_row = 10,
  main = "16 Common Genes"
)

dev.off()

top10 <- head(sig_genes, 10)

ggplot(
  top10,
  aes(
    reorder(SYMBOL, logFC),
    logFC
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top 10 DEGs",
    x = "Gene",
    y = "logFC"
  )
ggsave(
  "./Top10_DEGs.png",
  width = 8,
  height = 6,
  dpi = 300
)

library(clusterProfiler)
library(org.Hs.eg.db)

genes <- common_all

entrez <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_bp <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(go_bp)
barplot(
  go_bp,
  showCategory = 15,
  title = "GO Biological Process"
)

png(
  "./GO_BP.png",
  width=1400,
  height=1000,
  res=150
)

barplot(
  go_bp,
  showCategory=15
)

dev.off()
kegg <- enrichKEGG(
  gene = entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(kegg)
dotplot(
  kegg,
  showCategory = 15
)
png(
  "./KEGG_Pathways.png",
  width=1400,
  height=1000,
  res=150
)

dotplot(
  kegg,
  showCategory=15
)

dev.off()
length(top50_tcp)

sum(top50_tcp %in% rownames(expr_tcp))
mat_tcp <- expr_tcp[top50_tcp, ]

mat_tcp <- t(scale(t(mat_tcp)))

annotation_col <- data.frame(
  Group = group_tcp
)

rownames(annotation_col) <- colnames(mat_tcp)

library(pheatmap)

pheatmap(
  mat_tcp,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  fontsize_col = 10,
  main = "Tumor vs Chronic Pancreatitis"
)

png(
  "./Heatmap_Tumor_vs_CP.png",
  width = 1800,
  height = 1500,
  res = 200
)

pheatmap(
  mat_tcp,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  fontsize_col = 10,
  main = "Tumor vs Chronic Pancreatitis"
)

dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)

length(common_all)
head(common_all)

gene.df <- bitr(
  common_all,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene.df
nrow(gene.df)

ego <- enrichGO(
  gene          = gene.df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.2
)

as.data.frame(ego)

genes <- sig_tcp_genes$SYMBOL

genes <- genes[!is.na(genes)]
genes <- unique(genes)

length(genes)
gene.df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(as.data.frame(ego))
go_results <- as.data.frame(ego)

write.csv(
  go_results,
  "./GO_BP_Tumor_vs_CP.csv",
  row.names = FALSE
)
library(clusterProfiler)

png(
  "./GO_BP_Tumor_vs_CP.png",
  width = 1800,
  height = 1400,
  res = 200
)

dotplot(
  ego,
  showCategory = 20,
  font.size = 12
)

dev.off()
library(clusterProfiler)
library(org.Hs.eg.db)

genes <- sig_genes$SYMBOL

genes <- genes[!is.na(genes)]
genes <- unique(genes)

gene.df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_tn <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(as.data.frame(ego_tn))

png(
  "./GO_Tumor_vs_Normal.png",
  width = 1800,
  height = 1400,
  res = 200
)

dotplot(ego_tn, showCategory = 20)

dev.off()
genes <- sig_cp_genes$SYMBOL

genes <- genes[!is.na(genes)]
genes <- unique(genes)

gene.df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_cp <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(as.data.frame(ego_cp))
png(
  "./GO_CP_vs_Normal.png",
  width = 1800,
  height = 1400,
  res = 200
)

dotplot(ego_cp, showCategory = 20)

dev.off()
library(clusterProfiler)

kegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(as.data.frame(kegg))
options(timeout = 300)

kegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
  
  
)

write.table(
  sig_tcp_genes$SYMBOL,
  "./TCP_gene_symbols.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  common_genes,
  "./Common108_Genes.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

load("./Panc_Workspace.RData")
nrow(sig_genes)
nrow(sig_cp_genes)

gene_tn <- unique(sig_genes$SYMBOL)
gene_tn <- gene_tn[!is.na(gene_tn)]

gene_cp <- unique(sig_cp_genes$SYMBOL)
gene_cp <- gene_cp[!is.na(gene_cp)]

common_genes <- intersect(
  gene_tn,
  gene_cp
)

length(common_genes)

write.csv(
  data.frame(Gene = common_genes),
  "./Common108_Genes.csv",
  row.names = FALSE
)

exists("sig_tcp_genes")
exists("deg_tcp")

write.table(
  common_genes,
  file = "./STRING_108_Genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
head(common_genes)
length(common_genes)
table(pheno$`disease state:ch1`)

group_all <- pheno$`disease state:ch1`

keep_tcp <- group_all %in% c(
  "Chronic Pancreatitis",
  "Tumor"
)

table(keep_tcp)

expr_tcp <- expr[, keep_tcp]

group_tcp <- factor(
  group_all[keep_tcp],
  levels = c(
    "Chronic Pancreatitis",
    "Tumor"
  )
)

table(group_tcp)


length(group_tcp)
ncol(expr_tcp)


library(limma)

design_tcp <- model.matrix(~0 + group_tcp)

colnames(design_tcp) <- c(
  "CP",
  "Tumor"
)

design_tcp

fit_tcp <- lmFit(
  expr_tcp,
  design_tcp
)

cont_tcp <- makeContrasts(
  Tumor - CP,
  levels = design_tcp
)
fit_tcp2 <- contrasts.fit(
  fit_tcp,
  cont_tcp
)

fit_tcp2 <- eBayes(
  fit_tcp2
)

deg_tcp <- topTable(
  fit_tcp2,
  number = Inf,
  adjust.method = "BH"
)

head(deg_tcp)


sig_tcp <- subset(
  deg_tcp,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

nrow(sig_tcp)

table(
  ifelse(
    sig_tcp$logFC > 0,
    "Up",
    "Down"
  )
)
table(group_tcp)

nrow(sig_tcp)

table(ifelse(sig_tcp$logFC>0,"Up","Down"))

library(AnnotationDbi)
library(hta20transcriptcluster.db)

probe_ids <- rownames(sig_tcp)

annot_tcp <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = probe_ids,
  columns = c(
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "PROBEID"
)

sig_tcp$PROBEID <- rownames(sig_tcp)

sig_tcp_annot <- merge(
  sig_tcp,
  annot_tcp,
  by = "PROBEID",
  all.x = TRUE
)

sig_tcp_genes <- sig_tcp_annot[
  !is.na(sig_tcp_annot$SYMBOL),
]

sig_tcp_genes <- sig_tcp_genes[
  order(sig_tcp_genes$adj.P.Val),
]

sig_tcp_genes <- sig_tcp_genes[
  !duplicated(sig_tcp_genes$SYMBOL),
]

nrow(sig_tcp_genes)

write.csv(
  sig_tcp_genes,
  "./Tumor_vs_CP_DEGs.csv",
  row.names = FALSE
)
save.image(
  "./Panc_Workspace_v2.RData"
)



nrow(sig_tcp_genes)

head(
  sig_tcp_genes[
    ,
    c("SYMBOL","logFC","adj.P.Val")
  ],
  20
)

gene_cp <- unique(sig_cp_genes$SYMBOL)
gene_cp <- gene_cp[!is.na(gene_cp)]

gene_tcp <- unique(sig_tcp_genes$SYMBOL)
gene_tcp <- gene_tcp[!is.na(gene_tcp)]

malignant_switch <- setdiff(
  gene_tcp,
  gene_cp
)

length(malignant_switch)

head(malignant_switch, 20)


malignant_switch <- setdiff(
  unique(sig_tcp_genes$SYMBOL),
  unique(sig_cp_genes$SYMBOL)
)

length(malignant_switch)

write.csv(
  data.frame(Gene = malignant_switch),
  "./Malignant_Switch_628_Genes.csv",
  row.names = FALSE
)

writeLines(
  malignant_switch,
  "./Malignant_Switch_628_Genes.txt"
)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


colnames(sig_genes)
colnames(sig_cp_genes)
head(sig_genes[, c("SYMBOL","logFC")])
head(sig_cp_genes[, c("SYMBOL","logFC")])
nrow(sig_genes)
nrow(sig_cp_genes)
malignant_switch <- scan(
  "results/Malignant_Switch_628_Genes.txt",
  what = character()
)

length(malignant_switch)

gene.df <- bitr(
  malignant_switch,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene.df)

ego_bp <- enrichGO(
  gene          = gene.df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(as.data.frame(ego_bp))

write.csv(
  as.data.frame(ego_bp),
  "results/GO_BP_628.csv",
  row.names = FALSE
)

p1 <- dotplot(
  ego_bp,
  showCategory = 15
) +
  ggtitle("GO Biological Process Enrichment")

p1

ggsave(
  "results/GO_BP_Dotplot.png",
  p1,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "results/GO_BP_Dotplot.png",
  p1,
  width = 10,
  height = 7,
  dpi = 300
)

ekegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(as.data.frame(ekegg))

write.csv(
  as.data.frame(ekegg),
  "results/KEGG_628.csv",
  row.names = FALSE
)
p2 <- dotplot(
  ekegg,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment")

p2
ggsave(
  "results/KEGG_Dotplot.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)

options(timeout = 300)

ekegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
head(as.data.frame(ekegg))
write.csv(
  as.data.frame(ekegg),
  "results/KEGG_628.csv",
  row.names = FALSE
)

p2 <- dotplot(
  ekegg,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment")

p2

ggsave(
  "results/KEGG_Dotplot.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)
valid_genes <- scan(
  "results/valid_genes.txt",
  what = character()
)

length(valid_genes)
library(clusterProfiler)
library(org.Hs.eg.db)

gene.df80 <- bitr(
  valid_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene.df80)
options(timeout = 300)

ekegg80 <- enrichKEGG(
  gene = gene.df80$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(as.data.frame(ekegg80))

write.csv(
  as.data.frame(ekegg80),
  "results/KEGG_80_Genes.csv",
  row.names = FALSE
)

library(enrichplot)
library(ggplot2)

p80 <- dotplot(
  ekegg80,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment of 80 Spatially Validated Genes")

p80

ggsave(
  "results/KEGG_80_Genes_Dotplot.png",
  p80,
  width = 10,
  height = 7,
  dpi = 300
)

ego80 <- enrichGO(
  gene = gene.df80$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(
  ego80,
  showCategory = 15
) +
  ggtitle("GO Biological Process Enrichment of 80 Spatially Validated Genes")
pGO80 <- dotplot(
  ego80,
  showCategory = 15
) +
  ggtitle("GO Biological Process Enrichment of 80 Spatially Validated Genes")
ggsave(
  "results/GO_80_Genes_Dotplot.png",
  pGO80,
  width = 10,
  height = 7,
  dpi = 300
)

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

g <- grViz("
# put the workflow code here
")

svg <- export_svg(g)
writeLines(svg, "workflow.svg")

rsvg_png(
  charToRaw(svg),
  file = "Workflow_Figure.png",
  width = 3000,
  height = 4000
)

install.packages("DiagrammeR")

png("GO_80_Genes_Dotplot.png",
    width = 2800,
    height = 2400,
    res = 300)

dotplot(ego_bp,
        showCategory = 15,
        font.size = 12) +
  ggtitle("GO Biological Process Enrichment of 80 Spatially Validated Genes")

dev.off()
library(DiagrammeR)

grViz("
digraph workflow {

graph [layout = dot, rankdir = TB]

node [shape = box,
      style = filled,
      fillcolor = LightBlue,
      fontname = Helvetica]

A [label='GSE143754\nDiscovery Cohort']

B [label='PDAC vs Normal\n401 DEGs']

C [label='CP vs Normal\n365 DEGs']

D [label='PDAC vs CP\n1751 DEGs']

E [label='Malignant-Switch Filtering\n(PDAC vs CP) - (CP vs Normal)']

F [label='628 Malignant-Switch Genes']

G [label='Spatial Validation\nGSE208536']

H [label='80 Spatially Detected Genes']

I [label='ANOVA + Spearman\nProgression Analysis']

J [label='13 Significant\nProgression Genes']

A -> B
A -> C
A -> D

B -> E
C -> E
D -> E

E -> F
F -> G
G -> H
H -> I
I -> J
}
")
ls()
save.image(
  "./Panc_Workspace.RData"
)

