# Recovered effective analysis script
# Source: R_proj\Proj_Panc\.Rproj.user\A1A16783\sources\session-f90c0a1b\BE8B1EC3-contents
# Paths were normalized from the original local D:/R_proj layout to project-relative paths.
# The script is preserved from the working project history and may require the listed R packages
# plus public GEO downloads to rerun fully from scratch.

library(readxl)

ev <- read_excel(
  "data/source/GSE304572_processed_data.xlsx"
)

dim(ev)

head(ev[,1:10])

colnames(ev)[1:20]
tail(colnames(ev),20)

cand <- read.csv(
  "results/Candidate_miRNAs.csv"
)

head(cand)

# Convert miRNA names to match GSE304572 format
cand$miRNA <- gsub("^hsa-", "", cand$ID)

# miRNAs present in EV dataset
ev_mirnas <- ev$EnsgID

# Overlap
overlap <- intersect(
  cand$miRNA,
  ev_mirnas
)

length(overlap)

overlap

write.csv(
  data.frame(miRNA=overlap),
  "results/EV_Validated_miRNAs.csv",
  row.names=FALSE
)

length(overlap)

overlap

library(GEOquery)

gse <- getGEO(
  filename="data/source/GSE304572_series_matrix.txt.gz"
)

pheno <- pData(gse)

dim(pheno)

colnames(pheno)

colnames(pheno)

head(pheno[,1:15])
View(pheno)
grep("disease|group|diagnosis|status|condition",
     colnames(pheno),
     ignore.case = TRUE,
     value = TRUE)

table(grepl("PDAC", pheno$title))
head(pheno$title,20)
tail(pheno$title,20)
unique(sub(",.*","", pheno$title))
table(sub(",.*","", pheno$title))

group <- ifelse(
  grepl("PDAC", pheno$title),
  "PDAC",
  "Benign"
)

table(group)

validated <- overlap

ev17 <- ev[
  ev$EnsgID %in% validated,
]

dim(ev17)

expr <- as.matrix(ev17[, -1])

rownames(expr) <- ev17$EnsgID

expr <- log2(expr + 1)

dim(expr)

colnames(expr)[1:10]

pheno$title[1:10]

library(limma)

group <- factor(
  group,
  levels = c("Benign","PDAC")
)

design <- model.matrix(~group)

design

fit <- lmFit(expr, design)

fit <- eBayes(fit)

deg_miRNA <- topTable(
  fit,
  coef = 2,
  number = Inf
)

head(deg_miRNA)

write.csv(
  deg_miRNA,
  "results/EV_miRNA_DEG.csv"
)

sig_miRNA <- deg_miRNA[
  deg_miRNA$P.Value < 0.05,
]

sig_miRNA

write.csv(
  sig_miRNA,
  "results/EV_significant_miRNAs.csv"
)

library(ggplot2)

volc <- deg_miRNA

volc$miRNA <- rownames(volc)

volc$Significant <- ifelse(
  volc$P.Value < 0.05,
  "Significant",
  "Not Significant"
)

ggplot(volc,
       aes(x=logFC,
           y=-log10(P.Value),
           color=Significant)) +
  
  geom_point(size=4) +
  
  geom_vline(xintercept=c(-0.5,0.5),
             linetype="dashed") +
  
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed") +
  
  geom_text(
    data=subset(volc,P.Value<0.05),
    aes(label=miRNA),
    vjust=-1,
    size=4
  ) +
  
  theme_classic(base_size=14) +
  
  labs(
    title="EV miRNA Validation (PDAC vs Benign)",
    x="log2 Fold Change",
    y="-log10(P value)"
  )
png(
  "results/EV_miRNA_Volcano.png",
  width=2400,
  height=2000,
  res=300
)

ggplot(volc,
       aes(logFC,-log10(P.Value),
           color=Significant)) +
  geom_point(size=4) +
  geom_vline(xintercept=c(-0.5,0.5),
             linetype="dashed") +
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed") +
  geom_text(
    data=subset(volc,P.Value<0.05),
    aes(label=miRNA),
    vjust=-1,
    size=4
  ) +
  theme_classic(base_size=14)

dev.off

library(ggplot2)

box34 <- data.frame(
  Expression = as.numeric(expr["miR-34a-5p",]),
  Group = group
)

ggplot(box34,
       aes(Group,
           Expression,
           fill=Group)) +
  
  geom_boxplot(width=0.6,
               alpha=0.8) +
  
  geom_jitter(width=0.15,
              size=2) +
  
  theme_classic(base_size=14) +
  
  labs(
    title="miR-34a-5p",
    y="log2 Expression"
  )

png(
  "results/miR34a_Boxplot.png",
  width=1800,
  height=1600,
  res=300
)

ggplot(box34,
       aes(Group,
           Expression,
           fill=Group)) +
  geom_boxplot(width=0.6) +
  geom_jitter(width=0.15,size=2) +
  theme_classic(base_size=14)

dev.off()

box107 <- data.frame(
  Expression = as.numeric(expr["miR-107",]),
  Group = group
)

ggplot(box107,
       aes(Group,
           Expression,
           fill=Group)) +
  
  geom_boxplot(width=0.6,
               alpha=0.8) +
  
  geom_jitter(width=0.15,
              size=2) +
  
  theme_classic(base_size=14) +
  
  labs(
    title="miR-107",
    y="log2 Expression"
  )

png(
  "results/miR107_Boxplot.png",
  width=1800,
  height=1600,
  res=300
)

ggplot(box107,
       aes(Group,
           Expression,
           fill=Group)) +
  geom_boxplot(width=0.6) +
  geom_jitter(width=0.15,size=2) +
  theme_classic(base_size=14)

dev.off()
png(
  "results/miR107_Boxplot.png",
  width=1800,
  height=1600,
  res=300
)

ggplot(box107,
       aes(Group,
           Expression,
           fill=Group)) +
  geom_boxplot(width=0.6) +
  geom_jitter(width=0.15,size=2) +
  theme_classic(base_size=14)

dev.off()

# Check for missing values
sum(is.na(expr_sig))
sum(is.nan(as.matrix(expr_sig)))
sum(is.infinite(as.matrix(expr_sig)))
library(pheatmap)

# Significant miRNAs from limma result
sig_names <- rownames(sig_miRNA)

# Extract expression matrix
expr_sig <- expr[sig_names, ]

# Check
dim(expr_sig)
expr_sig

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(expr_sig)
expr_sig_scaled <- t(scale(t(expr_sig)))
expr_sig_scaled[is.na(expr_sig_scaled)] <- 0
expr_sig_scaled[is.na(expr_sig_scaled)] <- 0
pheatmap(
  expr_sig_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "EV-Validated Significant miRNAs"
)

png(
  "results/EV_miRNA_heatmap.png",
  width = 1800,
  height = 1200,
  res = 300
)

pheatmap(
  expr_sig_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "EV-Validated Significant miRNAs"
)

dev.off()
pdf(
  "results/EV_miRNA_heatmap.pdf",
  width = 8,
  height = 6
)

pheatmap(
  expr_sig_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "EV-Validated Significant miRNAs"
)

dev.off()
nrow(deg_miRNA)
deg_miRNA[, c("logFC","P.Value","adj.P.Val")]


pheatmap(
  expr_sig_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 14,
  cellwidth = 10,
  cellheight = 40,
  border_color = NA,
  main = "Differentially Expressed EV miRNAs in PDAC"
)

apply(expr, 1, range)
expr2 <- expr[
  apply(expr,1,var) > 0,
]

dim(expr)
dim(expr2)

fit <- lmFit(expr2, design)
fit <- eBayes(fit)

deg_miRNA <- topTable(
  fit,
  coef=2,
  number=Inf
)

deg_miRNA

wilcox_results <- data.frame()

for(i in 1:nrow(expr2)){
  
  x <- expr2[i, group=="PDAC"]
  y <- expr2[i, group=="Benign"]
  
  p <- wilcox.test(x,y)$p.value
  
  wilcox_results <- rbind(
    wilcox_results,
    data.frame(
      miRNA=rownames(expr2)[i],
      Pvalue=p
    )
  )
}

wilcox_results <- wilcox_results[
  order(wilcox_results$Pvalue),
]

wilcox_results

write.csv(
  wilcox_results,
  "Results/EV_Wilcoxon.csv",
  row.names=FALSE
)

library(EnhancedVolcano)

EnhancedVolcano(
  deg_miRNA,
  lab = rownames(deg_miRNA),
  x = "logFC",
  y = "P.Value",
  
  pCutoff = 0.05,
  FCcutoff = 0.3,
  
  pointSize = 4,
  labSize = 5,
  
  drawConnectors = TRUE,
  widthConnectors = 0.8,
  
  title = "EV miRNA Validation",
  subtitle = "PDAC vs Benign Plasma EVs"
)



library(pheatmap)

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(expr2)

pheatmap(
  expr2,
  scale="row",
  annotation_col=annotation_col,
  show_colnames=FALSE,
  fontsize_row=10
)

pdf(
  "Results/EV_heatmap.pdf",
  width=8,
  height=6
)

pheatmap(
  expr2,
  scale="row",
  annotation_col=annotation_col,
  show_colnames=FALSE
)

dev.off()


png(
  "Results/EV_heatmap.png",
  width=1800,
  height=1400,
  res=300
)

pheatmap(
  expr2,
  scale="row",
  annotation_col=annotation_col,
  show_colnames=FALSE
)

dev.off()

library(ggplot2)
library(ggrepel)

volcano_df <- deg_miRNA

volcano_df$miRNA <- rownames(volcano_df)

volcano_df$Significant <- ifelse(
  volcano_df$P.Value < 0.05 & abs(volcano_df$logFC) > 0.3,
  "Significant",
  "Not Significant"
)

p <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = -log10(P.Value),
    color = Significant
  )
) +
  geom_point(size = 4) +
  geom_text_repel(
    data = subset(volcano_df, P.Value < 0.1),
    aes(label = miRNA),
    size = 5,
    max.overlaps = 100
  ) +
  geom_vline(xintercept = c(-0.3,0.3),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_bw(base_size = 16) +
  labs(
    title = "EV miRNA Differential Expression",
    subtitle = "PDAC vs Benign",
    x = "Log2 Fold Change",
    y = "-log10(P-value)"
  )

p
ggsave(
  "Results/EV_Volcano.png",
  p,
  width=8,
  height=6,
  dpi=600
)

ggsave(
  "Results/EV_Volcano.pdf",
  p,
  width=8,
  height=6
)


library(ggplot2)

box_df <- data.frame(
  Expression = as.numeric(expr2["miR-107", ]),
  Group = group
)

p1 <- ggplot(
  box_df,
  aes(Group, Expression, fill=Group)
) +
  geom_boxplot(width=0.6) +
  geom_jitter(width=0.15,size=2) +
  theme_bw(base_size=16) +
  ggtitle("miR-107")

p1


ggsave(
  "Results/Boxplot_miR107.png",
  p1,
  width=5,
  height=5,
  dpi=600
)

ggsave(
  "Results/Boxplot_miR107.pdf",
  p1,
  width=5,
  height=5
)


box_df <- data.frame(
  Expression = as.numeric(expr2["miR-34a-5p", ]),
  Group = group
)

p2 <- ggplot(
  box_df,
  aes(Group, Expression, fill=Group)
) +
  geom_boxplot(width=0.6) +
  geom_jitter(width=0.15,size=2) +
  theme_bw(base_size=16) +
  ggtitle("miR-34a-5p")

p2


ggsave(
  "Results/Boxplot_miR34a.png",
  p2,
  width=5,
  height=5,
  dpi=600
)

ggsave(
  "Results/Boxplot_miR34a.pdf",
  p2,
  width=5,
  height=5
)

box_df <- data.frame(
  Expression = as.numeric(expr2["miR-15a-5p", ]),
  Group = group
)

p3 <- ggplot(
  box_df,
  aes(Group, Expression, fill=Group)
) +
  geom_boxplot(width=0.6) +
  geom_jitter(width=0.15,size=2) +
  theme_bw(base_size=16) +
  ggtitle("miR-15a-5p")

p3


ggsave(
  "Results/Boxplot_miR15a.png",
  p3,
  width=5,
  height=5,
  dpi=600
)

ggsave(
  "Results/Boxplot_miR15a.pdf",
  p3,
  width=5,
  height=5
)

library(pROC)

roc107 <- roc(
  response = group,
  predictor = as.numeric(expr2["miR-107", ]),
  levels = c("Benign","PDAC")
)

auc(roc107)

plot(
  roc107,
  print.auc = TRUE,
  main = "ROC - miR-107"
)


png(
  "Results/ROC_miR107.png",
  width=1800,
  height=1800,
  res=300
)

plot(
  roc107,
  print.auc=TRUE,
  main="ROC - miR-107"
)

dev.off()

panel_df <- data.frame(
  miR107 = as.numeric(expr2["miR-107", ]),
  miR34a = as.numeric(expr2["miR-34a-5p", ]),
  miR15a = as.numeric(expr2["miR-15a-5p", ]),
  Group = group
)

model <- glm(
  Group ~ miR107 + miR34a + miR15a,
  data = panel_df,
  family = binomial
)

prob <- predict(
  model,
  type="response"
)

roc_panel <- roc(
  panel_df$Group,
  prob,
  levels=c("Benign","PDAC")
)

auc(roc_panel)

plot(
  roc_panel,
  print.auc=TRUE,
  main="Combined EV miRNA Panel"
)


panel_df <- data.frame(
  t(expr2)
)

panel_df$Group <- group

model <- glm(
  Group ~ .,
  data = panel_df,
  family = binomial
)

prob <- predict(
  model,
  type="response"
)

library(pROC)

roc_all <- roc(
  panel_df$Group,
  prob,
  levels=c("Benign","PDAC")
)

auc(roc_all)
library(pROC)

# Example:
# score = combined biomarker score
# group = factor with Benign / PDAC

roc_obj <- roc(
  response = group,
  predictor = score,
  levels = c("Benign","PDAC"),
  direction = "<"
)

auc(roc_obj)
png(
  "Results/ROC_CombinedPanel.png",
  width = 2000,
  height = 1600,
  res = 300
)

plot(
  roc_obj,
  col = "red",
  lwd = 4,
  main = paste(
    "Combined EV miRNA Panel\nAUC =",
    round(auc(roc_obj),3)
  )
)

abline(a=0,b=1,lty=2,col="gray")

dev.off()
class(score)
str(score)
head(score)

sig4 <- c(
  "miR-107",
  "miR-15a-5p",
  "miR-34a-5p",
  "miR-1908-5p"
)

score_vec <- colMeans(expr2[sig4, ])

class(score_vec)
head(score_vec)


library(pROC)

roc4 <- roc(
  response = group,
  predictor = score_vec,
  levels = c("Benign","PDAC"),
  direction = "<"
)

auc(roc4)

roc4 <- roc(
  response = group,
  predictor = score_vec,
  levels = c("Benign","PDAC"),
  direction = ">"
)

auc(roc4)

tapply(score_vec, group, mean)

panel_df <- data.frame(
  t(expr2)
)

panel_df$Group <- group

model <- glm(
  Group ~ .,
  data = panel_df,
  family = binomial
)

prob <- predict(
  model,
  type="response"
)

library(pROC)

roc_all <- roc(
  panel_df$Group,
  prob,
  levels=c("Benign","PDAC")
)

auc(roc_all)

png(
  "Results/LogisticROC.png",
  width=2000,
  height=1800,
  res=300
)

plot(
  roc_all,
  col="blue",
  lwd=4,
  main=paste0(
    "11-miRNA Logistic Model (AUC = ",
    round(as.numeric(auc(roc_all)),3),
    ")"
  )
)

abline(a=0,b=1,lty=2)

dev.off()


top4 <- c(
  "miR-107",
  "miR-15a-5p",
  "miR-34a-5p",
  "miR-1908-5p"
)

panel4 <- data.frame(
  t(expr2[top4, ])
)

panel4$Group <- group

model4 <- glm(
  Group ~ .,
  data = panel4,
  family = binomial
)

prob4 <- predict(
  model4,
  type = "response"
)

roc4 <- roc(
  panel4$Group,
  prob4,
  levels = c("Benign","PDAC")
)

auc(roc4)

# Install if needed
install.packages("randomForest")

# Load package
library(randomForest)
library(pROC)

# Build dataframe
rf_df <- data.frame(
  t(expr2)
)

rf_df$Group <- group

# Check
dim(rf_df)
head(rf_df)

# Train RF
rf <- randomForest(
  Group ~ .,
  data = rf_df,
  ntree = 2000,
  importance = TRUE
)

# Predict probabilities
prob <- predict(
  rf,
  type = "prob"
)[,"PDAC"]

# ROC
roc_rf <- roc(
  rf_df$Group,
  prob,
  levels = c("Benign","PDAC")
)

auc(roc_rf)


panel3 <- data.frame(
  miR107 = as.numeric(expr2["miR-107",]),
  miR34a = as.numeric(expr2["miR-34a-5p",]),
  miR15a = as.numeric(expr2["miR-15a-5p",]),
  Group = group
)

model3 <- glm(
  Group ~ .,
  data = panel3,
  family = binomial
)

prob3 <- predict(
  model3,
  type="response"
)

library(pROC)

roc3 <- roc(
  panel3$Group,
  prob3,
  levels=c("Benign","PDAC")
)

auc(roc3)

library(pROC)

top3 <- c(
  "miR-107",
  "miR-15a-5p",
  "miR-34a-5p"
)

panel_df <- data.frame(
  t(expr2[top3, ])
)

panel_df$Group <- group

model3 <- glm(
  Group ~ .,
  data = panel_df,
  family = binomial
)

prob3 <- predict(
  model3,
  type = "response"
)

roc3 <- roc(
  panel_df$Group,
  prob3,
  levels = c("Benign","PDAC")
)

auc(roc3)
png(
  "Results/ROC_Top3_miRNA.png",
  width = 2000,
  height = 1800,
  res = 300
)

plot(
  roc3,
  col = "red",
  lwd = 4,
  main = "ROC Curve: Top 3 EV miRNAs"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

text(
  0.6,
  0.2,
  paste(
    "AUC =",
    round(auc(roc3),3)
  ),
  cex = 1.5
)

dev.off()
pdf(
  "Results/ROC_Top3_miRNA.pdf",
  width = 8,
  height = 7
)

plot(
  roc3,
  col = "red",
  lwd = 4,
  main = "ROC Curve: Top 3 EV miRNAs"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

text(
  0.6,
  0.2,
  paste(
    "AUC =",
    round(auc(roc3),3)
  ),
  cex = 1.2
)

dev.off()


library(GEOquery)

gse <- getGEO("GSE259327")

pheno <- pData(gse[[1]])

dim(pheno)

head(pheno)


library(GEOquery)

gse <- getGEO(
  "GSE259327",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

gse <- getGEO(
  "GSE259327",
  GSEMatrix = TRUE,
  AnnotGPL = FALSE,
  getGPL = FALSE
)

library(readxl)

gse259327 <- read_excel(
  "data/source/GSE259327_Raw_counts_PDAC_203samples_new.xlsx"
)

dim(gse259327)

head(gse259327)
colnames(gse259327)[1:20]

tail(colnames(gse259327),20)


sample_names <- as.character(unlist(gse259327[1,-1]))

tail(sample_names,50)

gse259327[20:60,1]
tail(sample_names,50)
print(
  gse259327[21:150,1],
  n = 130
)

library(dplyr)

expr259327 <- gse259327[-c(1,2), ]

rownames(expr259327) <- expr259327[[1]]

expr259327 <- expr259327[,-1]

expr259327 <- as.matrix(expr259327)

mode(expr259327) <- "numeric"

dim(expr259327)
sample_names <- colnames(expr259327)

group259327 <- ifelse(
  grepl("^CTRL", sample_names),
  "Control",
  "PDAC"
)

table(group259327)

validated
head(rownames(expr259327),50)

dim(gse259327)
head(gse259327[,1],20)
rownames(expr259327) <- gse259327[[1]]




expr259327 <- gse259327[-c(1,2), ]

mirna_names <- expr259327[[1]]

expr259327 <- expr259327[,-1]

expr259327 <- as.matrix(expr259327)

mode(expr259327) <- "numeric"

rownames(expr259327) <- mirna_names

dim(expr259327)

head(rownames(expr259327),20)


head(colnames(expr259327),20)



sample_names <- as.character(
  unlist(gse259327[1,-1])
)

head(sample_names,20)

tail(sample_names,20)


expr259327 <- gse259327[-c(1,2), ]

mirna_names <- expr259327[[1]]

expr259327 <- expr259327[,-1]

colnames(expr259327) <- sample_names

expr259327 <- as.matrix(expr259327)

mode(expr259327) <- "numeric"

rownames(expr259327) <- mirna_names


head(colnames(expr259327),20)

tail(colnames(expr259327),20)


group259327 <- ifelse(
  grepl("^CTRL", colnames(expr259327)),
  "Control",
  "PDAC"
)

table(group259327)


present <- intersect(
  validated,
  rownames(expr259327)
)

length(present)

present


library(limma)

expr_val <- expr259327[
  present,
]

expr_val <- log2(expr_val + 1)

group259327 <- factor(
  group259327,
  levels=c("Control","PDAC")
)

design <- model.matrix(~group259327)

fit <- lmFit(
  expr_val,
  design
)

fit <- eBayes(fit)

deg_val <- topTable(
  fit,
  coef=2,
  number=Inf
)

deg_val


write.csv(
  deg_val,
  "Results/GSE259327_Validation_DEG.csv"
)


library(EnhancedVolcano)

png(
  "Results/GSE259327_Validation_Volcano.png",
  width=2600,
  height=2000,
  res=300
)

EnhancedVolcano(
  deg_val,
  lab=rownames(deg_val),
  
  x="logFC",
  y="P.Value",
  
  pCutoff=0.05,
  FCcutoff=0.25,
  
  pointSize=5,
  labSize=6,
  
  drawConnectors=TRUE,
  
  title="External Validation Cohort",
  subtitle="GSE259327 (121 PDAC vs 82 Controls)"
)

dev.off()

library(pheatmap)

sig <- rownames(
  deg_val[
    deg_val$P.Value < 0.05,
  ]
)

expr_sig <- expr_val[sig,]

ann <- data.frame(
  Group=group259327
)

rownames(ann) <- colnames(expr_sig)

pheatmap(
  expr_sig,
  scale="row",
  annotation_col=ann,
  show_colnames=FALSE,
  fontsize_row=10,
  filename="Results/GSE259327_Validation_Heatmap.pdf"
)

deg_val[,c("logFC","P.Value","adj.P.Val")]

top_panel <- c(
  "miR-184",
  "miR-107",
  "miR-216b-5p",
  "miR-20a-5p",
  "miR-134-3p",
  "miR-215-3p"
)

panel_expr <- t(
  expr_val[top_panel,]
)

panel_df <- data.frame(panel_expr)

panel_df$Group <- group259327

model <- glm(
  Group ~ .,
  data=panel_df,
  family=binomial
)

prob <- predict(
  model,
  type="response"
)

library(pROC)

roc_panel <- roc(
  panel_df$Group,
  prob,
  levels=c("Control","PDAC")
)

auc(roc_panel)

png(
  "Results/Final_Biomarker_Panel_ROC.png",
  width=2200,
  height=1800,
  res=300
)

plot(
  roc_panel,
  col="red",
  lwd=4,
  main=paste(
    "Final Biomarker Panel\nAUC =",
    round(as.numeric(auc(roc_panel)),3)
  )
)

abline(
  a=0,
  b=1,
  lty=2
)

dev.off()

coef(summary(model))
abs(coef(model))


library(caret)
library(pROC)

set.seed(123)

ctrl <- trainControl(
  method="cv",
  number=10,
  classProbs=TRUE,
  summaryFunction=twoClassSummary,
  savePredictions="final"
)

panel_df$Group <- factor(
  panel_df$Group,
  levels=c("Control","PDAC")
)

cv_model <- train(
  Group ~ .,
  data=panel_df,
  method="glm",
  family="binomial",
  metric="ROC",
  trControl=ctrl
)

cv_model
roc_cv <- roc(
  cv_model$pred$obs,
  cv_model$pred$PDAC,
  levels=c("Control","PDAC")
)

auc_value <- auc(roc_cv)

auc_value <- auc(roc_cv)

plot(roc_cv, col="#2C7BB6", lwd=3,
     main="10-Fold CV ROC Curve")

abline(a=0, b=1, lty=2, col="gray")

legend("bottomright",
       legend=paste0("AUC = ", round(auc_value, 3)),
       bty="n", text.col="#2C7BB6")

png("ROC_CV_with_AUC.png", width=1200, height=900, res=150)

auc_value <- auc(roc_cv)

plot(roc_cv, col="#2C7BB6", lwd=3,
     main="10-Fold CV ROC Curve")

abline(a=0, b=1, lty=2, col="gray")

legend("bottomright",
       legend=paste0("AUC = ", round(auc_value, 3)),
       bty="n", text.col="#2C7BB6")

dev.off()
library(pROC)

opt <- coords(
  roc_panel,
  x = "best",
  best.method = "youden",
  ret = "threshold"
)

opt


# predicted probability
prob <- predict(model, type = "response")

# choose cutoff (default = 0.5)
pred_class <- ifelse(prob > 0.5, "PDAC", "Control")

pred_class <- factor(pred_class, levels = c("Control", "PDAC"))
panel_df$Group <- factor(panel_df$Group, levels = c("Control", "PDAC"))

library(caret)

conf_matrix <- confusionMatrix(
  pred_class,
  panel_df$Group,
  positive = "PDAC"
)

conf_matrix

sensitivity <- conf_matrix$byClass["Sensitivity"]
specificity <- conf_matrix$byClass["Specificity"]

sensitivity
specificity


library(pROC)

roc_panel <- roc(
  panel_df$Group,
  prob,
  levels = c("Control","PDAC")
)

auc_val <- auc(roc_panel)

png("Results/ROC_Final_Panel.png", width=2000, height=1800, res=300)

plot(
  roc_panel,
  col = "red",
  lwd = 4,
  main = paste0("EV miRNA Panel ROC\nAUC = ", round(auc_val, 3))
)

abline(a=0, b=1, lty=2, col="gray")

legend(
  "bottomright",
  legend = paste0("AUC = ", round(auc_val, 3)),
  bty = "n",
  text.col = "red"
)

dev.off()

library(ggplot2)

cm <- as.data.frame(conf_matrix$table)

ggplot(cm, aes(x=Reference, y=Prediction, fill=Freq)) +
  geom_tile(color="white") +
  geom_text(aes(label=Freq), size=6) +
  scale_fill_gradient(low="white", high="steelblue") +
  theme_minimal(base_size=14) +
  labs(title="Confusion Matrix - EV miRNA Panel")

ggsave(
  "Results/ConfusionMatrix.png",
  width=6,
  height=5,
  dpi=300
)

auc_cv <- auc(roc_cv)

png("Results/CV_ROC.png", width=2000, height=1800, res=300)

plot(
  roc_cv,
  col="#2C7BB6",
  lwd=4,
  main=paste0("10-Fold Cross Validation ROC\nAUC = ", round(auc_cv, 3))
)

abline(a=0, b=1, lty=2, col="gray")

legend(
  "bottomright",
  legend = paste0("AUC = ", round(auc_cv, 3)),
  bty="n",
  text.col="#2C7BB6"
)

dev.off()

library(ggplot2)

metrics <- data.frame(
  Metric = c("Sensitivity", "Specificity"),
  Value = c(
    conf_matrix$byClass["Sensitivity"],
    conf_matrix$byClass["Specificity"]
  )
)

ggplot(metrics, aes(x=Metric, y=Value, fill=Metric)) +
  geom_bar(stat="identity", width=0.6) +
  ylim(0,1) +
  theme_minimal(base_size=14) +
  labs(title="Model Performance Metrics")

ggsave(
  "Results/Sensitivity_Specificity.png",
  width=5,
  height=5,
  dpi=300
)

library(pheatmap)

pheatmap(
  expr_val[top_panel, ],
  scale = "row",
  annotation_col = data.frame(Group = group259327),
  show_colnames = FALSE,
  fontsize_row = 10,
  main = "Final EV miRNA Panel"
)

png("Results/Final_Heatmap.png", width=2000, height=1400, res=300)

pheatmap(
  expr_val[top_panel, ],
  scale = "row",
  annotation_col = data.frame(Group = group259327),
  show_colnames = FALSE,
  fontsize_row = 10,
  main = "Final EV miRNA Panel"
)

dev.off()
library(ggplot2)

miRNAs <- c("miR-184","miR-107","miR-216b-5p",
            "miR-20a-5p","miR-134-3p","miR-215-3p")

for(m in miRNAs){
  
  df <- data.frame(
    Expression = as.numeric(expr_val[m, ]),
    Group = group259327
  )
  
  p <- ggplot(df, aes(Group, Expression, fill=Group)) +
    geom_boxplot(width=0.6) +
    geom_jitter(width=0.15, size=2) +
    theme_classic(base_size=14) +
    ggtitle(m)
  
  ggsave(
    paste0("Results/Boxplot_", m, ".png"),
    p,
    width=5,
    height=5,
    dpi=300
  )
}


library(ggplot2)

# ----------------------------
# 1. Extract metrics
# ----------------------------
auc_val <- as.numeric(auc(roc_panel))

sens <- conf_matrix$byClass["Sensitivity"]
spec <- conf_matrix$byClass["Specificity"]
acc  <- conf_matrix$overall["Accuracy"]

metrics_df <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity"),
  Value  = c(auc_val, sens, spec)
)

# ----------------------------
# 2. Create summary plot
# ----------------------------
p <- ggplot(metrics_df, aes(x = Metric, y = Value, fill = Metric)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_text(
    aes(label = round(Value, 3)),
    vjust = -0.4,
    size = 5
  ) +
  
  ylim(0, 1) +
  
  theme_classic(base_size = 16) +
  
  labs(
    title = "Final EV-miRNA Biomarker Panel Performance",
    subtitle = "External Validation Cohort (GSE259327)",
    y = "Performance Metric (0â€“1)",
    x = NULL
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# ----------------------------
# 3. Save PNG (high resolution)
# ----------------------------
ggsave(
  "Results/Final_Biomarker_Summary.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 600
)

# ----------------------------
# 4. Save PDF (publication quality)
# ----------------------------
ggsave(
  "Results/Final_Biomarker_Summary.pdf",
  plot = p,
  width = 7,
  height = 5
)

# ----------------------------
# 5. Print plot
# ----------------------------
p



library(pROC)
library(ggplot2)
library(gridExtra)
library(grid)

# ---------------------------
# 1. Metrics
# ---------------------------
auc_val <- as.numeric(auc(roc_panel))

sens <- conf_matrix$byClass["Sensitivity"]
spec <- conf_matrix$byClass["Specificity"]

# ---------------------------
# 2. ROC base plot (ggplot-style)
# ---------------------------
roc_df <- data.frame(
  TPR = roc_panel$sensitivities,
  FPR = 1 - roc_panel$specificities
)

roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  
  geom_line(color = "red", size = 1.5) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Final EV-miRNA Diagnostic Model",
    subtitle = paste0("AUC = ", round(auc_val, 3)),
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  
  annotate(
    "text",
    x = 0.65,
    y = 0.25,
    label = paste0(
      "Sensitivity = ", round(sens, 3), "\n",
      "Specificity = ", round(spec, 3)
    ),
    size = 5,
    hjust = 0
  )

# ---------------------------
# 3. Metrics bar plot
# ---------------------------
metrics_df <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity"),
  Value = c(auc_val, sens, spec)
)

metrics_plot <- ggplot(metrics_df, aes(x = Metric, y = Value, fill = Metric)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 5) +
  
  ylim(0, 1) +
  
  theme_classic(base_size = 14) +
  
  theme(legend.position = "none") +
  
  labs(
    title = "Performance Summary",
    y = "Score"
  )

# ---------------------------
# 4. Combine both plots
# ---------------------------
final_fig <- grid.arrange(
  roc_plot,
  metrics_plot,
  ncol = 2
)

# ---------------------------
# 5. Save figure
# ---------------------------
ggsave(
  "Results/Final_ML_Figure.png",
  final_fig,
  width = 12,
  height = 5,
  dpi = 600
)

ggsave(
  "Results/Final_ML_Figure.pdf",
  final_fig,
  width = 12,
  height = 5
)
install.packages("gridExtra")

library(ggplot2)
library(gridExtra)

# ---------------------------
# miR-107
# ---------------------------
df1 <- data.frame(
  Expression = as.numeric(expr2["miR-107", ]),
  Group = group
)

p1 <- ggplot(df1, aes(Group, Expression, fill = Group)) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic(base_size = 14) +
  ggtitle("miR-107") +
  theme(legend.position = "none")

# ---------------------------
# miR-34a-5p
# ---------------------------
df2 <- data.frame(
  Expression = as.numeric(expr2["miR-34a-5p", ]),
  Group = group
)

p2 <- ggplot(df2, aes(Group, Expression, fill = Group)) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic(base_size = 14) +
  ggtitle("miR-34a-5p") +
  theme(legend.position = "none")

# ---------------------------
# miR-15a-5p
# ---------------------------
df3 <- data.frame(
  Expression = as.numeric(expr2["miR-15a-5p", ]),
  Group = group
)

p3 <- ggplot(df3, aes(Group, Expression, fill = Group)) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic(base_size = 14) +
  ggtitle("miR-15a-5p") +
  theme(legend.position = "none")

# ---------------------------
# Combine into one figure
# ---------------------------
final_boxplot <- grid.arrange(p1, p2, p3, ncol = 3)

# ---------------------------
# Save PNG
# ---------------------------
ggsave(
  "Results/EV_miRNA_All_Boxplots.png",
  final_boxplot,
  width = 12,
  height = 4,
  dpi = 600
)

# ---------------------------
# Save PDF
# ---------------------------
ggsave(
  "Results/EV_miRNA_All_Boxplots.pdf",
  final_boxplot,
  width = 12,
  height = 4
)

view(top_panel)

model <- glm(Group ~ ., data=panel_df, family=binomial)
coef_table <- summary(model)$coefficients

# remove intercept
coef_table <- coef_table[-1, ]

# sort by absolute effect size
top6 <- head(
  rownames(coef_table)[order(abs(coef_table[,1]), decreasing = TRUE)],
  6
)

top6

getwd()

list.files()

list.files(pattern = "268771")
library(data.table)

gse268771 <- fread(
  "data/source/GSE268771_miRNA_expression_RPM.txt.gz",
  data.table = FALSE
)

gse268771 <- read.delim(
  gzfile("data/source/GSE268771_miRNA_expression_RPM.txt.gz"),
  check.names = FALSE
)


dim(gse268771)

head(gse268771[,1:10])

colnames(gse268771)[1:20]

head(gse268771[,1])

series <- readLines(
  gzfile("data/source/GSE268771_series_matrix.txt.gz")
)

head(series,50)

grep("characteristics", series, value=TRUE)

grep("title", series, value=TRUE)

grep("source_name", series, value=TRUE)

panel6 <- c(
  "miR-184",
  "miR-107",
  "miR-216b-5p",
  "miR-20a-5p",
  "miR-134-3p",
  "miR-215-3p"
)
grep(
  "miR-184|miR-107|miR-216b|miR-20a|miR-134|miR-215",
  gse268771$miRNA,
  value=TRUE
)

coef(model)
prob <- plogis(score)

validated6 <- c(
  "hsa-miR-184",
  "hsa-miR-107",
  "hsa-miR-216b-5p",
  "hsa-miR-20a-5p",
  "hsa-miR-134-3p",
  "hsa-miR-215-5p"
)

val_panel <- gse268771[
  match(validated6, gse268771$miRNA),
]

val_panel

val_panel[,1]

expr_val <- val_panel[,-1]

expr_val <- as.data.frame(t(expr_val))

colnames(expr_val) <- c(
  "miR.184",
  "miR.107",
  "miR.216b.5p",
  "miR.20a.5p",
  "miR.134.3p",
  "miR.215.3p"
)

expr_val <- data.frame(
  lapply(expr_val, as.numeric)
)

head(expr_val)

colnames(gse268771)[1:67]

readLines(
  gzfile("data/source/GSE268771_series_matrix.txt.gz"),
  n = 100
)

score <- -11.8050305 +
  (-0.6064156)*expr_val$miR.184 +
  2.6258975 *expr_val$miR.107 +
  (-0.5090217)*expr_val$miR.216b.5p +
  (-1.1363684)*expr_val$miR.20a.5p +
  (-0.4045020)*expr_val$miR.134.3p +
  0.6704752 *expr_val$miR.215.3p

prob <- plogis(score)

summary(prob)

x <- readLines(
  gzfile("data/source/GSE268771_series_matrix.txt.gz")
)

grep("PDAC", x, value=TRUE)
grep("normal", x, value=TRUE)
grep("characteristics", x, value=TRUE)
prob <- plogis(score)

summary(prob)

cat(
  x[grep("characteristics", x)],
  sep="\n"
)
group268771 <- c(
  rep("Control",15),
  rep("PDAC",51)
)

library(pROC)

roc_external <- roc(
  group268771,
  prob,
  levels=c("Control","PDAC")
)

auc(roc_external)

plot(
  roc_external,
  col="blue",
  lwd=4,
  main=paste(
    "External Validation ROC\nAUC =",
    round(as.numeric(auc(roc_external)),3)
  )
)
abline(a=0,b=1,lty=2)

expr_val <- log2(expr_val + 1)
expr_val <- as.data.frame(scale(expr_val))

score <- predict(
  model,
  newdata = expr_val,
  type = "link"
)

prob <- plogis(score)

prob <- predict(
  model,
  newdata = expr_val,
  type = "response"
)
length(group268771)
ncol(expr_val)
roc_external <- roc(
  group268771,
  prob,
  levels = c("Control","PDAC")
)

x <- readLines(gzfile("data/source/GSE268771_series_matrix.txt.gz"))

grep("source_name", x, value = TRUE)
grep("disease", x, value = TRUE)
grep("PDAC", x, value = TRUE)
grep("normal", x, value = TRUE)

library(limma)

auc(roc_external)

miR-107
miR-184
miR-20a-5p
miR-134-3p
miR-216b-5p



title_line <- grep(
  "!Sample_title",
  x,
  value = TRUE
)

cat(title_line)



group268771 <- c(
  rep("PDAC",51),
  rep("Control",15)
)
candidate <- c(
  "miR-15a-5p",
  "miR-20a-5p",
  "miR-21-5p",
  "miR-107",
  "miR-34a-5p",
  "miR-182-5p",
  "miR-424-5p",
  "miR-184",
  "miR-630",
  "miR-216b-5p",
  "miR-509-3-5p",
  "miR-1908-5p",
  "miR-1306-5p",
  "miR-215-3p",
  "miR-134-3p"
)
grep(
  paste(candidate, collapse="|"),
  gse268771$miRNA,
  value=TRUE
)

design <- model.matrix(~group268771)

fit <- lmFit(expr268771, design)
fit <- eBayes(fit)

topTable(
  fit,
  coef=2,
  number=Inf
)


expr268771 <- as.matrix(gse268771[, -1])

mode(expr268771) <- "numeric"

rownames(expr268771) <- gse268771$miRNA


dim(expr268771)

head(rownames(expr268771))

group268771 <- c(
  rep("PDAC", 51),
  rep("Control", 15)
)

table(group268771)

candidate <- c(
  "hsa-miR-15a-5p",
  "hsa-miR-20a-5p",
  "hsa-miR-21-5p",
  "hsa-miR-107",
  "hsa-miR-34a-5p",
  "hsa-miR-182-5p",
  "hsa-miR-424-5p",
  "hsa-miR-184",
  "hsa-miR-630",
  "hsa-miR-216b-5p",
  "hsa-miR-1908-5p",
  "hsa-miR-1306-5p",
  "hsa-miR-134-3p"
)

present <- intersect(candidate, rownames(expr268771))

present

library(limma)

expr_sub <- expr268771[present, ]

design <- model.matrix(~ group268771)

fit <- lmFit(expr_sub, design)

fit <- eBayes(fit)

deg268771 <- topTable(
  fit,
  coef = 2,
  number = Inf
)

deg268771

deg268771[
  order(deg268771$adj.P.Val),
  c("logFC","P.Value","adj.P.Val")
]

prob <- predict(
  model,
  newdata = expr_val,
  type = "response"
)

roc_external <- roc(
  group268771,
  -prob,
  levels = c("Control","PDAC")
)

auc(roc_external)

candidate
rownames(expr268771)

common_mirs <- c(
  "hsa-miR-15a-5p",
  "hsa-miR-20a-5p",
  "hsa-miR-21-5p",
  "hsa-miR-107",
  "hsa-miR-34a-5p",
  "hsa-miR-182-5p",
  "hsa-miR-424-5p",
  "hsa-miR-184",
  "hsa-miR-216b-5p",
  "hsa-miR-1908-5p",
  "hsa-miR-1306-5p",
  "hsa-miR-134-3p"
)

expr268771_log <- log2(expr268771 + 1)
expr268771_scaled <- t(
  scale(
    t(expr268771_log)
  )
)

model_scaled <- glm(
  Group ~ .,
  data=train_df,
  family=binomial
)

table(group268771)

candidate <- c(
  "hsa-miR-15a-5p",
  "hsa-miR-20a-5p",
  "hsa-miR-21-5p",
  "hsa-miR-107",
  "hsa-miR-34a-5p",
  "hsa-miR-182-5p",
  "hsa-miR-424-5p",
  "hsa-miR-184",
  "hsa-miR-216b-5p",
  "hsa-miR-1908-5p",
  "hsa-miR-1306-5p",
  "hsa-miR-134-3p"
)

validation_summary <- data.frame(
  miRNA = rownames(deg268771),
  logFC = deg268771$logFC,
  P = deg268771$P.Value,
  FDR = deg268771$adj.P.Val
)

validation_summary

sum(sign(fc259327) == sign(fc268771))


# GSE259327
fc259327 <- deg_val$logFC
names(fc259327) <- rownames(deg_val)

# GSE268771
fc268771 <- deg268771$logFC
names(fc268771) <- rownames(deg268771)


names(fc268771) <- sub("^hsa-", "", names(fc268771))

common_mirs <- intersect(
  names(fc259327),
  names(fc268771)
)

common_mirs

sum(
  sign(fc259327[common_mirs]) ==
    sign(fc268771[common_mirs])
)


direction_table <- data.frame(
  miRNA = common_mirs,
  logFC_GSE259327 = fc259327[common_mirs],
  logFC_GSE268771 = fc268771[common_mirs]
)

direction_table$SameDirection <-
  sign(direction_table$logFC_GSE259327) ==
  sign(direction_table$logFC_GSE268771)

direction_table


library(pROC)

for(m in present){
  
  roc1 <- roc(
    group268771,
    expr268771[m,],
    levels=c("Control","PDAC")
  )
  
  cat(
    m,
    "AUC =",
    round(as.numeric(auc(roc1)),3),
    "\n"
  )
}

for(m in present){
  
  roc1 <- roc(
    group268771,
    -expr268771[m,],
    levels=c("Control","PDAC")
  )
  
  cat(
    m,
    "Reverse AUC =",
    round(as.numeric(auc(roc1)),3),
    "\n"
  )
}


panel <- c(
  "hsa-miR-184",
  "hsa-miR-107",
  "hsa-miR-216b-5p",
  "hsa-miR-20a-5p",
  "hsa-miR-134-3p",
  "hsa-miR-215-5p"
)

X <- as.data.frame(
  t(expr268771[panel,])
)

X$Group <- factor(group268771)

model268 <- glm(
  Group ~ .,
  data=X,
  family=binomial
)

prob268 <- predict(
  model268,
  type="response"
)

roc268 <- roc(
  X$Group,
  prob268,
  levels=c("Control","PDAC")
)

auc(roc268)


panel12 <- c(
  "hsa-miR-15a-5p",
  "hsa-miR-20a-5p",
  "hsa-miR-21-5p",
  "hsa-miR-107",
  "hsa-miR-34a-5p",
  "hsa-miR-182-5p",
  "hsa-miR-424-5p",
  "hsa-miR-184",
  "hsa-miR-216b-5p",
  "hsa-miR-1908-5p",
  "hsa-miR-1306-5p",
  "hsa-miR-134-3p"
)

X12 <- as.data.frame(t(expr268771[panel12,]))
X12$Group <- factor(group268771)

model12 <- glm(
  Group ~ .,
  data = X12,
  family = binomial
)

prob12 <- predict(model12, type="response")

roc12 <- pROC::roc(
  X12$Group,
  prob12,
  levels=c("Control","PDAC")
)

auc(roc12)
panel6 <- c(
  "miR.184",
  "miR.107",
  "miR.216b.5p",
  "miR.20a.5p",
  "miR.134.3p",
  "miR.215.3p"
)

expr_data
group


save.image(
  "./Panc_Workspace.RData"
)
expr259327
group259327

train_df

panel <- c(
  "hsa-miR-21-5p",
  "hsa-miR-424-5p",
  "hsa-miR-15a-5p",
  "hsa-miR-107",
  "hsa-miR-20a-5p",
  "hsa-miR-1306-5p"
)

panel %in% rownames(expr268771)

Is()
is()
ls()
library(reshape2)
library(ggplot2)

box_df <- data.frame(
  t(expr259327[panel, ])
)

box_df$Group <- group259327

long_df <- melt(
  box_df,
  id.vars = "Group",
  variable.name = "miRNA",
  value.name = "Expression"
)

p <- ggplot(
  long_df,
  aes(x = Group,
      y = Expression,
      fill = Group)
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 1.5
  ) +
  facet_wrap(
    ~ miRNA,
    scales = "free_y",
    ncol = 3
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    title = "Six-miRNA Biomarker Panel (GSE259327)",
    x = "",
    y = "Expression"
  )

p

ggsave(
  "Six_miRNA_Boxplots_GSE259327.png",
  p,
  width = 12,
  height = 8,
  dpi = 600
)
head(rownames(expr259327))
grep("107", rownames(expr259327), value=TRUE)
grep("184", rownames(expr259327), value=TRUE)
grep("216", rownames(expr259327), value=TRUE)
grep("20a", rownames(expr259327), value=TRUE)
grep("134", rownames(expr259327), value=TRUE)
grep("215", rownames(expr259327), value=TRUE)


panel <- c(
  "miR-184",
  "miR-107",
  "miR-216b-5p",
  "miR-20a-5p",
  "miR-134-3p",
  "miR-215-3p"
)

panel %in% rownames(expr259327)

box_df <- as.data.frame(
  t(expr259327[panel, ])
)

box_df$Group <- group259327

dim(box_df)
head(box_df)

long_df$Expression <- log2(long_df$Expression + 1)

library(reshape2)
library(ggplot2)

long_df <- melt(
  box_df,
  id.vars = "Group",
  variable.name = "miRNA",
  value.name = "Expression"
)

p <- ggplot(
  long_df,
  aes(Group, Expression, fill = Group)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 1.2
  ) +
  facet_wrap(
    ~ miRNA,
    scales = "free_y",
    ncol = 3
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face="bold"),
    plot.title = element_text(
      hjust=0.5,
      face="bold"
    )
  ) +
  labs(
    title = "Six-miRNA Biomarker Panel (GSE259327)",
    x = "",
    y = "Expression"
  )

p

ggsave(
  "Six_miRNA_Boxplots_GSE259327.png",
  p,
  width = 12,
  height = 8,
  dpi = 600
)
find_gene_objects <- function(gene_name) {
  
  results <- character(0)
  
  # Only list objects from your saved global workspace
  object_names <- ls(envir = .GlobalEnv)
  
  for (nm in object_names) {
    
    x <- get(nm, envir = .GlobalEnv, inherits = FALSE)
    
    # Search only data frames or matrices
    if (is.data.frame(x) || is.matrix(x)) {
      
      values <- c(
        as.character(unlist(x, use.names = FALSE)),
        rownames(x),
        colnames(x)
      )
      
      values <- trimws(values[!is.na(values)])
      
      # Exact gene-name match
      if (any(toupper(values) == toupper(gene_name))) {
        results <- c(results, nm)
      }
    }
  }
  
  unique(results)
}



