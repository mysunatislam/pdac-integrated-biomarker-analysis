#!/usr/bin/env Rscript
# Publication rendering of saved results only. No statistical models are refitted.
options(stringsAsFactors = FALSE)
.libPaths(c(".R-library", .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gridExtra); library(pROC) })
out <- "reanalysis/figures/publication_revision"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create("reanalysis/provenance", recursive = TRUE, showWarnings = FALSE)
blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"; amber <- "#E69F00"
ink <- "#202428"; gray <- "#A3ABB2"; light <- "#ECEFF1"
font <- "Arial"; base <- 9
requested <- commandArgs(trailingOnly = TRUE)
inputs <- character(); export_checks <- list()
read <- function(...) {
  path <- file.path("reanalysis/results", ...)
  inputs[[path]] <<- unname(tools::md5sum(path))
  read.csv(path, check.names = FALSE)
}
pub <- function() theme_classic(base_size = base, base_family = font) +
  theme(axis.text = element_text(size = 8.5, colour = ink), axis.title = element_text(size = 9),
        axis.line = element_line(linewidth = .3, colour = ink), axis.ticks = element_line(linewidth = .25),
        strip.background = element_blank(), strip.text = element_text(size = 9, face = "bold", hjust = 0),
        legend.title = element_text(size = 8.5), legend.text = element_text(size = 8.5),
        legend.key.size = unit(3.5, "mm"), legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5),
        plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(4, 5, 4, 4),
        panel.spacing = unit(6, "mm"), plot.background = element_rect(fill = "white", colour = NA))
save <- function(plot, stem, width = 180, height = 110) {
  if (length(requested) && !stem %in% requested) return(invisible(NULL))
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height, units = "mm",
         dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width, height = height, units = "mm",
         device = cairo_pdf, family = font, bg = "white")
  ggsave(file.path(out, paste0(stem, ".svg")), plot, width = width, height = height, units = "mm",
         device = svglite::svglite, bg = "white")
  export_checks[[stem]] <<- data.frame(Figure = stem, Width_mm = width, Height_mm = height,
                                     PNG_dpi = 600, Base_font_pt = base, Formats = "PNG;PDF;SVG", Refitted = FALSE)
  cat("Exported", stem, "\n")
}
arrange <- function(...) arrangeGrob(...)
short <- c(Historical_six = "Historical six", Published_three_refit = "Published three (refit)",
           EV_ranked_nested = "EV-ranked selection", Unrestricted_nested = "Unrestricted selection")
h12 <- read("integrated", "historical12_current_evidence.csv")
perf <- read("plasma", "diagnostic_to_PLCO_performance.csv")
ds <- read("plasma", "diagnostic_model_scores.csv")
ps <- read("plasma", "PLCO_locked_model_scores.csv")

# A deliberately sparse workflow. Detailed qualifications belong in the legend.
workflow <- function() {
  boxes <- data.frame(x = rep(c(2.5, 7.5), 3), y = rep(c(8.7, 5.65, 2.6), each = 2),
    title = c("Tissue candidates", "Circulating miRNAs", "Corrected associations", "Prediction evaluation", "Robustness checks", "Clinical benchmarks"),
    text = c("Historical 12-gene subset\nSaved STRING topology", "Historical six / published three\nExploratory EV candidates",
             "Spatial: 8 inferred profiles\nPaired tissue: 35 donors\nADM culture: 14 donors", "Diagnostic: 121 PDAC / 82 controls\nLocked PLCO: 48 cases / 48 controls",
             "PDAC / CP: 195 / 59 specimens\nBroad FDR and omission checks", "Serum: 88 PC / 27 controls\nBinary CA19-9, age and sex"),
    fill = rep(c("#EAF3F8", "#EAF6F1"), 3))
  p <- ggplot() + geom_rect(data = boxes, aes(xmin = x - 2.25, xmax = x + 2.25, ymin = y - 1.05,
           ymax = y + 1.05, fill = fill), colour = gray, linewidth = .4) + scale_fill_identity() +
    geom_text(data = boxes, aes(x, y = y + .61, label = title), size = 3.35, fontface = "bold", family = font) +
    geom_text(data = boxes, aes(x, y = y - .15, label = text), size = 3.05, family = font, lineheight = 1.3) +
    annotate("segment", x = rep(c(2.5, 7.5), 2), xend = rep(c(2.5, 7.5), 2),
             y = rep(c(7.6, 4.55), each = 2), yend = rep(c(6.78, 3.73), each = 2),
             linewidth = .4, arrow = arrow(length = unit(1.8, "mm"))) +
    annotate("text", x = 5, y = .63, label = "Joint interpretation: association support and clinical limitations", size = 3.25,
             fontface = "bold", family = font) + coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), expand = FALSE) +
    theme_void(base_family = font) + theme(plot.margin = margin(4, 4, 4, 4))
  p
}
save(workflow(), "Figure_01", 180, 128)

# Shared grid network renderer. Edges are undirected, and node sizes are constant.
netdir <- "figure_revision/networks"
node_col <- c("Broad association support" = blue, "Early spatial association" = amber,
              "Other historical candidate" = "#D0D4D8", "Network context" = "#F3F4F5", "Cached miRNA query" = green)
network <- function(nodes, edges, width, height, title, positions, annotation = FALSE, font_size = 8.5) {
  coords <- merge(nodes, positions, by = "Gene", sort = FALSE)
  stopifnot(nrow(coords) == nrow(nodes), !anyDuplicated(coords$Gene))
  children <- list(textGrob(title, x = unit(3, "mm"), y = unit(height - 4, "mm"), just = "left",
                            gp = gpar(fontfamily = font, fontsize = 10, fontface = "bold", col = ink)))
  for (i in seq_len(nrow(edges))) {
    a <- coords[match(edges$source[i], coords$Gene), ]; b <- coords[match(edges$target[i], coords$Gene), ]
    children[[length(children) + 1]] <- segmentsGrob(unit(a$x, "mm"), unit(a$y, "mm"), unit(b$x, "mm"), unit(b$y, "mm"),
      gp = gpar(col = if (annotation) "#8399A5" else "#AAB3BA", lwd = .65, lty = if (annotation) "dashed" else "solid"))
  }
  for (i in seq_len(nrow(coords))) {
    n <- coords[i, ]; w <- if (n$Node_type == "miRNA") 29 else max(14, nchar(n$Gene) * 1.7 + 3)
    fill <- node_col[[n$Evidence_class]]
    children[[length(children) + 1]] <- roundrectGrob(unit(n$x, "mm"), unit(n$y, "mm"), unit(w, "mm"), unit(6.5, "mm"),
      r = unit(1, "mm"), gp = gpar(fill = fill, col = ink, lwd = .65))
    children[[length(children) + 1]] <- textGrob(sub("hsa-", "", n$Gene), x = unit(n$x, "mm"), y = unit(n$y, "mm"),
      gp = gpar(fontfamily = font, fontsize = font_size, fontface = if (n$Node_type == "Gene") "italic" else "plain",
                col = if (n$Evidence_class %in% c("Broad association support", "Cached miRNA query")) "white" else ink))
  }
  legend <- intersect(c("Broad association support", "Early spatial association", "Other historical candidate", "Network context", "Cached miRNA query"), unique(nodes$Evidence_class))
  lg_labels <- c("Broad association support" = "Three-setting support", "Early spatial association" = "Early spatial association",
                 "Other historical candidate" = "Other historical candidate", "Network context" = "Network context", "Cached miRNA query" = "miRNA query")
  xpos <- c(5, 65, 5, 65, 125); ypos <- c(12, 12, 4, 4, 12)
  for (i in seq_along(legend)) {
    children[[length(children) + 1]] <- rectGrob(unit(xpos[i], "mm"), unit(ypos[i], "mm"), unit(3, "mm"), unit(3, "mm"),
      gp = gpar(fill = node_col[[legend[i]]], col = ink, lwd = .5))
    children[[length(children) + 1]] <- textGrob(lg_labels[[legend[i]]], x = unit(xpos[i] + 3, "mm"), y = unit(ypos[i], "mm"), just = "left",
      gp = gpar(fontfamily = font, fontsize = 8.5, col = ink))
  }
  gTree(children = do.call(gList, children), vp = viewport(width = unit(width, "mm"), height = unit(height, "mm")))
}
nn <- read(netdir, "historical12_nodes.csv"); ee <- read(netdir, "historical12_edges.csv")
cycle <- c("TOP2A", "PCNA", "FOXM1", "CENPF", "MKI67", "CCNB1", "ANLN", "CCNA2")
theta <- seq(0, 2*pi, length.out = 9)[-9]
pos12 <- rbind(data.frame(Gene = cycle, x = 118 + 39*cos(theta), y = 62 + 29*sin(theta)),
               data.frame(Gene = c("CASP3", "ACTB", "ITGA2", "STAT1"), x = c(54, 27, 14, 45), y = c(58, 47, 27, 86)))
gnet <- network(nn, ee, 180, 106, "A  Saved STRING subnetwork: 12 candidates, 29 associations", pos12)
keys <- c("Spatial_ADM_Normal", "Spatial_PDAC_ADM", "Spatial_PDAC_Normal", "Paired_tumor_normal", "Culture_D6_D0", "Chip_PDAC_CP")
labels <- c("Spatial\nADM / Normal", "Spatial\nPDAC / ADM", "Spatial\nPDAC / Normal", "Paired tissue\nTumour / Normal", "ADM culture\nDay 6 / Day 0", "Chip adjusted\nPDAC / CP")
long <- do.call(rbind, lapply(seq_along(keys), function(i) data.frame(Gene = h12$Gene, Contrast = labels[i],
               Effect = h12[[paste0(keys[i], "_effect")]], BH = h12[[paste0(keys[i], "_BH")]])))
long$Gene <- factor(long$Gene, levels = rev(h12$Gene)); long$Contrast <- factor(long$Contrast, levels = labels)
long$Label <- paste0(sprintf("%+.2f", long$Effect), ifelse(long$BH < .05, "*", ""))
p <- ggplot(long, aes(Contrast, Gene)) + geom_tile(fill = "#F1F3F5", colour = "white", linewidth = .5) +
  geom_tile(data = long[long$BH < .05, ], aes(fill = Effect), colour = "white", linewidth = .5) +
  geom_text(aes(label = Label), size = 2.95, family = font) +
  scale_fill_gradient2(low = "#A6CDE2", mid = "#FFFFFF", high = "#EDAA87", limits = c(-3, 3), name = "Log2 effect\n(BH < 0.05)") +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL, title = "B  Corrected expression associations") + pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(face = "italic"),
        legend.position = "right") + guides(fill = guide_colourbar(display = "rectangles", barheight = unit(50, "mm"), barwidth = unit(3.5, "mm")))
save(arrange(gnet, ggplotGrob(p), ncol = 1, heights = unit(c(106, 97), "mm")), "Figure_02", 180, 203)

sp <- read("spatial", "GSE208536_patient_blocked_all_targets.csv")
joint <- read("integrated_tissue", "jointly_measurable_spatial_associated_genes.csv")
top <- head(joint$Gene[order(joint$Patient_Blocked_F_BH_FDR_all_targets)], 20)
mp <- read("input_audit", "GSE208536_roi_mapping.csv")
mp$Profile <- sprintf("P%02d", match(mp$Clinical_Profile, unique(mp$Clinical_Profile)))
path <- "data/source/GSE208536_Processed_data.xlsx"; inputs[[path]] <- unname(tools::md5sum(path))
xl <- readxl::read_excel(path, sheet = "Sheet1"); raw <- as.matrix(xl[, -1]); storage.mode(raw) <- "double"
rownames(raw) <- as.character(xl[[1]]); stopifnot(identical(colnames(raw), mp$Expression_Column))
groups <- split(seq_len(nrow(mp)), paste(mp$Profile, mp$Condition))
means <- vapply(groups, function(ii) rowMeans(log2(raw[top, ii, drop = FALSE] + 1)), numeric(20))
cm <- do.call(rbind, strsplit(colnames(means), " "))
meta <- data.frame(Column = colnames(means), Profile = cm[, 1], State = factor(cm[, 2], levels = c("Normal", "ADM", "PDAC")))
meta <- meta[order(meta$State, meta$Profile), ]; z <- t(scale(t(means[, meta$Column])))
heat <- data.frame(Gene = rep(top, ncol(z)), Column = rep(colnames(z), each = 20), Z = as.vector(z))
heat$Gene <- factor(heat$Gene, levels = rev(top)); heat$Column <- factor(heat$Column, levels = meta$Column)
heat$State <- meta$State[match(heat$Column, meta$Column)]
p <- ggplot(heat, aes(Column, Gene, fill = Z)) + geom_tile(colour = "white", linewidth = .12) +
  facet_grid(. ~ State, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = blue, mid = "white", high = orange, limits = c(-2.5, 2.5), oob = scales::squish, name = "Row z-score") +
  scale_x_discrete(labels = function(x) sub(" .*", "", x)) + labs(x = "Inferred clinical profile (profile-state means)", y = NULL) + pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), panel.spacing.x = unit(3, "mm")) +
  guides(fill = guide_colourbar(display = "rectangles", barwidth = unit(30, "mm"), barheight = unit(3, "mm")))
save(p, "Figure_03", 180, 125)

cp <- read("publication", "E-MTAB-1791_frozen_spatial_family.csv")
cp <- cp[cp$Variant == "Unadjusted_for_composition", ]
cp$Model <- factor(cp$Model, levels = c("Unadjusted", "Chip_adjusted"), labels = c("A  No chip covariate", "B  Chip adjusted"))
cp$Evidence <- ifelse(cp$BH_frozen_family < .05 & cp$Same_direction, "BH < 0.05; concordant", "Other measured genes")
p <- ggplot(cp, aes(Spatial_log2FC, log2FC)) + geom_hline(yintercept = 0, colour = light, linewidth = .3) +
  geom_vline(xintercept = 0, colour = light, linewidth = .3) + geom_point(aes(colour = Evidence), size = .9, alpha = .65) +
  geom_point(data = cp[cp$Gene == "STAT1", ], shape = 21, fill = "white", colour = ink, size = 2.3, stroke = .5) +
  facet_wrap(~Model, nrow = 1) + scale_colour_manual(values = c("BH < 0.05; concordant" = blue, "Other measured genes" = gray), name = NULL) +
  labs(x = "Spatial PDAC / Normal (log2 effect)", y = "External PDAC / CP (log2 effect)") + pub()
save(p, "Figure_04", 180, 105)

sm <- read("publication", "STAT1_composition_sensitivity.csv")
vlabels <- c(Unadjusted_for_composition = "Base model", Plus_T_cells = "+ T-cell score", Plus_Fibroblasts = "+ Fibroblast score", Plus_both = "+ Both scores")
sm$Variant_label <- factor(vlabels[sm$Variant], levels = rev(vlabels))
sl <- c("GSE208536 ADM_vs_Normal" = "A  Spatial: ADM / Normal", "GSE208536 PDAC_vs_ADM" = "B  Spatial: PDAC / ADM",
        "GSE208536 PDAC_vs_Normal" = "C  Spatial: PDAC / Normal", "GSE15471 PDAC_vs_Normal" = "D  Paired tissue: tumour / normal",
        "GSE143754 PDAC_vs_CP" = "E  Bulk tissue: PDAC / CP", "GSE179248 Culture_D6_vs_D0" = "F  ADM culture: day 6 / day 0")
sm$Cohort <- factor(paste(sm$Dataset, sm$Comparison), levels = names(sl), labels = sl)
estimated <- sm[sm$Status == "Estimated", ]; missing <- sm[sm$Status != "Estimated", ]
p <- ggplot(estimated, aes(log2FC, Variant_label)) + geom_vline(xintercept = 0, colour = gray, linetype = 2, linewidth = .3) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .12, colour = ink, linewidth = .4) +
  geom_point(aes(colour = Variant), size = 2) +
  geom_text(data = missing, aes(x = .15, label = "Not estimable"), size = 2.95, family = font, hjust = 0, colour = "#666666") +
  facet_wrap(~Cohort, ncol = 2, drop = FALSE) + scale_y_discrete(drop = FALSE) +
  scale_colour_manual(values = c(Unadjusted_for_composition = ink, Plus_T_cells = blue, Plus_Fibroblasts = orange, Plus_both = green), guide = "none") +
  labs(x = "STAT1 log2 effect (pointwise 95% CI)", y = NULL) + pub()
save(p, "Figure_05", 180, 154)

curve <- function(y, score, model) {
  r <- pROC::roc(y, score, direction = "<", levels = c(0, 1), quiet = TRUE)
  data.frame(Model = model, FPR = 1 - r$specificities, TPR = r$sensitivities)
}
curves <- do.call(rbind, lapply(perf$Model, function(m) curve(ps$Outcome, ps[[m]], short[[m]])))
pl <- perf; pl$Model_label <- factor(short[pl$Model], levels = short[perf$Model])
pl$Text <- sprintf("AUC %.3f\n95%% CI %.3f-%.3f", pl$PLCO_AUC, pl$PLCO_AUC_CI_low, pl$PLCO_AUC_CI_high)
curves$Model_label <- factor(curves$Model, levels = levels(pl$Model_label))
p <- ggplot(curves, aes(FPR, TPR)) + geom_abline(slope = 1, intercept = 0, colour = gray, linetype = 2, linewidth = .3) +
  geom_path(colour = blue, linewidth = .6) + geom_text(data = pl, aes(x = .98, y = .06, label = Text), inherit.aes = FALSE,
    size = 3, family = font, hjust = 1, vjust = 0, lineheight = 1.2) + facet_wrap(~Model_label, ncol = 2) +
  scale_x_continuous(breaks = c(0, .5, 1)) + scale_y_continuous(breaks = c(0, .5, 1)) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "False-positive rate", y = "True-positive rate") + pub()
save(p, "Figure_06", 180, 166)

ss <- read("external", "GSE85589_CA19_9_out_of_fold_scores.csv")
sv <- read("external", "GSE85589_CA19_9_nested_CV_performance.csv")
sn <- c(CA19_9_age_sex = "Binary CA19-9, age and sex", CA19_9_plus_published3 = "Baseline + published three", CA19_9_plus_historical6 = "Baseline + historical six")
sc <- do.call(rbind, lapply(sv$Model, function(m) curve(ss$Case, ss[[m]], m)))
sc$Model <- factor(sc$Model, levels = sv$Model)
p <- ggplot(sc, aes(FPR, TPR, colour = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, colour = gray, linewidth = .3, linetype = 2) + geom_path(linewidth = .65) +
  scale_colour_manual(values = c(ink, blue, orange), labels = paste0(sn[sv$Model], " (AUC ", sprintf("%.3f", sv$AUC), ")"), name = NULL) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash"), labels = paste0(sn[sv$Model], " (AUC ", sprintf("%.3f", sv$AUC), ")"), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + labs(x = "False-positive rate", y = "True-positive rate") + pub() +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
save(p, "Figure_07", 135, 142)

volcano <- function(d, effect, pcol, qcol, tested, name) {
  stopifnot(nrow(d) == tested, all(is.finite(d[[effect]])), all(is.finite(d[[pcol]])), all(is.finite(d[[qcol]])))
  stopifnot(max(abs(p.adjust(d[[pcol]], "BH") - d[[qcol]])) < 1e-12)
  cl <- ifelse(d[[qcol]] < .05, "BH < 0.05", ifelse(d[[pcol]] < .05, "Nominal only", "P >= 0.05"))
  dd <- rbind(data.frame(Effect = d[[effect]], Height = -log10(pmax(d[[pcol]], 1e-300)), Class = cl, Panel = "A  Nominal P"),
              data.frame(Effect = d[[effect]], Height = -log10(pmax(d[[qcol]], 1e-300)), Class = cl, Panel = "B  BH-adjusted P"))
  a <- aggregate(Height ~ Panel, dd, max); a$Height <- pmax(a$Height, -log10(.05)) * 1.24
  a$Effect <- min(dd$Effect); a$Text <- c(sprintf("%s features\n%s nominal P < 0.05", format(tested, big.mark = ","), sum(d[[pcol]] < .05)),
                                       sprintf("%s / %s BH < 0.05", sum(d[[qcol]] < .05), format(tested, big.mark = ",")))
  ggplot(dd, aes(Effect, Height)) + geom_blank(data = a) + geom_hline(yintercept = -log10(.05), linetype = 2, colour = ink, linewidth = .3) +
    geom_vline(xintercept = 0, colour = light, linewidth = .3) + geom_point(aes(colour = Class), size = .65, alpha = .6) +
    geom_text(data = a, aes(label = Text), hjust = 0, vjust = 1, family = font, size = 3, lineheight = 1.2) +
    facet_wrap(~Panel, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = c("P >= 0.05" = gray, "Nominal only" = orange, "BH < 0.05" = blue), name = NULL) +
    labs(x = if (name == "bulk") "PDAC / CP (log2 fold change)" else "PDAC - pooled CP/IPMN (median log2 difference)",
         y = "-log10(P or BH-adjusted P)") + pub() + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
}
bk <- read("bulk", "GSE143754_gene_level_contrasts.csv")
ev <- read("ev", "GSE304572_EV_differential.csv")
save(volcano(bk, "PDAC_vs_CP_log2FC", "PDAC_vs_CP_p", "PDAC_vs_CP_BH_FDR", 23207, "bulk"), "Figure_S01", 180, 104)
save(volcano(ev, "Median_Difference_log2CPM", "Wilcoxon_p", "BH_FDR", 440, "ev"), "Figure_S02", 180, 104)

age <- read("bulk", "GSE143754_sample_metadata.csv")
age$Age <- as.numeric(sub("^age: ", "", age$Sample_characteristics_ch1)); age$Sex <- sub("^Sex: ", "", age$Sample_characteristics_ch1.1)
counts <- table(age$Group)
p <- ggplot(age, aes(Group, Age)) + geom_boxplot(width = .4, outlier.shape = NA, colour = gray, fill = light, linewidth = .4) +
  geom_point(aes(colour = Sex, shape = Sex), size = 2, position = position_jitter(width = .11, height = 0, seed = 203)) +
  scale_colour_manual(values = c(Female = orange, Male = blue), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", counts[x], ")")) + labs(x = NULL, y = "Age (years)") + pub()
save(p, "Figure_S03", 135, 95)
qc <- read("plasma", "sample_assay_qc.csv")
qc$Group <- factor(qc$Outcome, levels = c(0, 1), labels = c("Control", "PDAC"))
p <- ggplot(qc, aes(Cohort, log10(Total_Counts), fill = Group)) + geom_boxplot(outlier.shape = NA, width = .6, linewidth = .4) +
  geom_point(aes(group = Group), position = position_jitterdodge(jitter.width = .18, dodge.width = .75, seed = 42),
             size = .55, alpha = .45) + scale_fill_manual(values = c(Control = "#B4D9ED", PDAC = "#EFAF91"), name = NULL) +
  labs(x = NULL, y = "Log10 deposited total counts") + pub()
save(p, "Figure_S04", 135, 97)

cal <- read("qc", "PLCO_calibration_bins.csv"); cal$Model <- factor(cal$Model, levels = gsub("_", " ", names(short)), labels = short)
p <- ggplot(cal, aes(Mean_Score, Observed_Fraction)) + geom_abline(slope = 1, intercept = 0, colour = gray, linetype = 2, linewidth = .3) +
  geom_errorbar(aes(ymin = Observed_CI_low, ymax = Observed_CI_high), width = .012, colour = blue, linewidth = .4) +
  geom_point(colour = orange, size = 1.7) + facet_wrap(~Model, ncol = 2) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean score in rank bin", y = "Observed case fraction (95% CI)") + pub()
save(p, "Figure_S05", 180, 140)
cc <- read("secondary_ev", "cohort_composition.csv"); cc$Group <- factor(cc$Group, levels = cc$Group)
p <- ggplot(cc, aes(Group, Samples)) + geom_col(fill = blue, width = .65) +
  geom_text(aes(label = Samples), vjust = -.4, size = 3.1, family = font) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) + labs(x = NULL, y = "Participants") + pub()
save(p, "Figure_S06", 135, 84)

fam <- read("external", "tissue_frozen_spatial_family.csv")
fam <- fam[fam$Comparison %in% c("PDAC_vs_paired_Normal", "PDAC_vs_normal", "Culture_D6_vs_D0"), ]
fl <- c("GSE15471 PDAC_vs_paired_Normal" = "A  Paired tissue\n35 donor pairs", "GSE91035 PDAC_vs_normal" = "B  Commercial normal\n25 PDAC / 8 normal", "GSE179248 Culture_D6_vs_D0" = "C  ADM culture\n14 donor pairs")
fam$Cohort <- factor(paste(fam$Dataset, fam$Comparison), levels = names(fl), labels = fl)
fam$Direction <- ifelse(fam$Same_direction, "Same direction", "Opposite direction")
p <- ggplot(fam, aes(Spatial_log2FC, log2FC)) + geom_hline(yintercept = 0, colour = light, linewidth = .3) +
  geom_vline(xintercept = 0, colour = light, linewidth = .3) + geom_point(aes(colour = Direction), size = .8, alpha = .6) +
  facet_wrap(~Cohort, nrow = 1, scales = "free") + scale_colour_manual(values = c("Same direction" = blue, "Opposite direction" = orange), name = NULL) +
  labs(x = "Spatial reference contrast (log2 effect)", y = "External contrast (log2 effect)") + pub()
save(p, "Figure_S07", 180, 95)

ext <- read("external", "tissue_all_gene_contrasts.csv")
ext <- ext[ext$Gene == "STAT1" & ext$Dataset %in% c("GSE15471", "GSE179248"), ]
s <- sp[sp$Gene == "STAT1", ]; b <- bk[bk$Gene == "STAT1", ]
forest <- rbind(data.frame(Label = c("Spatial: ADM / Normal", "Spatial: PDAC / ADM"), Effect = c(s$Normal_to_ADM_log2, s$ADM_to_PDAC_log2),
                          Low = c(s$Normal_to_ADM_CI_low, s$ADM_to_PDAC_CI_low), High = c(s$Normal_to_ADM_CI_high, s$ADM_to_PDAC_CI_high)),
                data.frame(Label = ifelse(ext$Dataset == "GSE15471", "Paired tissue: tumour / normal", "ADM culture: day 6 / day 0"),
                           Effect = ext$log2FC, Low = ext$CI_low, High = ext$CI_high),
                data.frame(Label = "Bulk tissue: PDAC / CP", Effect = b$PDAC_vs_CP_log2FC, Low = b$PDAC_vs_CP_CI_low, High = b$PDAC_vs_CP_CI_high))
forest$Label <- factor(forest$Label, levels = rev(forest$Label))
p <- ggplot(forest, aes(Effect, Label)) + geom_vline(xintercept = 0, colour = gray, linetype = 2, linewidth = .3) +
  geom_errorbar(aes(xmin = Low, xmax = High), orientation = "y", width = .15, linewidth = .4) + geom_point(size = 2, colour = blue) +
  labs(x = "STAT1 log2 effect (pointwise 95% CI)", y = NULL) + pub()
save(p, "Figure_S08", 180, 83)
stat <- read("spatial", "STAT1_patient_stage_means.csv")
stat$Condition <- factor(stat$Condition, levels = c("Normal", "ADM", "PDAC"))
p <- ggplot(stat, aes(Condition, Expression, group = Patient)) + geom_line(colour = gray, linewidth = .4) +
  geom_point(aes(colour = Condition), size = 2) + scale_colour_manual(values = c(Normal = blue, ADM = amber, PDAC = orange), guide = "none") +
  labs(x = NULL, y = "Mean log2(normalized expression + 1)") + pub()
save(p, "Figure_S09", 135, 94)

score <- do.call(rbind, lapply(names(short), function(m) rbind(
  data.frame(Model = short[[m]], Cohort = "Diagnostic", Outcome = ds$Outcome, Score = ds[[paste0(m, "_Fitted")]]),
  data.frame(Model = short[[m]], Cohort = "PLCO", Outcome = ps$Outcome, Score = ps[[m]]))))
score$Outcome <- factor(score$Outcome, levels = c(0, 1), labels = c("Control", "PDAC"))
threshold <- data.frame(Model = short[perf$Model], Threshold = perf$Training_Youden_Threshold)
p <- ggplot(score, aes(Cohort, log10(pmax(Score, 1e-10)), fill = Outcome)) +
  geom_boxplot(outlier.size = .6, linewidth = .35, width = .6) +
  geom_hline(data = threshold, aes(yintercept = log10(Threshold)), inherit.aes = FALSE, linetype = 2, colour = ink, linewidth = .35) +
  facet_wrap(~Model, ncol = 2) + scale_fill_manual(values = c(Control = "#B4D9ED", PDAC = "#EFAF91"), name = NULL) +
  labs(x = NULL, y = "Log10 model probability") + pub()
save(p, "Figure_S10", 180, 125)
aucs <- read("external", "serum_fixed_direction_marker_AUC.csv")
ref <- read("external", "frozen_miRNA_reference.csv"); nine <- ref$miRNA[ref$Historical6 | ref$Published3]
aucs <- aucs[aucs$miRNA %in% nine & aucs$Subgroup == "All", ]; aucs$miRNA <- factor(aucs$miRNA, levels = rev(nine))
aucs$Comparison <- factor(aucs$Comparison, levels = c("PC_vs_Healthy", "PC_vs_Benign", "PC_vs_Other_cancer"), labels = c("PC / healthy", "PC / benign disease", "PC / other cancers"))
p <- ggplot(aucs, aes(AUC, miRNA)) + geom_vline(xintercept = .5, colour = gray, linetype = 2, linewidth = .3) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .12, linewidth = .35, colour = gray) +
  geom_point(size = 1.5, colour = blue) + facet_grid(Dataset ~ Comparison) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) + scale_y_discrete(drop = FALSE) +
  labs(x = "Fixed-direction marker AUC (95% CI)", y = NULL) + pub()
save(p, "Figure_S11", 180, 145)

ta <- read("integrated", "historical_miRNA_target_evidence.csv")
ta$Gene <- factor(ta$Gene, levels = rev(h12$Gene)); ta$MicroRNA <- sub("hsa-", "", ta$MicroRNA)
ta$Support <- ifelse(ta$Validated_Record_Count > 0, "Cached record present", "No cached record")
p <- ggplot(ta, aes(MicroRNA, Gene, fill = Support)) + geom_tile(colour = "white", linewidth = .5) +
  scale_fill_manual(values = c("Cached record present" = "#8DCCB8", "No cached record" = light), name = NULL) +
  labs(x = NULL, y = NULL) + pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(face = "italic"))
save(p, "Figure_S12", 100, 124)

# Full recovered network, not a freshly queried or recomputed hub-ranking graph.
full_n <- read(netdir, "ppi_nodes.csv"); full_e <- read(netdir, "ppi_edges.csv")
g <- igraph::graph_from_data_frame(full_e[, c("source", "target")], directed = FALSE, vertices = full_n$Gene)
set.seed(208536); xy <- igraph::layout_with_kk(g)
normalize <- function(v, lo, hi) lo + (v - min(v)) / diff(range(v)) * (hi - lo)
positions <- data.frame(Gene = igraph::V(g)$name, x = normalize(xy[, 1], 17, 163), y = normalize(xy[, 2], 30, 145))
widths <- pmax(14, nchar(positions$Gene) * 1.7 + 3)
# Resolve rectangular label overlaps while retaining the library-derived layout.
for (iteration in seq_len(1000)) {
  collisions <- 0
  for (i in seq_len(nrow(positions) - 1)) for (j in (i + 1):nrow(positions)) {
    dx <- positions$x[j] - positions$x[i]; dy <- positions$y[j] - positions$y[i]
    ox <- (widths[i] + widths[j])/2 + 1.2 - abs(dx); oy <- 8.0 - abs(dy)
    if (ox > 0 && oy > 0) {
      collisions <- collisions + 1
      if (ox < oy) { delta <- (ox/2 + .01)*ifelse(dx >= 0, 1, -1); positions$x[i] <- positions$x[i] - delta; positions$x[j] <- positions$x[j] + delta
      } else { delta <- (oy/2 + .01)*ifelse(dy >= 0, 1, -1); positions$y[i] <- positions$y[i] - delta; positions$y[j] <- positions$y[j] + delta }
    }
  }
  positions$x <- pmax(widths/2 + 3, pmin(180 - widths/2 - 3, positions$x)); positions$y <- pmax(27, pmin(145, positions$y))
  if (collisions == 0) break
}
if (collisions > 0) {
  # Snap to well-spaced sites when continuous label repulsion cannot resolve a dense component.
  sites <- expand.grid(x = seq(16.5, 163.5, length.out = 8), y = seq(28, 143, length.out = 9))
  for (i in order(igraph::degree(g), decreasing = TRUE)) {
    distance <- (sites$x - positions$x[i])^2 + (sites$y - positions$y[i])^2
    k <- which.min(distance); positions$x[i] <- sites$x[k]; positions$y[i] <- sites$y[k]
    sites <- sites[-k, , drop = FALSE]
  }
}
for (i in seq_len(nrow(positions) - 1)) for (j in (i + 1):nrow(positions)) {
  stopifnot(abs(positions$x[i] - positions$x[j]) >= (widths[i] + widths[j])/2 + 1 || abs(positions$y[i] - positions$y[j]) >= 8)
}
save(network(full_n, full_e, 180, 159, "Saved STRING network: 56 nodes, 112 associations", positions), "Figure_S13", 180, 159)
write.csv(positions, file.path("reanalysis/results", netdir, "ppi_publication_positions.csv"), row.names = FALSE)
write.csv(pos12, file.path("reanalysis/results", netdir, "historical12_publication_positions.csv"), row.names = FALSE)

mir_n <- read(netdir, "mirna_nodes.csv"); mir_e <- read(netdir, "mirna_edges.csv")
# Bipartite placement keeps every annotated target visible and every edge traceable.
mirpos <- rbind(data.frame(Gene = h12$Gene, x = 127, y = seq(145, 30, length.out = 12)),
                data.frame(Gene = c("hsa-miR-107", "hsa-miR-20a-5p"), x = 38, y = c(118, 58)))
save(network(mir_n, mir_e, 180, 164, "Cached target annotations: 2 miRNAs, 12 genes, 22 links", mirpos, annotation = TRUE), "Figure_S14", 180, 164)
write.csv(mirpos, file.path("reanalysis/results", netdir, "mirna_publication_positions.csv"), row.names = FALSE)

go <- read("enrichment", "GO_BP_jointly_measurable_background.csv")
go <- head(go[order(go$BH_FDR, go$GO_ID), ], 12)
go$Label <- vapply(paste0(go$Term, " (", go$GO_ID, ")"), function(s) paste(strwrap(s, width = 50), collapse = "\n"), character(1))
go$Label <- factor(go$Label, levels = rev(go$Label))
p <- ggplot(go, aes(Enrichment_Ratio, Label)) + geom_segment(aes(x = 1, xend = Enrichment_Ratio, yend = Label), colour = light, linewidth = .8) +
  geom_point(aes(size = Selected_Hits, colour = -log10(BH_FDR)), alpha = .95) +
  scale_colour_gradient(low = "#58A7CD", high = "#00557F", name = "-log10(BH FDR)") +
  scale_size_continuous(range = c(2, 4.5), name = "Selected genes") + labs(x = "Enrichment ratio", y = NULL) + pub() +
  guides(colour = guide_colourbar(display = "rectangles", barwidth = unit(30, "mm"), barheight = unit(3, "mm")))
save(p, "Figure_S15", 180, 145)

stopifnot(all(vapply(names(inputs), function(p) identical(unname(tools::md5sum(p)), inputs[[p]]), logical(1))))
design_path <- "reanalysis/provenance/publication_figure_design_checks.csv"
design <- do.call(rbind, export_checks)
if (length(requested) && file.exists(design_path)) {
  previous <- read.csv(design_path)
  design <- rbind(previous[!previous$Figure %in% design$Figure, ], design)
  design <- design[order(design$Figure), ]
}
write.csv(design, design_path, row.names = FALSE)
jsonlite::write_json(as.list(inputs), "reanalysis/provenance/figure_revision_input_md5.json", pretty = TRUE, auto_unbox = TRUE)
cat(length(export_checks), "figure sets exported; saved analysis inputs unchanged.\n")
