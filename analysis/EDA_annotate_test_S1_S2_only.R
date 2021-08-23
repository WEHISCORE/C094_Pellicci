# EDA: annotate (S1 and S2 only)
# William Ho
# 2021-08-12

library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(cowplot)
library(edgeR)
library(Glimma)
library(BiocParallel)
library(patchwork)
library(janitor)
library(pheatmap)
library(batchelor)
library(rmarkdown)
library(BiocStyle)
library(readxl)
library(dplyr)
library(tidyr)
library(ggrepel)
library(magrittr)

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.single-cell.merged.S1_S2_only.SCE.rds"))

# pre-create directories for saving export, or error (dir not exists)
dir.create(here("data", "marker_genes", "S1_S2_only"), recursive = TRUE)
dir.create(here("output", "marker_genes", "S1_S2_only"), recursive = TRUE)

# Some useful colours
plate_number_colours <- setNames(
  unique(sce$colours$plate_number_colours),
  unique(names(sce$colours$plate_number_colours)))
plate_number_colours <- plate_number_colours[levels(sce$plate_number)]
tissue_colours <- setNames(
  unique(sce$colours$tissue_colours),
  unique(names(sce$colours$tissue_colours)))
tissue_colours <- tissue_colours[levels(sce$tissue)]
donor_colours <- setNames(
  unique(sce$colours$donor_colours),
  unique(names(sce$colours$donor_colours)))
donor_colours <- donor_colours[levels(sce$donor)]
stage_colours <- setNames(
  unique(sce$colours$stage_colours),
  unique(names(sce$colours$stage_colours)))
stage_colours <- stage_colours[levels(sce$stage)]
group_colours <- setNames(
  unique(sce$colours$group_colours),
  unique(names(sce$colours$group_colours)))
group_colours <- group_colours[levels(sce$group)]
cluster_colours <- setNames(
  unique(sce$colours$cluster_colours),
  unique(names(sce$colours$cluster_colours)))
cluster_colours <- cluster_colours[levels(sce$cluster)]

# Some useful gene sets
mito_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME == "MT")]
ribo_set <- grep("^RP(S|L)", rownames(sce), value = TRUE)
# NOTE: A more curated approach for identifying ribosomal protein genes
#       (https://github.com/Bioconductor/OrchestratingSingleCellAnalysis-base/blob/ae201bf26e3e4fa82d9165d8abf4f4dc4b8e5a68/feature-selection.Rmd#L376-L380)
library(msigdbr)
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_set <- union(
  ribo_set,
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$gene_symbol)
ribo_set <- intersect(ribo_set, rownames(sce))
sex_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))]
pseudogene_set <- rownames(sce)[
  any(grepl("pseudogene", rowData(sce)$ENSEMBL.GENEBIOTYPE))]
# NOTE: get rid of psuedogene seems not to be good enough for HVG determination of this dataset
protein_coding_gene_set <- rownames(sce)[
  any(grepl("protein_coding", rowData(sce)$ENSEMBL.GENEBIOTYPE))]

# include FACS data
facs <- t(assays(altExp(sce, "FACS"))$pseudolog)
facs_markers <- grep("V525_50_A_CD4_BV510|B530_30_A_CD161_FITC", colnames(facs), value = TRUE)
facs_selected <- facs[,facs_markers]
head(facs_selected)
colnames(facs_selected) <- c("CD161", "CD4")
colData(sce) <- cbind(colData(sce), facs_selected)

# summary - UMAP
p1 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "cluster", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p2 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "stage", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = stage_colours, name = "stage")
p3 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p4 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p5 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p6 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "group", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = group_colours, name = "group")
(p1 | p2) / (p3 | p4) / (p5 | p6)

# summary - stacked barplot
p1 <- ggcells(sce) +
  geom_bar(aes(x = cluster, fill = cluster)) +
  coord_flip() +
  ylab("Number of samples") +
  theme_cowplot(font_size = 8) +
  scale_fill_manual(values = cluster_colours) +
  geom_text(stat='count', aes(x = cluster, label=..count..), hjust=1.5, size=2)
p2 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = stage),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = stage_colours) +
  theme_cowplot(font_size = 8)
p3 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = plate_number),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = plate_number_colours) +
  theme_cowplot(font_size = 8)
p4 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = tissue),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = tissue_colours) +
  theme_cowplot(font_size = 8)
p5 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = donor),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = donor_colours) +
  theme_cowplot(font_size = 8)
p6 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = group),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = group_colours) +
  theme_cowplot(font_size = 8)
(p1 | p2) / (p3 | p4) / (p5 | p6)

# block on plate
sce$block <- paste0(sce$plate_number)





# NOTE: To illustrate what would be the most relevant number of clusters for this subset, we try to tweak the k values

# backup
sce0 <- sce

# clustering (different k values)
# number of cluster = 6
set.seed(4759)
snn_gr_1 <- buildSNNGraph(sce0, use.dimred = "corrected", k=5)
clusters_1 <- igraph::cluster_louvain(snn_gr_1)
sce0$cluster_1 <- factor(clusters_1$membership)
cluster_colours_1 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_1))],
  levels(sce0$cluster_1))
sce0$colours$cluster_colours_1 <- cluster_colours_1[sce0$cluster_1]
# number of cluster = 5
set.seed(4759)
snn_gr_2 <- buildSNNGraph(sce0, use.dimred = "corrected", k=8)
clusters_2 <- igraph::cluster_louvain(snn_gr_2)
sce0$cluster_2 <- factor(clusters_2$membership)
cluster_colours_2 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_2))],
  levels(sce0$cluster_2))
sce0$colours$cluster_colours_2 <- cluster_colours_2[sce0$cluster_2]
# number of cluster = 4 (default)
set.seed(4759)
snn_gr_3 <- buildSNNGraph(sce0, use.dimred = "corrected", k=10)
clusters_3 <- igraph::cluster_louvain(snn_gr_3)
sce0$cluster_3 <- factor(clusters_3$membership)
cluster_colours_3 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_3))],
  levels(sce0$cluster_3))
sce0$colours$cluster_colours_3 <- cluster_colours_3[sce0$cluster_3]
# number of cluster = 3
set.seed(4759)
snn_gr_4 <- buildSNNGraph(sce0, use.dimred = "corrected", k=15)
clusters_4 <- igraph::cluster_louvain(snn_gr_4)
sce0$cluster_4 <- factor(clusters_4$membership)
cluster_colours_4 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_4))],
  levels(sce0$cluster_4))
sce0$colours$cluster_colours_4 <- cluster_colours_4[sce0$cluster_4]

# UMAP plot
p1 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p2 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p3 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p4 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")

p5 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "cluster_1", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_1, name = "cluster_1")
p6 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "cluster_2", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_2, name = "cluster_2")
p7 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "cluster_3", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_3, name = "cluster_3")
p8 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "cluster_4", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_4, name = "cluster_4")

p9 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "stage", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = stage_colours, name = "stage")
p10 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "stage", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = stage_colours, name = "stage")
p11 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "stage", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = stage_colours, name = "stage")
p12 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "stage", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = stage_colours, name = "stage")

p13 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p14 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p15 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p16 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")

p17 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p18 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p19 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p20 <- plotReducedDim(sce0, "UMAP_corrected", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")

p21 <- ggcells(sce0) + geom_bar(aes(x = cluster_1, fill = group), position = position_fill(reverse = TRUE)) +
  coord_flip() + ylab("Frequency") + scale_fill_manual(values = group_colours) + theme_cowplot(font_size = 8)
p22 <- ggcells(sce0) + geom_bar(aes(x = cluster_2, fill = group), position = position_fill(reverse = TRUE)) +
  coord_flip() + ylab("Frequency") + scale_fill_manual(values = group_colours) + theme_cowplot(font_size = 8)
p23 <- ggcells(sce0) + geom_bar(aes(x = cluster_3, fill = group), position = position_fill(reverse = TRUE)) +
  coord_flip() + ylab("Frequency") + scale_fill_manual(values = group_colours) + theme_cowplot(font_size = 8)
p24 <- ggcells(sce0) + geom_bar(aes(x = cluster_4, fill = group), position = position_fill(reverse = TRUE)) +
  coord_flip() + ylab("Frequency") + scale_fill_manual(values = group_colours) + theme_cowplot(font_size = 8)

p25 <- ggcells(sce0) + geom_bar(aes(x = cluster_1, fill = cluster_1)) + coord_flip() + ylab("Number of samples") +
  theme_cowplot(font_size = 8) + scale_fill_manual(values = cluster_colours_1) +
  geom_text(stat='count', aes(x = cluster_1, label=..count..), hjust=1.5, size=2)
p26 <- ggcells(sce0) + geom_bar(aes(x = cluster_2, fill = cluster_2)) + coord_flip() + ylab("Number of samples") +
  theme_cowplot(font_size = 8) + scale_fill_manual(values = cluster_colours_2) +
  geom_text(stat='count', aes(x = cluster_2, label=..count..), hjust=1.5, size=2)
p27 <- ggcells(sce0) + geom_bar(aes(x = cluster_3, fill = cluster_3)) + coord_flip() + ylab("Number of samples") +
  theme_cowplot(font_size = 8) + scale_fill_manual(values = cluster_colours_3) +
  geom_text(stat='count', aes(x = cluster_3, label=..count..), hjust=1.5, size=2)
p28 <- ggcells(sce0) + geom_bar(aes(x = cluster_4, fill = cluster_4)) + coord_flip() + ylab("Number of samples") +
  theme_cowplot(font_size = 8) + scale_fill_manual(values = cluster_colours_4) +
  geom_text(stat='count', aes(x = cluster_4, label=..count..), hjust=1.5, size=2)

p1 + p2 + p3 + p4 +
  p5 +p6 + p7 + p8 +
  p9 + p10 +p11 + p12 +
  p13 + p14 + p15 +p16 +
  p17 + p18 + p19 + p20 +
  p21 + p22 + p23 + p24 +
  p25 +p26 + p27 + p28 +
  plot_layout(ncol = 4, guides = "collect")

# comment:
# we expect there are at least 4 groups in this subset: thymus.s1, thymus.s2, blood.s1, blood.s2
# it seems like it can divide the Blood.s1/s2 from thymus.s1/s2
# 4 cluster is nice, but considering keeping the statistical power
3

###############################
# selected the best `number of cluster`
sce$cluster <- factor(clusters_4$membership)
# re-numbering of cluster according to collaborators' request
sce$cluster <- factor(
  dplyr::case_when(
    sce$cluster == "1" ~ "6",
    sce$cluster == "2" ~ "7",
    sce$cluster == "3" ~ "8"))

cluster_colours <- setNames(
  palette.colors(nlevels(sce$cluster), "Okabe-Ito"),
  levels(sce$cluster))










###################################
# (M1) raw unique
#
# cluster 6 (i.e. pure thymus.S1.S2) vs 7 (i.e. mostly thymus.S1.S2, more blood) vs 8 (i.e. mostly thymus.S1.S2, less blood)

# find unique DE ./. clusters
uniquely_up <- findMarkers(
  sce,
  groups = sce$cluster,
  block = sce$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  uniquely_up,
  here("data", "marker_genes", "S1_S2_only", "C094_Pellicci.uniquely_up.cluster_6_vs_7_vs_8.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S1_S2_only", "uniquely_up", "cluster_6_vs_7_vs_8"), recursive = TRUE)

vs_pair <- c("6", "7", "8")

message("Writing 'uniquely_up (cluster_6_vs_7_vs_8)' marker genes to file.")
for (n in names(uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S1_S2_only",
      "uniquely_up",
      "cluster_6_vs_7_vs_8",
      paste0("cluster_",
             vs_pair[which(names(uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(uniquely_up) %in% n)][1],
             "_vs_",
             vs_pair[-which(names(uniquely_up) %in% n)][2],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = uniquely_up[[n]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_ID"),
    file = gzout,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE)
  close(gzout)
}

##########################################
# look at cluster 6 (i.e. pure thymus.S1.S2)
chosen <- "6"
cluster6_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(pure thymus.S1.S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster6_uniquely_up <- cluster6_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster6_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster6_uniquely_up_noiseR <- cluster6_uniquely_up[setdiff(rownames(cluster6_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster6_uniquely_up_noiseR) %in% "CD4"),
       cluster6_uniquely_up_noiseR[which(rownames(cluster6_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster6_uniquely_up_noiseR) %in% "KLRB1"),
       cluster6_uniquely_up_noiseR[which(rownames(cluster6_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster6_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number,
    sce$CD4,
    sce$CD161),
  colour_columns_by = c(
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################
# look at cluster 7 (i.e. mostly thymus.S1.S2, more blood)
chosen <- "7"
cluster7_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly thymus.S1.S2, more blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster7_uniquely_up <- cluster7_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster7_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster7_uniquely_up_noiseR <- cluster7_uniquely_up[setdiff(rownames(cluster7_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster7_uniquely_up_noiseR) %in% "CD4"),
       cluster7_uniquely_up_noiseR[which(rownames(cluster7_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster7_uniquely_up_noiseR) %in% "KLRB1"),
       cluster7_uniquely_up_noiseR[which(rownames(cluster7_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster7_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number,
    sce$CD4,
    sce$CD161),
  colour_columns_by = c(
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################
# look at cluster 8 (i.e. mostly thymus.S1.S2, less blood)
chosen <- "8"
cluster8_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly thymus.S1.S2, less blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster8_uniquely_up <- cluster8_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster8_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster8_uniquely_up_noiseR <- cluster8_uniquely_up[setdiff(rownames(cluster8_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster8_uniquely_up_noiseR) %in% "CD4"),
       cluster8_uniquely_up_noiseR[which(rownames(cluster8_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster8_uniquely_up_noiseR) %in% "KLRB1"),
       cluster8_uniquely_up_noiseR[which(rownames(cluster8_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster8_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number,
    sce$CD4,
    sce$CD161),
  colour_columns_by = c(
    "cluster",
    "CD4",
    "CD161",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



# COMMENT:

# cluster 6 (i.e. pure thymus.s1.s2
# cluster 7 (i.e. mostly thymus.s1.s2, more blood
# cluster 8 (i.e. mostly thymus.s1.s2, less blood

# for global unique markers
# cluster 6: clear ACTG1, ACTB, GAPDH, MIR1244-3 and H3F3A expression as global unique markers
# cluster 7: clear HLA, CCL5, NKG7 expression as global unique markers (mostly independent of tissue)
# cluster 8: statistically no, but visually yes
# although they are mostly indistinguishable thymus.s1s2 cells with more/less/none blood.s1s2 cells
# this point out clearly at least cluster 6 and 7 are different clusters (and cluster 3 uniqueness may also show in pairwise unique)













###################################
# (M2) selected pairwise comparison

#########
# A vs B
#########

################################################################################
# cluster 6 (i.e. pure thymus.S1.S2) vs 7 (i.e. mostly thymus.S1.S2, more blood)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "6" | cp$cluster == "7"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs1 <- factor(ifelse(cp$cluster == 6, "A", "B"))

# set vs colours
vs1_colours <- setNames(
  palette.colors(nlevels(cp$vs1), "Set1"),
  levels(cp$vs1))
cp$colours$vs1_colours <- vs1_colours[cp$vs1]

# find unique DE ./. cluster-groups
vs1_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs1,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs1_uniquely_up,
  here("data", "marker_genes", "S1_S2_only", "C094_Pellicci.uniquely_up.cluster_6_vs_7.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S1_S2_only", "uniquely_up", "cluster_6_vs_7"), recursive = TRUE)

vs_pair <- c("6", "7")

message("Writing 'uniquely_up (cluster_6_vs_7)' marker genes to file.")
for (n in names(vs1_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S1_S2_only",
      "uniquely_up",
      "cluster_6_vs_7",
      paste0("cluster_",
             vs_pair[which(names(vs1_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs1_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs1_uniquely_up[[n]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_ID"),
    file = gzout,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE)
  close(gzout)
}

##########################################################
# look at cluster-group A / cluster 6 (i.e. pure thymus.S1.S2)
chosen <- "A"
A_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 6; pure thymus.S1.S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# A_uniquely_up_pcg <- A_uniquely_up[intersect(protein_coding_gene_set, rownames(A_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
A_uniquely_up_noiseR <- A_uniquely_up[setdiff(rownames(A_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(A_uniquely_up_noiseR) %in% "CD4"),
       A_uniquely_up_noiseR[which(rownames(A_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(A_uniquely_up_noiseR) %in% "KLRB1"),
       A_uniquely_up_noiseR[which(rownames(A_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- A_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs1,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs1",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs1 = vs1_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################################
# look at cluster-group B / cluster 7 (i.e. mostly thymus.S1.S2, more blood)
chosen <- "B"
B_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 7; mostly thymus.S1.S2, more blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# B_uniquely_up_pcg <- B_uniquely_up[intersect(protein_coding_gene_set, rownames(B_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
B_uniquely_up_noiseR <- B_uniquely_up[setdiff(rownames(B_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(B_uniquely_up_noiseR) %in% "CD4"),
       B_uniquely_up_noiseR[which(rownames(B_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(B_uniquely_up_noiseR) %in% "KLRB1"),
       B_uniquely_up_noiseR[which(rownames(B_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- B_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs1,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs1",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs1 = vs1_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



#########
# C vs D
#########

################################################################################
# cluster 6 (i.e. pure thymus.S1.S2) vs 8 (i.e. mostly thymus.S1.S2, less blood)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "6" | cp$cluster == "8"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs2 <- factor(ifelse(cp$cluster == 6, "C", "D"))

# set vs colours
vs2_colours <- setNames(
  palette.colors(nlevels(cp$vs2), "Set1"),
  levels(cp$vs2))
cp$colours$vs2_colours <- vs2_colours[cp$vs2]

# find unique DE ./. cluster-groups
vs2_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs2,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs2_uniquely_up,
  here("data", "marker_genes", "S1_S2_only", "C094_Pellicci.uniquely_up.cluster_6_vs_8.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S1_S2_only", "uniquely_up", "cluster_6_vs_8"), recursive = TRUE)

vs_pair <- c("1", "3")

message("Writing 'uniquely_up (cluster_6_vs_8)' marker genes to file.")
for (n in names(vs2_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S1_S2_only",
      "uniquely_up",
      "cluster_6_vs_8",
      paste0("cluster_",
             vs_pair[which(names(vs2_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs2_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs2_uniquely_up[[n]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_ID"),
    file = gzout,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE)
  close(gzout)
}

##########################################################
# look at cluster-group C / cluster 6 (i.e. pure thymus.S1.S2)
chosen <- "C"
C_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 6; pure thymus.S1.S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# C_uniquely_up_pcg <- C_uniquely_up[intersect(protein_coding_gene_set, rownames(C_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
C_uniquely_up_noiseR <- C_uniquely_up[setdiff(rownames(C_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(C_uniquely_up_noiseR) %in% "CD4"),
       C_uniquely_up_noiseR[which(rownames(C_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(C_uniquely_up_noiseR) %in% "KLRB1"),
       C_uniquely_up_noiseR[which(rownames(C_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- C_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs2,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs2",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs2 = vs2_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################################
# look at cluster-group D / cluster 8 (i.e. mostly thymus.S1.S2, less blood)
chosen <- "D"
D_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 8; mostly thymus.S1.S2, less blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# D_uniquely_up_pcg <- D_uniquely_up[intersect(protein_coding_gene_set, rownames(D_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
D_uniquely_up_noiseR <- D_uniquely_up[setdiff(rownames(D_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(D_uniquely_up_noiseR) %in% "CD4"),
       D_uniquely_up_noiseR[which(rownames(D_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(D_uniquely_up_noiseR) %in% "KLRB1"),
       D_uniquely_up_noiseR[which(rownames(D_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- D_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs2,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs2",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs2 = vs2_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)




#########
# E vs F
#########

##########################################################################
# cluster 6 (i.e. pure thymus.s1.s2) vs cluster 7_8 (i.e. mostly thymus.s1.s2 with blood.s1.s2)

# checkpoint
cp <- sce

# classify cluster-group for comparison
cp$vs3 <- factor(ifelse(cp$cluster == 6, "E", "F"))

# set vs colours
vs3_colours <- setNames(
  palette.colors(nlevels(cp$vs3), "Set1"),
  levels(cp$vs3))
cp$colours$vs3_colours <- vs3_colours[cp$vs3]

# find unique DE ./. cluster-groups
vs3_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs3,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs3_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_6_vs_7_8.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_6_vs_7_8"), recursive = TRUE)

vs_pair <- c("1", "2_3")

message("Writing 'uniquely_up (cluster_6_vs_7_8)' marker genes to file.")
for (n in names(vs3_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_6_vs_7_8",
      paste0("cluster_",
             vs_pair[which(names(vs3_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs3_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs3_uniquely_up[[n]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_ID"),
    file = gzout,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE)
  close(gzout)
}

###############################################################
# look at cluster-group E / cluster 6 (i.e. pure thymus.s1.s2)
chosen <- "E"
E_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 6; pure thymus.s1.s2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# E_uniquely_up_pcg <- E_uniquely_up[intersect(protein_coding_gene_set, rownames(E_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
E_uniquely_up_noiseR <- E_uniquely_up[setdiff(rownames(E_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(E_uniquely_up_noiseR) %in% "CD4"),
       E_uniquely_up_noiseR[which(rownames(E_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(E_uniquely_up_noiseR) %in% "KLRB1"),
       E_uniquely_up_noiseR[which(rownames(E_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- E_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs3,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs3",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs3 = vs3_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group F / cluster 7_8 (i.e. mostly thymus.s1.s2 with blood.s1.s2)
chosen <- "F"
F_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 7_8; mostly thymus.s1.s2 with blood.s1.s2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# F_uniquely_up_pcg <- F_uniquely_up[intersect(protein_coding_gene_set, rownames(F_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
F_uniquely_up_noiseR <- F_uniquely_up[setdiff(rownames(F_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(F_uniquely_up_noiseR) %in% "CD4"),
       F_uniquely_up_noiseR[which(rownames(F_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(F_uniquely_up_noiseR) %in% "KLRB1"),
       F_uniquely_up_noiseR[which(rownames(F_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- F_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs3,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs3",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs3 = vs3_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)






#########
# G vs H
#########

################################################################################
# cluster 7 (i.e. mostly thymus.S1.S2, more blood) vs 8 (i.e. mostly thymus.S1.S2, less blood)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "7" | cp$cluster == "8"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs4 <- factor(ifelse(cp$cluster == 7, "G", "H"))

# set vs colours
vs4_colours <- setNames(
  palette.colors(nlevels(cp$vs4), "Set1"),
  levels(cp$vs4))
cp$colours$vs4_colours <- vs4_colours[cp$vs4]

# find unique DE ./. cluster-groups
vs4_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs4,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs4_uniquely_up,
  here("data", "marker_genes", "S1_S2_only", "C094_Pellicci.uniquely_up.cluster_7_vs_8.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S1_S2_only", "uniquely_up", "cluster_7_vs_8"), recursive = TRUE)

vs_pair <- c("7", "8")

message("Writing 'uniquely_up (cluster_7_vs_8)' marker genes to file.")
for (n in names(vs4_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S1_S2_only",
      "uniquely_up",
      "cluster_2_vs_3",
      paste0("cluster_",
             vs_pair[which(names(vs4_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs4_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs4_uniquely_up[[n]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_ID"),
    file = gzout,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE)
  close(gzout)
}

##########################################################
# look at cluster-group G / cluster 7 (i.e. mostly thymus.S1.S2, more blood)
chosen <- "G"
G_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 7; mostly thymus.S1.S2, more blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# G_uniquely_up_pcg <- G_uniquely_up[intersect(protein_coding_gene_set, rownames(G_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
G_uniquely_up_noiseR <- G_uniquely_up[setdiff(rownames(G_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(G_uniquely_up_noiseR) %in% "CD4"),
       G_uniquely_up_noiseR[which(rownames(G_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(G_uniquely_up_noiseR) %in% "KLRB1"),
       G_uniquely_up_noiseR[which(rownames(G_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- G_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs4,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs4",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs4 = vs4_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################################
# look at cluster-group H / cluster 8 (i.e. mostly thymus.S1.S2, less blood)
chosen <- "H"
H_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 8; mostly thymus.S1.S2, less blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# H_uniquely_up_pcg <- H_uniquely_up[intersect(protein_coding_gene_set, rownames(H_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
H_uniquely_up_noiseR <- H_uniquely_up[setdiff(rownames(H_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(H_uniquely_up_noiseR) %in% "CD4"),
       H_uniquely_up_noiseR[which(rownames(H_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(H_uniquely_up_noiseR) %in% "KLRB1"),
       H_uniquely_up_noiseR[which(rownames(H_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- H_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs4,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number,
    cp$CD4,
    cp$CD161),
  colour_columns_by = c(
    "vs4",
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number",
    "CD4",
    "CD161"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      # TODO: temp trick to deal with the row-colouring problem
      # levels = c("Yes", "No")),
      levels = c("Yes")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs4 = vs4_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



# COMMENT:
#
# fx: how is thymus.s1s2.alike.by.more.blood is different with pure thymus.s1.s2 cells
#
# 1 vs 2 (A vs B)
# cluster 1 (pure thymus.s1.s2 >>> lots of marker up-regulated in cluster 1; in which, say STMN1 and MIR1244-4, seems to be more highly expression in S1 than S2
# cluster 2 (mostly thymus.s1.s2, more blood >>> number of unique marker found (say, IL7R, , LTB, CCL5, NKG7. etc)
# COMMENT: cluster 1 and 2 are 2 separated clusters
#
# 1 vs 3 (C vs D)
# cluster 1 (pure thymus.s1.s2 >>> 7 markers only, i.e. PTMA, MIR1244-2, and GAPDH, etc; maybe cluster 1 is more similar to 3, than to 2 ; related to amount of blood cells? 3 = less blood, 2 = more blood ???)
# cluster 3 (mostly thymus.s1.s2, less blood >>> with markers like LTB, and the more frequently expressed DGKA, EPB41
# COMMENT: except for LTB, there seems to have little number of genes in common between the pairwisely unique of 2 and 3, 2 and 3 supposed not to be the same
#
# 1 vs 2_3 (E vs F)
# cluster 1 (pure thymus.s1.s2 >>> got number of markers, eg ACTG1, GAPDH, PTMA, etc
# cluster 2_3 (mostly thymus.s1.s2 with more or less blood >>> also got lots of DE, eg.HLA, LTB, TXNIP; other markers like up-regulation of CCL5 and ILR7 clearly associated with cluster 2 only
# COMMENT: cluster 1 is clearly distinctive from cluster 2_3
#
# fx: how is thymus.s1s2.alike.by.more.blood and thymus.s1s2.alike.by.less.blood different from one and other
#
# 2 vs 3 (G vs H)
# cluster 2 (mostly thymus.s1.s2, more blood >>> got number of good markers, e.g. class 1 HLA, IFITM family, NKG7, which are all more frequently expressed in 2
# cluster 3 (mostly thymus.s1.s2, less blood >>> got another set of unique markers, eg. BCL11B, DGKA, FYB1, etc. for cluster 3
# COMMENT: there are clearly 2 distinct subtypes (though both are mostly thymus.s1.s2 with blood),





# CONCLUSION
# by default SNN setting, there are 4 clusters, but after testing, I found out 2 cluster show no sig difference, thus I pick 3 cluster for this subset
#
# S1 and S2 cells are basically inseparable as in terms of transcriptomes, they are highly overlapped
# having said so, we do find some clues about the key marker driven the S1 and S2 apart (e.g. STMN1 and MIR1244-4 seems to be more highly expression in S1 than S2; eg2 up-regulation of CCL5 and ILR7 clearly associated with cluster 2 only)
# also, we are able to spot the 3 unique transcriptional subtype of cells within the S1-S2 mix
# from cluster 1, we can identify markers true for pure thymus.s1.s2 (i.e. found in thymus only)
# for the thymus.s1.s2.alike.blood.bulk, it can be subdivided into 2 subtypes, of which their unique markers can be found in cluster2 and 3 (ie found in both thymus and blood)





























































###############################################
# Heatmap using minibulk sig markers as feature

# laod package for read in csv.gz
library(data.table)
library(R.utils)

# read in
a <- fread(here("output", "DEGs", "excluding_blood_1-3", "Thymus.S2_vs_Thymus.S1.aggregated_tech_reps.DEGs.csv.gz"))
b <- fread(here("output", "DEGs", "excluding_blood_1-3", "Thymus.S3_vs_Blood.S3.aggregated_tech_reps.DEGs.csv.gz"))
c <- fread(here("output", "DEGs", "excluding_blood_1-3", "Thymus.S3_vs_Thymus.S1.aggregated_tech_reps.DEGs.csv.gz"))
d <- fread(here("output", "DEGs", "excluding_blood_1-3", "Thymus.S3_vs_Thymus.S2.aggregated_tech_reps.DEGs.csv.gz"))

# extract DEGlist (FDR < 0.05)
minibulkDEG.a <- a$ENSEMBL.GENENAME[a$FDR<0.05]
minibulkDEG.b <- b$ENSEMBL.GENENAME[b$FDR<0.05]
minibulkDEG.c <- c$ENSEMBL.GENENAME[c$FDR<0.05]
minibulkDEG.d <- d$ENSEMBL.GENENAME[d$FDR<0.05]

# keep only unique markers
uniq.minibulkDEG.a <- Reduce(setdiff, list(minibulkDEG.a,
                                           minibulkDEG.b,
                                           minibulkDEG.c,
                                           minibulkDEG.d))
uniq.minibulkDEG.b <- Reduce(setdiff, list(minibulkDEG.b,
                                           minibulkDEG.a,
                                           minibulkDEG.c,
                                           minibulkDEG.d))
uniq.minibulkDEG.c <- Reduce(setdiff, list(minibulkDEG.c,
                                           minibulkDEG.a,
                                           minibulkDEG.b,
                                           minibulkDEG.d))
uniq.minibulkDEG.d <- Reduce(setdiff, list(minibulkDEG.d,
                                           minibulkDEG.a,
                                           minibulkDEG.b,
                                           minibulkDEG.c))

# check number of unique minibulkDEG in each
length(uniq.minibulkDEG.a)
length(uniq.minibulkDEG.b)
length(uniq.minibulkDEG.c)
length(uniq.minibulkDEG.d)

# keep only top50
top.uniq.minibulkDEG.a <- if(length(uniq.minibulkDEG.a) >=50){uniq.minibulkDEG.a[1:50]} else {uniq.minibulkDEG.a}
top.uniq.minibulkDEG.b <- if(length(uniq.minibulkDEG.b) >=50){uniq.minibulkDEG.b[1:50]} else {uniq.minibulkDEG.b}
top.uniq.minibulkDEG.c <- if(length(uniq.minibulkDEG.c) >=50){uniq.minibulkDEG.c[1:50]} else {uniq.minibulkDEG.c}
top.uniq.minibulkDEG.d <- if(length(uniq.minibulkDEG.d) >=50){uniq.minibulkDEG.d[1:50]} else {uniq.minibulkDEG.d}

# feature
minibulk_markers <- c(top.uniq.minibulkDEG.a,
                      top.uniq.minibulkDEG.b,
                      top.uniq.minibulkDEG.c,
                      top.uniq.minibulkDEG.d)

# plot heatmap
library(scater)
plotHeatmap(
  sce,
  features = minibulk_markers,
  columns = order(
    sce$cluster,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "stage",
    "tissue",
    "donor",
    "group",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  # TODO: temp trick to deal with the row-colouring problem
  annotation_row = data.frame(
    thymus.s2.vs.thymus.s1 = factor(ifelse(minibulk_markers %in% top.uniq.minibulkDEG.a, "DE", "not DE"), levels = c("DE")),
    thymus.s3.vs.blood.s3 = factor(ifelse(minibulk_markers %in% top.uniq.minibulkDEG.b, "DE", "not DE"), levels = c("DE")),
    thymus.s3.vs.thymus.s1 = factor(ifelse(minibulk_markers %in% top.uniq.minibulkDEG.c, "DE", "not DE"), levels = c("DE")),
    thymus.s3.vs.thymus.s2 = factor(ifelse(minibulk_markers %in% top.uniq.minibulkDEG.d, "DE", "not DE"), levels = c("DE")),
    row.names = minibulk_markers),
  main = "Row-normalized log expression of top unique markers from minibulk against different clusters",
  column_annotation_colors = list(
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)









