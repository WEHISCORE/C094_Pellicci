# EDA: annotate (thymus only)
# William Ho
# 2021-08-06

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

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.single-cell.merged.thymus_only.SCE.rds"))

# pre-create directories for saving export, or error (dir not exists)
dir.create(here("data", "marker_genes", "thymus_only"), recursive = TRUE)
dir.create(here("output", "marker_genes", "thymus_only"), recursive = TRUE)

# remove "Unknown" (as it is not informative at all)
sce <- sce[, sce$stage != "Unknown"]
colData(sce) <- droplevels(colData(sce))

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
# number of cluster = 2
set.seed(4759)
snn_gr_1 <- buildSNNGraph(sce0, use.dimred = "corrected", k=70)
clusters_1 <- igraph::cluster_louvain(snn_gr_1)
sce0$cluster_1 <- factor(clusters_1$membership)
cluster_colours_1 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_1))],
  levels(sce0$cluster_1))
sce0$colours$cluster_colours_1 <- cluster_colours_1[sce0$cluster_1]
# number of cluster = 3 (default)
set.seed(4759)
snn_gr_2 <- buildSNNGraph(sce0, use.dimred = "corrected", k=10)
clusters_2 <- igraph::cluster_louvain(snn_gr_2)
sce0$cluster_2 <- factor(clusters_2$membership)
cluster_colours_2 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_2))],
  levels(sce0$cluster_2))
sce0$colours$cluster_colours_2 <- cluster_colours_2[sce0$cluster_2]
# number of cluster = 4
set.seed(4759)
snn_gr_3 <- buildSNNGraph(sce0, use.dimred = "corrected", k=4)
clusters_3 <- igraph::cluster_louvain(snn_gr_3)
sce0$cluster_3 <- factor(clusters_3$membership)
cluster_colours_3 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_3))],
  levels(sce0$cluster_3))
sce0$colours$cluster_colours_3 <- cluster_colours_3[sce0$cluster_3]
# number of cluster = 8
set.seed(4759)
snn_gr_4 <- buildSNNGraph(sce0, use.dimred = "corrected", k=3)
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
  p5 + p6 + p7 + p8 +
  p9 + p10 + p11 + p12 +
  p13 + p14 + p15 + p16 +
  p17 + p18 + p19 + p20 +
  p21 + p22 + p23 + p24 +
  p25 + p26 + p27 + p28 +
  plot_layout(ncol = 4, guides = "collect")

# comment:
# we expect there are at least 3 groups in this subset: thymus.s1, thymus.s2, thymus.s3
# we hope to separate: thymus.s1 and thymus.s2, but increasing the number of cluster cannot achieve this purpose, but only further divide cells in group thymus.s3
# so it seems we can only compare thymus.s1+thymus.s2 and thymus.s3
# if go for 3 clusters, cluster 2 only has 8 cells, not useful at all
# although when we go for 4 clusters, it seems only subdivide thymus.s3 cells into 3, and cluster 3 and 4 seems to be the same based on the composition of experimental groups,
# cluster in between S1 and S2 could be of interest (novel ?)
# thus I decide to try with 4 clusters (but maybe merge cluster 3 and 4)
3

# selected the best `number of cluster` (colours) and replace the ones saved in SCE
sce$cluster <- factor(clusters_3$membership)
cluster_colours <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce$cluster))],
  levels(sce$cluster))
sce$colours$cluster_colours <- cluster_colours[sce$cluster]



###################################
# (M1) raw unique
#
# cluster 1 (i.e. mostly thymus S3, with S1 and S2) vs 2 (i.e. mostly thymus S1, S2) vs 3 (i.e. mostly thymus S3, with almost only S2, set 1) vs 4 (i.e. mostly thymus S3, with almost only S2, set 2)

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
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_1_vs_2_vs_3_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_1_vs_2_vs_3_vs_4"), recursive = TRUE)

vs_pair <- c("1", "2", "3", "4")

message("Writing 'uniquely_up (cluster_1_vs_2_vs_3_vs_4)' marker genes to file.")
for (n in names(uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_1_vs_2_vs_3_vs_4",
      paste0("cluster_",
             vs_pair[which(names(uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(uniquely_up) %in% n)][1],
             "_vs_",
             vs_pair[-which(names(uniquely_up) %in% n)][2],
             "_vs_",
             vs_pair[-which(names(uniquely_up) %in% n)][3],
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
# look at cluster 1 (i.e. mostly thymus S3, with S1 and S2)
chosen <- "1"
cluster1_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly thymus S3, with S1 and S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster1_uniquely_up <- cluster1_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster1_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster1_uniquely_up_noiseR <- cluster1_uniquely_up[setdiff(rownames(cluster1_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster1_uniquely_up_noiseR) %in% "CD4"),
       cluster1_uniquely_up_noiseR[which(rownames(cluster1_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster1_uniquely_up_noiseR) %in% "KLRB1"),
       cluster1_uniquely_up_noiseR[which(rownames(cluster1_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster1_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$CD4,
    sce$CD161,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "CD4",
    "CD161",
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

##########################################
# look at cluster 2 (i.e. mostly thymus S1, S2)
chosen <- "2"
cluster2_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly thymus S1, S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster2_uniquely_up <- cluster2_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster2_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster2_uniquely_up_noiseR <- cluster2_uniquely_up[setdiff(rownames(cluster2_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster2_uniquely_up_noiseR) %in% "CD4"),
       cluster2_uniquely_up_noiseR[which(rownames(cluster2_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster2_uniquely_up_noiseR) %in% "KLRB1"),
       cluster2_uniquely_up_noiseR[which(rownames(cluster2_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster2_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$CD4,
    sce$CD161,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "CD4",
    "CD161",
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
# look at cluster 3 (i.e. mostly thymus S3, with almost only S2, set 1)
chosen <- "3"
cluster3_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly thymus S3, with almost only S2, set 1)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster3_uniquely_up <- cluster3_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster3_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster3_uniquely_up_noiseR <- cluster3_uniquely_up[setdiff(rownames(cluster3_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster3_uniquely_up_noiseR) %in% "CD4"),
       cluster3_uniquely_up_noiseR[which(rownames(cluster3_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster3_uniquely_up_noiseR) %in% "KLRB1"),
       cluster3_uniquely_up_noiseR[which(rownames(cluster3_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster3_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$CD4,
    sce$CD161,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "CD4",
    "CD161",
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

##########################################
# look at cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2)
chosen <- "4"
cluster4_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly thymus S3, with almost only S2, set 2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster4_uniquely_up <- cluster4_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster4_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster4_uniquely_up_noiseR <- cluster4_uniquely_up[setdiff(rownames(cluster4_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster4_uniquely_up_noiseR) %in% "CD4"),
       cluster4_uniquely_up_noiseR[which(rownames(cluster4_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster4_uniquely_up_noiseR) %in% "KLRB1"),
       cluster4_uniquely_up_noiseR[which(rownames(cluster4_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only
best_set <- cluster4_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$CD4,
    sce$CD161,
    sce$stage,
    sce$tissue,
    sce$donor,
    sce$group,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "CD4",
    "CD161",
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

# cluster 1 (i.e. mostly thymus S3, with S1 and S2
# cluster 2 (i.e. mostly thymus S1, S2
# cluster 3 (i.e. mostly thymus S3, with almost only S2, set 1
# cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2

# COMMENT: cluster 1 and 3 don't show any unique marker, while marker of cluster 4 is only a gene family
# for cluster 1, it is the point of interest and should be left untouched, may show if look for pairwise unique, rather than global unique
# for cluster 2, it is kind of expected that it got tonnes of unique genes when compare to the rest (that are all thymus.s3 cells); this also goes in line with the minibulk
# no marker for cluster 3, and cluster 4 only have NPIB family gene as marker; should be combined if DE test tell they have no difference

# DE test to be performed:

# fx: show how thymus.s3 differed from S1-S2-mix
# 1 vs 2
# 3 vs 2
# 4 vs 2

# fx: show if 3 and 4 are different
# 3 vs 4

# fx: if yes, compare to the remaining thymus.s3
# 1 vs 3 (x)
# 1 vs 4 (o) ?

# fx: if no, combine and compare to the rest
# 1 vs 3_4 (o) ?
# 2 vs 3_4 (o)


















###################################
# (M2) selected pairwise comparison

#########
# A vs B
#########

################################################################################
# cluster 1 (i.e. mostly thymus S3, with S1 and S2) vs 2 (i.e. mostly thymus S1, S2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "2"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs1 <- factor(ifelse(cp$cluster == 1, "A", "B"))

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
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_1_vs_2.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_1_vs_2"), recursive = TRUE)

vs_pair <- c("1", "2")

message("Writing 'uniquely_up (cluster_1_vs_2)' marker genes to file.")
for (n in names(vs1_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_1_vs_2",
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
# look at cluster-group A / cluster 1 (i.e. mostly thymus S3, with S1 and S2)
chosen <- "A"
A_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; mostly thymus S3, with S1 and S2)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs1",
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
# look at cluster-group B / cluster 2 (i.e. mostly thymus S1, S2)
chosen <- "B"
B_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; mostly thymus S1, S2)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs1",
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
# cluster 2 (i.e. mostly thymus S1, S2) vs 3 (i.e. mostly thymus S3, with almost only S2, set 1)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs2 <- factor(ifelse(cp$cluster == 2, "C", "D"))

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
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_2_vs_3.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_2_vs_3"), recursive = TRUE)

vs_pair <- c("2", "3")

message("Writing 'uniquely_up (cluster_2_vs_3)' marker genes to file.")
for (n in names(vs2_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_2_vs_3",
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
# look at cluster-group C / cluster 2 (i.e. mostly thymus S1, S2)
chosen <- "C"
C_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2;  mostly thymus S1, S2)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs2",
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
# look at cluster-group D / cluster 3 (i.e. mostly thymus S3, with almost only S2, set 1)
chosen <- "D"
D_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly thymus S3, with almost only S2, set 1)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs2",
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

################################################################################
# cluster 2 (i.e. mostly thymus S1, S2) vs 4 (i.e. mostly thymus S3, with almost only S2, set 2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs3 <- factor(ifelse(cp$cluster == 2, "E", "F"))

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
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_2_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_2_vs_4"), recursive = TRUE)

vs_pair <- c("2", "4")

message("Writing 'uniquely_up (cluster_2_vs_4)' marker genes to file.")
for (n in names(vs3_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_2_vs_4",
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

##########################################################
# look at cluster-group E / cluster 2 (i.e. mostly thymus S1, S2)
chosen <- "E"
E_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2;  mostly thymus S1, S2)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs3",
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

##########################################################
# look at cluster-group F / cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2)
chosen <- "F"
F_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly thymus S3, with almost only S2, set 2)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs3",
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
# cluster 3 (i.e. mostly thymus S3, with almost only S2, set 1) vs 4 (i.e. mostly thymus S3, with almost only S2, set 2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs4 <- factor(ifelse(cp$cluster == 3, "G", "H"))

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
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_3_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_3_vs_4"), recursive = TRUE)

vs_pair <- c("3", "4")

message("Writing 'uniquely_up (cluster_3_vs_4)' marker genes to file.")
for (n in names(vs4_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_3_vs_4",
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
# look at cluster-group G / cluster 3 (i.e. mostly thymus S1, S2)
chosen <- "G"
G_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3;  mostly thymus S3, with almost only S2, set 1)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs4",
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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
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
# look at cluster-group H / cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2)
chosen <- "H"
H_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly thymus S3, with almost only S2, set 2)"

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
    cp$plate_number),
  colour_columns_by = c(
    "vs4",
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















#########
# I vs J
#########

################################################################################
# cluster 1 (i.e. mostly thymus S3, with S1 and S2) vs 3 (i.e. mostly thymus S3, with almost only S2, set 1)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs5 <- factor(ifelse(cp$cluster == 1, "I", "J"))

# set vs colours
vs5_colours <- setNames(
  palette.colors(nlevels(cp$vs5), "Set1"),
  levels(cp$vs5))
cp$colours$vs5_colours <- vs5_colours[cp$vs5]

# find unique DE ./. cluster-groups
vs5_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs5,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs5_uniquely_up,
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_1_vs_3.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_1_vs_3"), recursive = TRUE)

vs_pair <- c("1", "3")

message("Writing 'uniquely_up (cluster_1_vs_3)' marker genes to file.")
for (n in names(vs5_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_1_vs_3",
      paste0("cluster_",
             vs_pair[which(names(vs5_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs5_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs5_uniquely_up[[n]] %>%
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
# look at cluster-group I / cluster 1 (i.e. mostly thymus S3, with S1 and S2)
chosen <- "I"
I_uniquely_up <- vs5_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1;  mostly thymus S3, with S1 and S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# I_uniquely_up_pcg <- I_uniquely_up[intersect(protein_coding_gene_set, rownames(I_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
I_uniquely_up_noiseR <- I_uniquely_up[setdiff(rownames(I_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(I_uniquely_up_noiseR) %in% "CD4"),
       I_uniquely_up_noiseR[which(rownames(I_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(I_uniquely_up_noiseR) %in% "KLRB1"),
       I_uniquely_up_noiseR[which(rownames(I_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- I_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs5,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs5",
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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs5 = vs5_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################################
# look at cluster-group J / cluster 3 (i.e. mostly thymus S3, with almost only S2, set 1)
chosen <- "J"
J_uniquely_up <- vs5_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly thymus S3, with almost only S2, set 1)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# J_uniquely_up_pcg <- J_uniquely_up[intersect(protein_coding_gene_set, rownames(J_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
J_uniquely_up_noiseR <- J_uniquely_up[setdiff(rownames(J_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(J_uniquely_up_noiseR) %in% "CD4"),
       J_uniquely_up_noiseR[which(rownames(J_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(J_uniquely_up_noiseR) %in% "KLRB1"),
       J_uniquely_up_noiseR[which(rownames(J_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- J_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs5,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs5",
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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs5 = vs5_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



#########
# K vs L
#########

################################################################################
# cluster 1 (i.e. mostly thymus S3, with S1 and S2) vs 4 (i.e. mostly thymus S3, with almost only S2, set 2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs6 <- factor(ifelse(cp$cluster == 1, "K", "L"))

# set vs colours
vs6_colours <- setNames(
  palette.colors(nlevels(cp$vs6), "Set1"),
  levels(cp$vs6))
cp$colours$vs6_colours <- vs6_colours[cp$vs6]

# find unique DE ./. cluster-groups
vs6_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs6,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs6_uniquely_up,
  here("data", "marker_genes", "thymus_only", "C094_Pellicci.uniquely_up.cluster_1_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "thymus_only", "uniquely_up", "cluster_1_vs_4"), recursive = TRUE)

vs_pair <- c("1", "4")

message("Writing 'uniquely_up (cluster_1_vs_4)' marker genes to file.")
for (n in names(vs6_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "thymus_only",
      "uniquely_up",
      "cluster_1_vs_4",
      paste0("cluster_",
             vs_pair[which(names(vs6_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs6_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs6_uniquely_up[[n]] %>%
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
# look at cluster-group K / cluster 1 (i.e. mostly thymus S3, with S1 and S2)
chosen <- "K"
K_uniquely_up <- vs6_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1;  mostly thymus S3, with S1 and S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# K_uniquely_up_pcg <- K_uniquely_up[intersect(protein_coding_gene_set, rownames(K_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
K_uniquely_up_noiseR <- K_uniquely_up[setdiff(rownames(K_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(K_uniquely_up_noiseR) %in% "CD4"),
       K_uniquely_up_noiseR[which(rownames(K_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(K_uniquely_up_noiseR) %in% "KLRB1"),
       K_uniquely_up_noiseR[which(rownames(K_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- K_uniquely_up_noiseR[1:25, ]


# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs6,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs6",
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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs6 = vs6_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##########################################################
# look at cluster-group L / cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2)
chosen <- "L"
L_uniquely_up <- vs6_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly thymus S3, with almost only S2, set 2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# L_uniquely_up_pcg <- L_uniquely_up[intersect(protein_coding_gene_set, rownames(L_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
L_uniquely_up_noiseR <- L_uniquely_up[setdiff(rownames(L_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(L_uniquely_up_noiseR) %in% "CD4"),
       L_uniquely_up_noiseR[which(rownames(L_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(L_uniquely_up_noiseR) %in% "KLRB1"),
       L_uniquely_up_noiseR[which(rownames(L_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- L_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs6,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs6",
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
    vs6 = vs6_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)







#########
# M vs M
#########

##########################################################################
# cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2) vs cluster 1_3 (i.e. mostly thymus S3, S1S2.alike.and.S2.alike-mix)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs7 <- factor(ifelse(cp$cluster == 4, "M", "N"))

# set vs colours
vs7_colours <- setNames(
  palette.colors(nlevels(cp$vs7), "Set1"),
  levels(cp$vs7))
cp$colours$vs7_colours <- vs7_colours[cp$vs7]

# find unique DE ./. cluster-groups
vs7_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs7,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs7_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_4_vs_1_3.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_4_vs_1_3"), recursive = TRUE)

vs_pair <- c("4", "1_3")

message("Writing 'uniquely_up (cluster_4_vs_1_3)' marker genes to file.")
for (n in names(vs7_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_4_vs_1_3",
      paste0("cluster_",
             vs_pair[which(names(vs7_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs7_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs7_uniquely_up[[n]] %>%
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
# look at cluster-group M / cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2)
chosen <- "M"
M_uniquely_up <- vs7_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly thymus S3, with almost only S2, set 2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# M_uniquely_up_pcg <- M_uniquely_up[intersect(protein_coding_gene_set, rownames(M_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
M_uniquely_up_noiseR <- M_uniquely_up[setdiff(rownames(M_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(M_uniquely_up_noiseR) %in% "CD4"),
       M_uniquely_up_noiseR[which(rownames(M_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(M_uniquely_up_noiseR) %in% "KLRB1"),
       M_uniquely_up_noiseR[which(rownames(M_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- M_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs7,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs7",
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
    vs7 = vs7_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group N / cluster 1_3 (i.e. mostly thymus S3, S1S2.alike.and.S2.alike-mix)
chosen <- "N"
N_uniquely_up <- vs7_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1_3; mostly thymus S3, S1S2.alike.and.S2.alike-mix)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# N_uniquely_up_pcg <- N_uniquely_up[intersect(protein_coding_gene_set, rownames(N_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
N_uniquely_up_noiseR <- N_uniquely_up[setdiff(rownames(N_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(N_uniquely_up_noiseR) %in% "CD4"),
       N_uniquely_up_noiseR[which(rownames(N_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(N_uniquely_up_noiseR) %in% "KLRB1"),
       N_uniquely_up_noiseR[which(rownames(N_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- N_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs7,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs7",
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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs7 = vs7_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



#########
# O vs P
#########

##########################################################################
# cluster 1 (i.e. mostly thymus S3, with S1 and S2) vs cluster 3_4 (i.e. mostly thymus S3, with almost only S2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs8 <- factor(ifelse(cp$cluster == 1, "O", "P"))

# set vs colours
vs8_colours <- setNames(
  palette.colors(nlevels(cp$vs8), "Set1"),
  levels(cp$vs8))
cp$colours$vs8_colours <- vs8_colours[cp$vs8]

# find unique DE ./. cluster-groups
vs8_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs8,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs8_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_1_vs_3_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_1_vs_3_4"), recursive = TRUE)

vs_pair <- c("1", "3_4")

message("Writing 'uniquely_up (cluster_1_vs_3_4)' marker genes to file.")
for (n in names(vs8_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_1_vs_3_4",
      paste0("cluster_",
             vs_pair[which(names(vs8_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs8_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs8_uniquely_up[[n]] %>%
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
# look at cluster-group O / cluster 1 (i.e. mostly thymus S3, with S1 and S2)
chosen <- "O"
O_uniquely_up <- vs8_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; mostly thymus S3, with S1 and S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# O_uniquely_up_pcg <- O_uniquely_up[intersect(protein_coding_gene_set, rownames(O_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
O_uniquely_up_noiseR <- O_uniquely_up[setdiff(rownames(O_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(O_uniquely_up_noiseR) %in% "CD4"),
       O_uniquely_up_noiseR[which(rownames(O_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(O_uniquely_up_noiseR) %in% "KLRB1"),
       O_uniquely_up_noiseR[which(rownames(O_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- O_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs8,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs8",
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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("Yes")),
  #   row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - \n",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; \n",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    # Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs8 = vs8_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group P / cluster 3_4 (i.e. mostly thymus S3, with almost only S2)
chosen <- "P"
P_uniquely_up <- vs8_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3_4; mostly thymus S3, with almost only S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# P_uniquely_up_pcg <- P_uniquely_up[intersect(protein_coding_gene_set, rownames(P_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
P_uniquely_up_noiseR <- P_uniquely_up[setdiff(rownames(P_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(P_uniquely_up_noiseR) %in% "CD4"),
       P_uniquely_up_noiseR[which(rownames(P_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(P_uniquely_up_noiseR) %in% "KLRB1"),
       P_uniquely_up_noiseR[which(rownames(P_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- P_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs8,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs8",
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
    vs8 = vs8_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)




#########
# Q vs R
#########

##########################################################################
# cluster 2 (i.e. mostly thymus S1, S2) vs cluster 3_4 (i.e. mostly thymus S3, with almost only S2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs9 <- factor(ifelse(cp$cluster == 2, "Q", "R"))

# set vs colours
vs9_colours <- setNames(
  palette.colors(nlevels(cp$vs9), "Set1"),
  levels(cp$vs9))
cp$colours$vs9_colours <- vs9_colours[cp$vs9]

# find unique DE ./. cluster-groups
vs9_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs9,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs9_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_2_vs_3_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_2_vs_3_4"), recursive = TRUE)

vs_pair <- c("2", "3_4")

message("Writing 'uniquely_up (cluster_2_vs_3_4)' marker genes to file.")
for (n in names(vs9_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_2_vs_3_4",
      paste0("cluster_",
             vs_pair[which(names(vs9_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs9_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs9_uniquely_up[[n]] %>%
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
# look at cluster-group Q / cluster 2 (i.e. mostly thymus S1, S2)
chosen <- "Q"
Q_uniquely_up <- vs9_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; mostly thymus S1, S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# Q_uniquely_up_pcg <- Q_uniquely_up[intersect(protein_coding_gene_set, rownames(Q_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
Q_uniquely_up_noiseR <- Q_uniquely_up[setdiff(rownames(Q_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(Q_uniquely_up_noiseR) %in% "CD4"),
       Q_uniquely_up_noiseR[which(rownames(Q_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(Q_uniquely_up_noiseR) %in% "KLRB1"),
       Q_uniquely_up_noiseR[which(rownames(Q_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- Q_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs9,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs9",
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
    vs9 = vs9_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group R / cluster 3_4 (i.e. mostly thymus S3, with almost only S2)
chosen <- "R"
R_uniquely_up <- vs9_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3_4; mostly thymus S3, with almost only S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# R_uniquely_up_pcg <- R_uniquely_up[intersect(protein_coding_gene_set, rownames(R_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
R_uniquely_up_noiseR <- R_uniquely_up[setdiff(rownames(R_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(R_uniquely_up_noiseR) %in% "CD4"),
       R_uniquely_up_noiseR[which(rownames(R_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(R_uniquely_up_noiseR) %in% "KLRB1"),
       R_uniquely_up_noiseR[which(rownames(R_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- R_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs9,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs9",
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
    vs9 = vs9_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)


# cluster 1 (i.e. mostly thymus S3, with S1 and S2
# cluster 2 (i.e. mostly thymus S1, S2
# cluster 3 (i.e. mostly thymus S3, with almost only S2, set 1
# cluster 4 (i.e. mostly thymus S3, with almost only S2, set 2
# cluster 3_4 (i.e. mostly thymus S3, with almost only S2



# fx: show how thymus.s3 differed from S1-S2-mix

# 1 vs 2 (A vs B)
# cluster 1 (mostly thymus S3, with S1 and S2 >>> quite some strong markers unique in thymus.s3.alike.s1s2
# cluster 2 (mostly thymus S1, S2 >>> lots of markers expressed more frequently in thymus.s1.s2.mix, but mostly not in thymus.s3.alike.s1s2
# COMMENT: thymus.s1.s2.mix and thymus.s3.alike.s1s2 are different from each other

# 2 vs 3 (C vs D)
# cluster 3 (mostly thymus S3, with almost only S2, set 1 >>> thymus.s3.alike.s2.set1 got lots of strong markers (including HLA)
# cluster 2 (mostly thymus S1, S2 >>> lots of markers expressed more frequently in thymus.s1.s2.mix, but mostly not in thymus.s3.alike.s2.set1
# COMMENT: thymus.s1.s2.mix and thymus.s3.alike.s2.set1 are also different from each other

# 2 vs 4 (E vs F)
# cluster 2 (mostly thymus S1, S2 >>> lot of more frequently expressed markers in thymus.s1.s2.mix
# cluster 4 (mostly thymus S3, with almost only S2 >>> also have lot of markers more frequently expressed in thymus.s3.alike.s2.set2
# COMMENT: thymus.s1.s2.mix and thymus.s3.alike.s2.set2 are also different from each other



# fx: show if 3 and 4 are different

# 3 vs 4 (G vs H)
# cluster 3 (mostly thymus S3, with almost only S2, set 1 >>> no unique markers detected for thymus.s3.alike.s2.set1 (though visually, may have more frequently expressed genes; sample size not enough ?)
# cluster 4 (mostly thymus S3, with almost only S2, set 2 >>> except NPIPB, no unique markers detected for thymus.s3.alike.s2.set2 neither (though visually, may have more frequently expressed genes; sample size not enough)
# COMMENT: although the 2 sets maybe different based on the frequency of expression in heatmap; based on statistics only, there should have no sig difference between cluster 3 and 4; should be combined !



# fx: if yes, compare to the remaining thymus.s3

# 1 vs 3 (I vs J)
# cluster 1 (mostly thymus S3, with S1 and S2 >>> visually yes; statistically no
# cluster 3 (mostly thymus S3, with almost only S2, set 1 >>> visually yes; statistically no
# COMMENT: no significant difference between thymus.s3.alike.s1s2 and thymus.s3.alike.s2.set1

# 1 vs 4 (K vs L)
# cluster 1 (mostly thymus S3, with S1 and S2 >>> visually yes; statistically no
# cluster 4 (mostly thymus S3, with almost only S2, set 2 >>> besides NPIPB family, CCL5, MRPS21, SIT1 and SDHB shown to be sig different
# COMMENT: thymus.s3.alike.s1s2 does not have anything unique, but thymus.s3.alike.s2.set2 does show several marker unique to the cluster
# although 3 and 4 are not sig differed, when reference to 1, 4 with unique markers, which may explain why 4 was in different cluster with 3
# it may be clearer if there are more cells > try combing 1 and 3 vs 4



# fx: if no, combine and compare to the rest

# 4 vs 1_3 (M vs N)
# cluster 4 (mostly thymus S3, with almost only S2, set 2 >>> visually yes; statistically no
# cluster 1_3 (mostly thymus S3, S1S2.alike.and.S2.alike-mix) >>> besides NPIPB, there seems to be AC009022.1 as a unique marker
# remark: AC009022.1 = a LncRNA enhances cells proliferation, migration, and invasion in colorectal cancer (J Cell Biochem. 2020 Feb;121(2):1934-1944. doi: 10.1002/jcb.29428)
# COMMENT: it doesn't help to separate, but even further mask the unique marker in 4, though revealing AC009022.1 (not invasive, but proliferating and migration, relevant ?)

# 1 vs 3_4 (O vs P)
# cluster 1 (mostly thymus S3, with S1 and S2 >>> visually yes; statistically no unique marker found on thymus.s3.alike.s1s2
# cluster 3_4 (mostly thymus S3, with almost only S2 >>> lots of more frequently expressed markers (even excluding NPIPB)
# COMMENT: cluster 3 and 4 together, they does have number of unique markers in common that separate them out from cluster 1 (including AC009022.1)

# 2 vs 3_4 (Q vs R)
# cluster 2 (mostly thymus S1, S2 >>> lot of clear markers detected
# cluster 3_4 (mostly thymus S3, with almost only S2 >>> lot of clear markers detected
# COMMENT: cluster 3_4 is also strongly differed from cluster 2



# CONCLUSION:
#
# consistent with minibulk, thymus.s1.s2.mix is strongly different from any component of the thymus.s3.bulk
#
# by pairwise DE detection, cluster 3 and 4 shows no difference and should be grouped
# then when comparing the grouped 3_4 to 1, we see a number of markers that are sig more frequently express in 3_4 only
# on one hand, as 3_4 does contain a small amount of S2, that may explain why the gene expression of those markers are not consistently expressed throughout all cells in the cluster
# on the other hand, it may also indicate cluster 1 is a transition or interim cell group (new subgroup !) found with thymus that transit from S1-S2 to S3, with increasing expression of these markers found on 3_4/S3 only in thymus
#
# also one important point to note is that
# although gene family (like the NPIPB here) appeared as sig marker may not be optimal (as the "variance" detected could be due to mapping artifact due to their sequence conservation of this gene family 3')
# but there are tonnes of gene family in human and I believe number of them have 3' sequence conservancy,
# if it is a common problem true for all gene family, then they should all appear as DE, but not only NPIPB - which show a strong and clear expression in cluster 4 only, and this must not be neglected ! (novelty !)
# but one point that is true is, which marker is the actual DE is not definite and I agree that we need to be cautious !























































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


