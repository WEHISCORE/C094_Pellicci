# EDA: annotate (S3 only)
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

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.single-cell.merged.S3_only.SCE.rds"))

# pre-create directories for saving export, or error (dir not exists)
dir.create(here("data", "marker_genes", "S3_only"), recursive = TRUE)
dir.create(here("output", "marker_genes", "S3_only"), recursive = TRUE)

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
# number of cluster = 6 (default)
set.seed(4759)
snn_gr_1 <- buildSNNGraph(sce0, use.dimred = "corrected", k=10)
clusters_1 <- igraph::cluster_louvain(snn_gr_1)
sce0$cluster_1 <- factor(clusters_1$membership)
cluster_colours_1 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_1))],
  levels(sce0$cluster_1))
sce0$colours$cluster_colours_1 <- cluster_colours_1[sce0$cluster_1]
# number of cluster = 5
set.seed(4759)
snn_gr_2 <- buildSNNGraph(sce0, use.dimred = "corrected", k=15)
clusters_2 <- igraph::cluster_louvain(snn_gr_2)
sce0$cluster_2 <- factor(clusters_2$membership)
cluster_colours_2 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_2))],
  levels(sce0$cluster_2))
sce0$colours$cluster_colours_2 <- cluster_colours_2[sce0$cluster_2]
# number of cluster = 4
set.seed(4759)
snn_gr_3 <- buildSNNGraph(sce0, use.dimred = "corrected", k=20)
clusters_3 <- igraph::cluster_louvain(snn_gr_3)
sce0$cluster_3 <- factor(clusters_3$membership)
cluster_colours_3 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_3))],
  levels(sce0$cluster_3))
sce0$colours$cluster_colours_3 <- cluster_colours_3[sce0$cluster_3]
# number of cluster = 3
set.seed(4759)
snn_gr_4 <- buildSNNGraph(sce0, use.dimred = "corrected", k=40)
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
# we expect there are at least 2 groups in this subset: thymus.s3, blood.s3
# in default setting, we got 6 clusters, almost all clusters are dominant by cells in thymus.s3
# I cannot see how it could tell a difference unless we group cluster by ourselves based on 6 clusters
# NOTE: tried to group: i.e. cluster 1, 2, 3 as 1 group, 4 and 5 form another group, and cluster 6 on its own > but no significant result
# as the clustering algorithm seems to be struggle in forming cluster except the "3 cluster"
# + need around 100 cells for GLM
# also for the 4 cluster setting, found IL7R that is quite clear and unique in cluster 1, so pick 4 clusters !
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
# cluster 1 (i.e.  S3-mix, set 1) vs 2 (i.e.  S3-mix, set 2) vs 3 (i.e.  S3-mix, set 3) vs 4 (i.e.  S3-mix, set 4)

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
# look at cluster 1 (i.e. S3-mix, set 1)
chosen <- "1"
cluster1_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3-mix, set 1)"

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
# look at cluster 2 (i.e. S3-mix, set 2)
chosen <- "2"
cluster2_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3-mix, set 2)"

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
# look at cluster 3 (i.e. S3-mix, set 3)
chosen <- "3"
cluster3_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3-mix, set 3)"

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
# look at cluster 4 (i.e. S3-mix, set 4)
chosen <- "4"
cluster4_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3-mix, set 3)"

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
# cluster 1: quite clear to have unique IL7R expression
# cluster 2: got a few markers (excluding NPBIB, got AC009022.1); also visually, CCL5 seems also be a marker
# cluster 3: no statistically sig marker; but viusally, SOX4 could be a candidate if pairwise compare
# cluster 4: no sig marker, but perhaps because it is at center

# thus, they can be annotated as:
# cluster 1 (i.e. S3-mix more thymus, peripheral IL7R driven
# cluster 2 (i.e. S3-mix more thymus, peripheral may CCL5 driven
# cluster 3 (i.e. S3-mix more thymus, peripheral may SOX4 drive
# cluster 4 (i.e. S3-mix more thymus, center






###################################
# (M2) selected pairwise comparison




#########
# A vs B
#########

################################################################################
# cluster 1 (i.e. S3-mix, set 1) vs 2 (i.e. S3-mix, set 2)

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
  here("data", "marker_genes", "S3_only", "C094_Pellicci.uniquely_up.cluster_1_vs_2.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S3_only", "uniquely_up", "cluster_1_vs_2"), recursive = TRUE)

vs_pair <- c("1", "2")

message("Writing 'uniquely_up (cluster_1_vs_2)' marker genes to file.")
for (n in names(vs1_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S3_only",
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
# look at cluster-group A / cluster 1 (i.e. S3-mix, set 1)
chosen <- "A"
A_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix, set 1)"

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
# look at cluster-group B / cluster 2 (i.e. S3-mix, set 2)
chosen <- "B"
B_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix, set 2)"

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
# cluster 1 (i.e. S3-mix, set 1) vs 3 (i.e. S3-mix, set 3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs2 <- factor(ifelse(cp$cluster == 1, "C", "D"))

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
  here("data", "marker_genes", "S3_only", "C094_Pellicci.uniquely_up.cluster_1_vs_3.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S3_only", "uniquely_up", "cluster_1_vs_3"), recursive = TRUE)

vs_pair <- c("1", "3")

message("Writing 'uniquely_up (cluster_1_vs_3)' marker genes to file.")
for (n in names(vs2_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S3_only",
      "uniquely_up",
      "cluster_1_vs_3",
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
# look at cluster-group C / cluster 1 (i.e. S3-mix, set 1)
chosen <- "C"
C_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix, set 1)"

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
# look at cluster-group D / cluster 3 (i.e. S3-mix, set 3)
chosen <- "D"
D_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix, set 2)"

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
# cluster 1 (i.e. S3-mix, set 1) vs 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs3 <- factor(ifelse(cp$cluster == 1, "E", "F"))

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
  here("data", "marker_genes", "S3_only", "C094_Pellicci.uniquely_up.cluster_1_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S3_only", "uniquely_up", "cluster_1_vs_4"), recursive = TRUE)

vs_pair <- c("1", "4")

message("Writing 'uniquely_up (cluster_1_vs_4)' marker genes to file.")
for (n in names(vs3_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S3_only",
      "uniquely_up",
      "cluster_1_vs_4",
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
# look at cluster-group E / cluster 1 (i.e. S3-mix, set 1)
chosen <- "E"
E_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix, set 1)"

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
# look at cluster-group F / cluster 4 (i.e. S3-mix, set 4)
chosen <- "F"
F_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

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
# cluster 2 (i.e. S3-mix, set 2) vs 3 (i.e. S3-mix, set 3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs4 <- factor(ifelse(cp$cluster == 2, "G", "H"))

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
  here("data", "marker_genes", "S3_only", "C094_Pellicci.uniquely_up.cluster_2_vs_3.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S3_only", "uniquely_up", "cluster_2_vs_3"), recursive = TRUE)

vs_pair <- c("2", "3")

message("Writing 'uniquely_up (cluster_2_vs_3)' marker genes to file.")
for (n in names(vs4_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S3_only",
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
# look at cluster-group G / cluster 2 (i.e. S3-mix, set 2)
chosen <- "G"
G_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix, set 2)"

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
# look at cluster-group H / cluster 3 (i.e. S3-mix, set 3)
chosen <- "H"
H_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

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




#########
# I vs J
#########

################################################################################
# cluster 2 (i.e. S3-mix, set 2) vs 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs5 <- factor(ifelse(cp$cluster == 2, "I", "J"))

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
  here("data", "marker_genes", "S3_only", "C094_Pellicci.uniquely_up.cluster_2_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S3_only", "uniquely_up", "cluster_2_vs_4"), recursive = TRUE)

vs_pair <- c("2", "4")

message("Writing 'uniquely_up (cluster_2_vs_4)' marker genes to file.")
for (n in names(vs5_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S3_only",
      "uniquely_up",
      "cluster_2_vs_4",
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
# look at cluster-group I / cluster 2 (i.e. S3-mix, set 2)
chosen <- "I"
I_uniquely_up <- vs5_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix, set 2)"

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
# look at cluster-group J / cluster 4 (i.e. S3-mix, set 4)
chosen <- "J"
J_uniquely_up <- vs5_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

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
# cluster 3 (i.e. S3-mix, set 3) vs 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs6 <- factor(ifelse(cp$cluster == 3, "K", "L"))

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
  here("data", "marker_genes", "S3_only", "C094_Pellicci.uniquely_up.cluster_3_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "S3_only", "uniquely_up", "cluster_3_vs_4"), recursive = TRUE)

vs_pair <- c("3", "4")

message("Writing 'uniquely_up (cluster_3_vs_4)' marker genes to file.")
for (n in names(vs6_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "S3_only",
      "uniquely_up",
      "cluster_3_vs_4",
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
# look at cluster-group K / cluster 3 (i.e. S3-mix, set 3)
chosen <- "K"
K_uniquely_up <- vs6_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; S3-mix, set 3)"

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

##########################################################
# look at cluster-group L / cluster 4 (i.e. S3-mix, set 4)
chosen <- "L"
L_uniquely_up <- vs6_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

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







#########
# M vs N
#########

##########################################################################
# cluster 1_2 (i.e. S3-mix, set 1_2) vs cluster 3 (i.e. S3-mix, set 3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "2" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs7 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 2, "M", "N"))

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
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_1_2_vs_3.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_1_2_vs_3"), recursive = TRUE)

vs_pair <- c("1_2", "3")

message("Writing 'uniquely_up (cluster_1_2_vs_3)' marker genes to file.")
for (n in names(vs7_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_1_2_vs_3",
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
# look at cluster-group M / cluster 1_2 (i.e. S3-mix, set 1_2)
chosen <- "M"
M_uniquely_up <- vs7_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1_2; S3-mix, set 1_2)"

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

##############################################################
# look at cluster-group N / cluster 3 (i.e. S3-mix, set 3)
chosen <- "N"
N_uniquely_up <- vs7_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; S3-mix, set 3)"

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
  # row.names = rownames(best_set)),
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
# cluster 2_3 (i.e. S3-mix, set 2_3) vs cluster 1 (i.e. S3-mix, set 1)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "2" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs8 <- factor(ifelse(cp$cluster == 2 | cp$cluster == 3, "O", "P"))

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
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_2_3_vs_1.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_2_3_vs_1"), recursive = TRUE)

vs_pair <- c("2_3", "1")

message("Writing 'uniquely_up (cluster_2_3_vs_1)' marker genes to file.")
for (n in names(vs8_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_2_3_vs_1",
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
# look at cluster-group O / cluster 2_3 (i.e. S3-mix, set 2_3)
chosen <- "O"
O_uniquely_up <- vs8_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2_3; S3-mix, set 2_3)"

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

##############################################################
# look at cluster-group P / cluster 1 (i.e. S3-mix, set 1)
chosen <- "P"
P_uniquely_up <- vs8_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix, set 1)"

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
# cluster 1_3 (i.e. S3-mix, set 1_3) vs cluster 2 (i.e. S3-mix, set 2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "2" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs9 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 3, "Q", "R"))

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
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_1_3_vs_2.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_1_3_vs_2"), recursive = TRUE)

vs_pair <- c("1_3", "2")

message("Writing 'uniquely_up (cluster_1_3_vs_2)' marker genes to file.")
for (n in names(vs9_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_1_3_vs_2",
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
# look at cluster-group Q / cluster 1_3 (i.e. S3-mix, set 1_3)
chosen <- "Q"
Q_uniquely_up <- vs9_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1_3; S3-mix, set 1_3)"

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
# look at cluster-group R / cluster 2 (i.e. S3-mix, set 2)
chosen <- "R"
R_uniquely_up <- vs9_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix, set 2)"

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







#########
# S vs T
#########

##########################################################################
# cluster 1_2 (i.e. S3-mix, set 1_2) vs cluster 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "2" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs10 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 2, "S", "T"))

# set vs colours
vs10_colours <- setNames(
  palette.colors(nlevels(cp$vs10), "Set1"),
  levels(cp$vs10))
cp$colours$vs10_colours <- vs10_colours[cp$vs10]

# find unique DE ./. cluster-groups
vs10_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs10,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs10_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_1_2_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_1_2_vs_4"), recursive = TRUE)

vs_pair <- c("1_2", "4")

message("Writing 'uniquely_up (cluster_1_2_vs_4)' marker genes to file.")
for (n in names(vs10_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_1_2_vs_4",
      paste0("cluster_",
             vs_pair[which(names(vs10_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs10_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs10_uniquely_up[[n]] %>%
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
# look at cluster-group S / cluster 1_2 (i.e. S3-mix, set 1_2)
chosen <- "S"
S_uniquely_up <- vs10_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1_2; S3-mix, set 1_2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# S_uniquely_up_pcg <- S_uniquely_up[intersect(protein_coding_gene_set, rownames(S_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
S_uniquely_up_noiseR <- S_uniquely_up[setdiff(rownames(S_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(S_uniquely_up_noiseR) %in% "CD4"),
       S_uniquely_up_noiseR[which(rownames(S_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(S_uniquely_up_noiseR) %in% "KLRB1"),
       S_uniquely_up_noiseR[which(rownames(S_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- S_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs10,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs10",
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
    vs10 = vs10_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group T / cluster 4 (i.e. S3-mix, set 4)
chosen <- "T"
T_uniquely_up <- vs10_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# T_uniquely_up_pcg <- T_uniquely_up[intersect(protein_coding_gene_set, rownames(T_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
T_uniquely_up_noiseR <- T_uniquely_up[setdiff(rownames(T_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(T_uniquely_up_noiseR) %in% "CD4"),
       T_uniquely_up_noiseR[which(rownames(T_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(T_uniquely_up_noiseR) %in% "KLRB1"),
       T_uniquely_up_noiseR[which(rownames(T_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- T_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs10,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs10",
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
    vs10 = vs10_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)







#########
# U vs V
#########

##########################################################################
# cluster 2_3 (i.e. S3-mix, set 2_3) vs cluster 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs11 <- factor(ifelse(cp$cluster == 2 | cp$cluster == 3, "U", "V"))

# set vs colours
vs11_colours <- setNames(
  palette.colors(nlevels(cp$vs11), "Set1"),
  levels(cp$vs11))
cp$colours$vs11_colours <- vs11_colours[cp$vs11]

# find unique DE ./. cluster-groups
vs11_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs11,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs11_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_2_3_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_2_3_vs_4"), recursive = TRUE)

vs_pair <- c("2_3", "4")

message("Writing 'uniquely_up (cluster_2_3_vs_4)' marker genes to file.")
for (n in names(vs11_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_2_3_vs_4",
      paste0("cluster_",
             vs_pair[which(names(vs11_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs11_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs11_uniquely_up[[n]] %>%
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
# look at cluster-group U / cluster 2_3 (i.e. S3-mix, set 2_3)
chosen <- "U"
U_uniquely_up <- vs11_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2_3; S3-mix, set 2_3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# U_uniquely_up_pcg <- U_uniquely_up[intersect(protein_coding_gene_set, rownames(U_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
U_uniquely_up_noiseR <- U_uniquely_up[setdiff(rownames(U_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(U_uniquely_up_noiseR) %in% "CD4"),
       U_uniquely_up_noiseR[which(rownames(U_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(U_uniquely_up_noiseR) %in% "KLRB1"),
       U_uniquely_up_noiseR[which(rownames(U_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- U_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs11,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs11",
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
    vs11 = vs11_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group V / cluster 4 (i.e. S3-mix, set 4)
chosen <- "V"
V_uniquely_up <- vs11_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# V_uniquely_up_pcg <- V_uniquely_up[intersect(protein_coding_gene_set, rownames(V_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
V_uniquely_up_noiseR <- V_uniquely_up[setdiff(rownames(V_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(V_uniquely_up_noiseR) %in% "CD4"),
       V_uniquely_up_noiseR[which(rownames(V_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(V_uniquely_up_noiseR) %in% "KLRB1"),
       V_uniquely_up_noiseR[which(rownames(V_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- V_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs11,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs11",
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
    vs11 = vs11_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)







#########
# W vs X
#########

##########################################################################
# cluster 1_3 (i.e. S3-mix, set 1_3) vs cluster 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs12 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 3, "W", "X"))

# set vs colours
vs12_colours <- setNames(
  palette.colors(nlevels(cp$vs12), "Set1"),
  levels(cp$vs12))
cp$colours$vs12_colours <- vs12_colours[cp$vs12]

# find unique DE ./. cluster-groups
vs12_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs12,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs12_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_1_3_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_1_3_vs_4"), recursive = TRUE)

vs_pair <- c("1_3", "4")

message("Writing 'uniquely_up (cluster_1_3_vs_4)' marker genes to file.")
for (n in names(vs12_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_1_3_vs_4",
      paste0("cluster_",
             vs_pair[which(names(vs12_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs12_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs12_uniquely_up[[n]] %>%
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
# look at cluster-group W / cluster 1_3 (i.e. S3-mix, set 1_3)
chosen <- "W"
W_uniquely_up <- vs12_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1_3; S3-mix, set 1_3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# W_uniquely_up_pcg <- W_uniquely_up[intersect(protein_coding_gene_set, rownames(W_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
W_uniquely_up_noiseR <- W_uniquely_up[setdiff(rownames(W_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(W_uniquely_up_noiseR) %in% "CD4"),
       W_uniquely_up_noiseR[which(rownames(W_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(W_uniquely_up_noiseR) %in% "KLRB1"),
       W_uniquely_up_noiseR[which(rownames(W_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- W_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs12,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs12",
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
    vs12 = vs12_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group X / cluster 4 (i.e. S3-mix, set 4)
chosen <- "X"
X_uniquely_up <- vs12_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# X_uniquely_up_pcg <- X_uniquely_up[intersect(protein_coding_gene_set, rownames(X_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
X_uniquely_up_noiseR <- X_uniquely_up[setdiff(rownames(X_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(X_uniquely_up_noiseR) %in% "CD4"),
       X_uniquely_up_noiseR[which(rownames(X_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(X_uniquely_up_noiseR) %in% "KLRB1"),
       X_uniquely_up_noiseR[which(rownames(X_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- X_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs12,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs12",
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
    vs12 = vs12_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)







#########
# Y vs Z
#########

##########################################################################
# cluster 1_2_3 (i.e. S3-mix, set 1_2_3) vs cluster 4 (i.e. S3-mix, set 4)

# checkpoint
cp <- sce

# classify cluster-group for comparison
cp$vs13 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 2 | cp$cluster == 3, "Y", "Z"))

# set vs colours
vs13_colours <- setNames(
  palette.colors(nlevels(cp$vs13), "Set1"),
  levels(cp$vs13))
cp$colours$vs13_colours <- vs13_colours[cp$vs13]

# find unique DE ./. cluster-groups
vs13_uniquely_up <- findMarkers(
  cp,
  groups = cp$vs13,
  block = cp$block,
  pval.type = "all",
  direction = "up")

# export DGE lists
saveRDS(
  vs13_uniquely_up,
  here("data", "marker_genes", "whole_cell", "C094_Pellicci.uniquely_up.cluster_1_2_3_vs_4.rds"),
  compress = "xz")

dir.create(here("output", "marker_genes", "whole_cell", "uniquely_up", "cluster_1_2_3_vs_4"), recursive = TRUE)

vs_pair <- c("1_2_3", "4")

message("Writing 'uniquely_up (cluster_1_2_3_vs_4)' marker genes to file.")
for (n in names(vs13_uniquely_up)) {
  message(n)
  gzout <- gzfile(
    description = here(
      "output",
      "marker_genes",
      "whole_cell",
      "uniquely_up",
      "cluster_1_2_3_vs_4",
      paste0("cluster_",
             vs_pair[which(names(vs13_uniquely_up) %in% n)],
             "_vs_",
             vs_pair[-which(names(vs13_uniquely_up) %in% n)][1],
             ".uniquely_up.csv.gz")),
    open = "wb")
  write.table(
    x = vs13_uniquely_up[[n]] %>%
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
# look at cluster-group Y / cluster 1_2_3 (i.e. S3-mix, set 1_2_3)
chosen <- "Y"
Y_uniquely_up <- vs13_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1_2_3; S3-mix, set 1_2_3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# Y_uniquely_up_pcg <- Y_uniquely_up[intersect(protein_coding_gene_set, rownames(Y_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
Y_uniquely_up_noiseR <- Y_uniquely_up[setdiff(rownames(Y_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(Y_uniquely_up_noiseR) %in% "CD4"),
       Y_uniquely_up_noiseR[which(rownames(Y_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(Y_uniquely_up_noiseR) %in% "KLRB1"),
       Y_uniquely_up_noiseR[which(rownames(Y_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- Y_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs13,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs13",
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
    vs13 = vs13_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

##############################################################
# look at cluster-group Z / cluster 4 (i.e. S3-mix, set 4)
chosen <- "Z"
Z_uniquely_up <- vs13_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; S3-mix, set 4)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# Z_uniquely_up_pcg <- Z_uniquely_up[intersect(protein_coding_gene_set, rownames(Z_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
Z_uniquely_up_noiseR <- Z_uniquely_up[setdiff(rownames(Z_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(Z_uniquely_up_noiseR) %in% "CD4"),
       Z_uniquely_up_noiseR[which(rownames(Z_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(Z_uniquely_up_noiseR) %in% "KLRB1"),
       Z_uniquely_up_noiseR[which(rownames(Z_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- Z_uniquely_up_noiseR[1:25, ]

# heatmap
plotHeatmap(
  cp,
  features = rownames(best_set),
  columns = order(
    cp$vs13,
    cp$cluster,
    cp$stage,
    cp$tissue,
    cp$donor,
    cp$group,
    cp$plate_number),
  colour_columns_by = c(
    "vs13",
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
    vs13 = vs13_colours,
    cluster = cluster_colours,
    stage = stage_colours,
    tissue = tissue_colours,
    donor = donor_colours,
    group = group_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)










# cluster 1 (i.e. S3-mix more thymus, peripheral IL7R driven
# cluster 2 (i.e. S3-mix more thymus, peripheral may CCL5 driven
# cluster 3 (i.e. S3-mix more thymus, peripheral may SOX4 drive
# cluster 4 (i.e. S3-mix more thymus, center

# as there is no standout feature in terms of experimental group for each cluster,
# all of them are S3-mix (with mostly thymus.s3 with some blood.s3)
# though, as in the minibulk with 20 DE only and we cannot tell difference between thymus.s3 and blood.s3
# it may indicate in both thymus and blood, they got maybe 4 sub-stages of S3 that are in common
#
# thus I start comparing all combination of single cluster to narrow down

# fx: per cluster compare

# 1 vs 2 (A vs B)
# cluster 1 (S3-mix more thymus, peripheral IL7R driven >>> IL7R and LTB as markers; DE not associated with tissue
# cluster 2 (S3-mix more thymus, peripheral may CCL5 driven >>> besides NPIPB family, number of markers found (e.g. CCL5 for both tissues; KLRD1 and EFHD2 mostly for blood)
# COMMENT: with ref to cluster 1, these are 2 separated subtype thymus.s3 (esp cluster 1: IL7R driven)

# 1 vs 3 (C vs D)
# cluster 1 (S3-mix more thymus, peripheral IL7R driven >>> besides NPIPB family, IL7R, IFITM3, etc shown as clear markers
# cluster 3 (S3-mix more thymus, peripheral may SOX4 drive >>> CCL5 as only unique makers (may look more highly expressed in "blood" cells)
# COMMENT: with ref to cluster 1, these are 2 separated subtype thymus.s3 (esp cluster 1 is IL7R IFITM3 driven, cluster 3 is CCL5 driven

# 1 vs 4 (E vs F)
# cluster 1 (S3-mix more thymus, peripheral IL7R driven >>> cluster 1 does show a number of cells more frequently expressed with lots of markers (not assocaited with tissue)
# cluster 4 (S3-mix more thymus, center >>> no sig marker found, though visually, there maybe some
# COMMENT: compared to cluster4, cluster 1 does have a number of genes more frequently express in 1 only, but not in 4; needed to keep them apart

# 2 vs 3 (G vs H)
# cluster 2 (S3-mix more thymus, peripheral may CCL5 driven >>> besides NPIPB, SYNE2 and AC009022.1 are clear marker
# cluster 3 (S3-mix more thymus, peripheral may SOX4 drive >>> statistically no, visually yes *may have some more frequently expressed genes
# COMMENT: SYNE2 and AC009022.1 (and maybe NPIPB) may driven cluster 2 from 3

# 2 vs 4 (I vs J)
# cluster 2 (S3-mix more thymus, peripheral may CCL5 driven >>> besides NPIPB, got number of markers (esp ACTG1)
# cluster 4 (S3-mix more thymus, center >>> statistically no, visually yes
# COMMENT: ACTG1 etc gene drove cluster 2 away from 4

# 3 vs 4 (K vs L)
# cluster 3 (S3-mix more thymus, peripheral may SOX4 drive >>> got lots of genes statistically sig express in cluster 3 only, e.g. ACTG1
# cluster 4 (S3-mix more thymus, center >>> statistically no, visually yes
# COMMENT: cluster 3 got number of unique markers (e.g. ACTG1) more frequently express and drove cluster 3 from 4






# fx: combine peripheral cluster compare with each peripheral cluster
# I: not necessary !?

# 1_2 vs 3 (M vs N)
# cluster 1_2 (S3-mix more thymus, peripheral, IL7R driven + may CCL5 driven >>> statistically no; also no visual marker with expression level common to both cluster 1 2, and 3
# cluster 3 (S3-mix more thymus, peripheral may SOX4 drive >>> statistically no, visually yes
# COMMENT: cluster 1 and 2 are too different; when compare something that is internally different with another cluster will not have enough statistical power ?

# 2_3 vs 1 (O vs P)
# cluster 2_3 (S3-mix more thymus, peripheral, may CCL5 and SOX4 driven >>> got several marker more frequently express in cluster 2 and 3 (e.g CCL5, EFHD2, CST7)
# cluster 1 (S3-mix more thymus, peripheral IL7R driven >>> got several marker as well (e.g. IL7R, LRRC75A; for LTB, IFITM3 and IFITM1, they seems to have higher expression in donor 3)
# COMMENT: cluster 1 clearly separated from 2_3, where 2 and 3 could be more similar; also explain why cluster 1 partner with 2 show no difference when compare to 3

# 1_3 vs 2 (Q vs R)
# cluster 1_3 (S3-mix more thymus, peripheral, IL7R driven + may SOX4 driven >>> lots of markers more frequent express in 1_3
# cluster 2 (S3-mix more thymus, peripheral may CCL5 driven >>> besides NPIPB, got CCL5, SYNE2 and AC009022.1 as markers
# COMMENT: if cluster 3 is in the combination, there will be DE (?); this means ... cluster 3 is more different from cluster 1 or 2 ?





# fx: combine peripheral cluster compare with core/center cluster

# 1_2 vs 4 (S vs T)
# cluster 1_2 (S3-mix more thymus, peripheral, IL7R driven + may CCL5 driven >>> lot of markers more frequently express in 1_2
# cluster 4 (i.e. S3-mix more thymus, center >>> statistically no, visually yes
# COMMENT: consistently, 4 show no marker, uniqueness is from cluster 1_2

# 2_3 vs 4 (U vs V)
# cluster 2_3 (S3-mix more thymus, peripheral, may CCL5 and SOX4 driven >>> lots of marker more frequently expressed in 2_3
# cluster 4 (i.e. S3-mix more thymus, center >>> statistically no, visually yes
# COMMENT: consistently, 4 show no marker, uniqueness is from cluster 2_3

# 1_3 vs 4 (W vs X)
# cluster 1_3 (S3-mix more thymus, peripheral, IL7R driven + may SOX4 driven >>> lots of marker more frequently expressed in 1_3
# cluster 4 (i.e. S3-mix more thymus, center >>> statistically no, visually yes
# COMMENT: consistently, 4 show no marker, uniqueness is from cluster 1_3

# 1_2_3 vs 4 (Y vs Z)
# cluster 1_2_3 (S3-mix more thymus, peripheral, IL7R driven + may SOX4, CCL5 driven >>> lots of marker more frequently expressed in 1_2_3
# cluster 4 (i.e. S3-mix more thymus, center >>> statistically no, visually yes
# COMMENT: consistently, 4 show no marker, uniqueness is from cluster 1_2_3


# CONCLUSION:
# as there is no standout feature in terms of experimental group for each cluster,
# all of them are S3-mix (with mostly thymus.s3 with some blood.s3)
# though, as in the minibulk with 20 DE only and we cannot tell difference between thymus.s3 and blood.s3
# it may indicate in both thymus and blood, they got maybe 4 sub-stages of S3 that are in common
#
# based on pairwise DE detection, the peripheral cluster (i.e. cluster 1, 2, 3) do have unique makers up-regulated and drive them apart
# e.g. cluster 1 (IL7R and LTB), cluster 2 (CCL5 for both tissues; KLRD1 and EFHD2 mostly for blood), cluster 3 (CCL5 for blood),
# whilst for the centre cluster, i.e. cluster 4
# it does not have any marker when compared to to any other peripheral cluster (1,2, or 3)
# it means in could be a cluster have all feature common to cluster 1,2 and 3;  could be a multipotent stem cells/ S3 cells lineage parent in stage S3 ??
#
# and this hypothesis is being double confirmed when compare either 1_2, 2_3, 1_3 and 1_2_3 with 4

















