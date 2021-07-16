
# NOTE: not suggest to narrow down into protein coding genes (pcg) as it remove all significant candidate in most of the comparison !!!
# TODO: is there any way to change the Sig colour using but still using plotHeatmap
# TODO: cluster 1, 2, 4 on their own don't show any unique marker at all !
# TODO: cluster 1,2 as a group, then compare cluster 1+2 (i.e. mostly Blood.S3) vs 3 (i.e. mostly Thymus.S1+.S2) vs 4 (i.e. mostly Thymus.S3)

library(here)
library(SingleCellExperiment)

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.cells_selected.SCE.rds"))

# dir.create(here("data", "marker_genes"), recursive = TRUE)
# dir.create(here("output", "marker_genes"), recursive = TRUE)

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
sample_colours <- setNames(
  unique(sce$colours$sample_colours),
  unique(names(sce$colours$sample_colours)))
sample_colours <- sample_colours[levels(sce$sample)]
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

# summary - UMAP
p1 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "cluster", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p2 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "sample", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_colours, name = "sample")
p3 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "stage", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = stage_colours, name = "stage")
p4 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p5 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p6 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p7 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "group", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = group_colours, name = "group")
(p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | plot_spacer())

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
    aes(x = cluster, fill = sample),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = sample_colours) +
  theme_cowplot(font_size = 8)
p3 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = stage),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = stage_colours) +
  theme_cowplot(font_size = 8)
p4 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = plate_number),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = plate_number_colours) +
  theme_cowplot(font_size = 8)
p5 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = tissue),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = tissue_colours) +
  theme_cowplot(font_size = 8)
p6 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = donor),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = donor_colours) +
  theme_cowplot(font_size = 8)
p7 <- ggcells(sce) +
  geom_bar(
    aes(x = cluster, fill = group),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = group_colours) +
  theme_cowplot(font_size = 8)
(p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | plot_spacer())

# block on plate
sce$block <- paste0(sce$plate_number)
















###################################
# cluster 1 (i.e. S3 more Blood) vs 2 (i.e. S3 more Thymus) vs 3 (i.e. mostly Thymus.S1 + .S2) vs 4 (i.e. mostly Thymus.S3)

# find unique DE ./. clusters
uniquely_up <- findMarkers(
  sce,
  groups = sce$cluster,
  block = sce$block,
  pval.type = "all",
  direction = "up")

##########################################
# look at cluster 1 (i.e. S3 more Blood)
chosen <- "1"
cluster1_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3 more Blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster1_uniquely_up <- cluster1_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster1_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
cluster1_uniquely_up_noiseR <- cluster1_uniquely_up[setdiff(rownames(cluster1_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster1_uniquely_up_noiseR) %in% "CD4"),
       cluster_uniquely_up_noiseR[which(rownames(cluster1_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster1_uniquely_up_noiseR) %in% "KLRB1"),
       cluster_uniquely_up_noiseR[which(rownames(cluster1_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- cluster1_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "group",
    "sample",
    "stage",
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
  main = paste0("Cluster: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



##########################################
# look at cluster 2 (i.e. S3 more Thymus)
chosen <- "2"
cluster2_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3 more Thymus)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster2_uniquely_up <- cluster2_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster2_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
cluster2_uniquely_up_noiseR <- cluster2_uniquely_up[setdiff(rownames(cluster2_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster2_uniquely_up_noiseR) %in% "CD4"),
       cluster2_uniquely_up_noiseR[which(rownames(cluster2_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster2_uniquely_up_noiseR) %in% "KLRB1"),
       cluster2_uniquely_up_noiseR[which(rownames(cluster2_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- cluster2_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "group",
    "sample",
    "stage",
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
      levels = c("No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



##########################################
# look at cluster 3 (i.e. mostly Thymus.S1 + .S2)
chosen <- "3"
cluster3_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Thymus.S1 + .S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster3_uniquely_up <- cluster3_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster3_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
cluster3_uniquely_up_noiseR <- cluster3_uniquely_up[setdiff(rownames(cluster3_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster3_uniquely_up_noiseR) %in% "CD4"),
       cluster3_uniquely_up_noiseR[which(rownames(cluster3_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster3_uniquely_up_noiseR) %in% "KLRB1"),
       cluster3_uniquely_up_noiseR[which(rownames(cluster3_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- cluster3_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "group",
    "sample",
    "stage",
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
  main = paste0("Cluster: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



##########################################
# look at cluster 4 (i.e. mostly Thymus.S3)
chosen <- "4"
cluster4_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster4_uniquely_up <- cluster4_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster4_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
cluster4_uniquely_up_noiseR <- cluster4_uniquely_up[setdiff(rownames(cluster4_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster4_uniquely_up_noiseR) %in% "CD4"),
       cluster4_uniquely_up_noiseR[which(rownames(cluster4_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster4_uniquely_up_noiseR) %in% "KLRB1"),
       cluster4_uniquely_up_noiseR[which(rownames(cluster4_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- cluster4_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
    "group",
    "sample",
    "stage",
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
  main = paste0("Cluster: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)












##########################################################################
# cluster 1 + 2 + 4 (i.e. mostly S3) vs cluster 3 (i.e. mostly S1 and S2)

# classify cluster-group for comparison
sce$vs1 <- factor(ifelse(sce$cluster == 1 | sce$cluster == 2 | sce$cluster == 4, "A", "B"))

# set vs colours
vs1_colours <- setNames(
  palette.colors(nlevels(sce$vs1), "Set1"),
  levels(sce$vs1))
sce$colours$vs1_colours <- vs1_colours[sce$vs1]

# find unique DE ./. cluster-groups
vs1_uniquely_up <- findMarkers(
  sce,
  groups = sce$vs1,
  block = sce$block,
  pval.type = "all",
  direction = "up")


##########################################
# look at cluster-group A (i.e. mostly S3)
chosen <- "A"
A_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# A_uniquely_up_pcg <- A_uniquely_up[intersect(protein_coding_gene_set, rownames(A_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
A_uniquely_up_noiseR <- A_uniquely_up[setdiff(rownames(A_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
  which(rownames(A_uniquely_up_noiseR) %in% "CD4"),
  A_uniquely_up_noiseR[which(rownames(A_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
  which(rownames(A_uniquely_up_noiseR) %in% "KLRB1"),
  A_uniquely_up_noiseR[which(rownames(A_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- A_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$vs1,
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "vs1",
    "cluster",
    "group",
    "sample",
    "stage",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      levels = c("Yes", "No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs1 = vs1_colours,
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)


###############################################
# look at cluster-group B (i.e. mostly S1 + S2)
chosen <- "B"
B_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly S1 and S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# B_uniquely_up_pcg <- B_uniquely_up[intersect(protein_coding_gene_set, rownames(B_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
B_uniquely_up_noiseR <- B_uniquely_up[setdiff(rownames(B_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(B_uniquely_up_noiseR) %in% "CD4"),
       B_uniquely_up_noiseR[which(rownames(B_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(B_uniquely_up_noiseR) %in% "KLRB1"),
       B_uniquely_up_noiseR[which(rownames(B_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- B_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$vs1,
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "vs1",
    "cluster",
    "group",
    "sample",
    "stage",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      levels = c("Yes", "No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs1 = vs1_colours,
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



##########################################################################
# cluster 1 + 2 (i.e. mostly Blood.S3) vs cluster 3 (i.e. mostly Thymus.S1 and Thymus.S2)

# classify cluster-group for comparison
# NOTE: NI = not interested
sce$vs2 <- factor(ifelse(sce$cluster == 1 | sce$cluster == 2 , "C",
                         ifelse(sce$cluster == 3 , "D", "NI")))

# set vs colours
vs2_colours <- setNames(
  palette.colors(nlevels(sce$vs2), "Set1"),
  levels(sce$vs2))
sce$colours$vs2_colours <- vs2_colours[sce$vs2]

# find unique DE ./. cluster-groups
vs2_uniquely_up <- findMarkers(
  sce,
  groups = sce$vs2,
  block = sce$block,
  pval.type = "all",
  direction = "up")


#################################################
# look at cluster-group C  (i.e. mostly Blood.S3)
chosen <- "C"
C_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Blood.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# C_uniquely_up_pcg <- C_uniquely_up[intersect(protein_coding_gene_set, rownames(C_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
C_uniquely_up_noiseR <- C_uniquely_up[setdiff(rownames(C_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(C_uniquely_up_noiseR) %in% "CD4"),
       C_uniquely_up_noiseR[which(rownames(C_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(C_uniquely_up_noiseR) %in% "KLRB1"),
       C_uniquely_up_noiseR[which(rownames(C_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- C_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$vs2,
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "vs2",
    "cluster",
    "group",
    "sample",
    "stage",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      levels = c("Yes", "No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen,  " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs2 = vs2_colours,
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)


##################################################
# look at cluster-group D  (i.e. mostly Thymus.S1 and Thymus.S2)
chosen <- "D"
D_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Thymus.S1 and Thymus.S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# D_uniquely_up_pcg <- D_uniquely_up[intersect(protein_coding_gene_set, rownames(D_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
D_uniquely_up_noiseR <- D_uniquely_up[setdiff(rownames(D_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(D_uniquely_up_noiseR) %in% "CD4"),
       D_uniquely_up_noiseR[which(rownames(D_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(D_uniquely_up_noiseR) %in% "KLRB1"),
       D_uniquely_up_noiseR[which(rownames(D_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- D_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$vs2,
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "vs2",
    "cluster",
    "group",
    "sample",
    "stage",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      levels = c("Yes", "No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    vs2 = vs2_colours,
    cluster = cluster_colours,
    group = group_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



##########################################################################
# cluster 1 + 2 (i.e. mostly Blood.S3) vs cluster 4 (i.e. mostly Thymus.S3)

# classify cluster-group for comparison
# NOTE: NI = not interested
sce$vs3 <- factor(ifelse(sce$cluster == 1 | sce$cluster == 2 , "E",
                         ifelse(sce$cluster == 4 , "F", "NI")))

# set vs colours
vs3_colours <- setNames(
  palette.colors(nlevels(sce$vs3), "Set1"),
  levels(sce$vs3))
sce$colours$vs3_colours <- vs3_colours[sce$vs3]

# find unique DE ./. cluster-groups
vs3_uniquely_up <- findMarkers(
  sce,
  groups = sce$vs3,
  block = sce$block,
  pval.type = "all",
  direction = "up")


#################################################
# look at cluster-group E  (i.e. mostly Blood.S3)
chosen <- "E"
E_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Blood.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# E_uniquely_up_pcg <- E_uniquely_up[intersect(protein_coding_gene_set, rownames(E_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
E_uniquely_up_noiseR <- E_uniquely_up[setdiff(rownames(E_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(E_uniquely_up_noiseR) %in% "CD4"),
       E_uniquely_up_noiseR[which(rownames(E_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(E_uniquely_up_noiseR) %in% "KLRB1"),
       E_uniquely_up_noiseR[which(rownames(E_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- E_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$vs3,
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "vs3",
    "cluster",
    "group",
    "sample",
    "stage",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      levels = c("Yes", "No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen,  " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),  column_annotation_colors = list(
                  Sig = c("Yes" = "red", "No" = "lightgrey"),
                  vs3 = vs3_colours,
                  cluster = cluster_colours,
                  group = group_colours,
                  sample = sample_colours,
                  stage = stage_colours,
                  plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)


##################################################
# look at cluster-group F  (i.e. mostly Thymus.S3)
chosen <- "F"
F_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# F_uniquely_up_pcg <- F_uniquely_up[intersect(protein_coding_gene_set, rownames(F_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito) that collaborator not interested in
F_uniquely_up_noiseR <- F_uniquely_up[setdiff(rownames(F_uniquely_up), c(pseudogene_set, mito_set, ribo_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(F_uniquely_up_noiseR) %in% "CD4"),
       F_uniquely_up_noiseR[which(rownames(F_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(F_uniquely_up_noiseR) %in% "KLRB1"),
       F_uniquely_up_noiseR[which(rownames(F_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top50 only
best_set <- F_uniquely_up_noiseR[1:50, ]

# heatmap
plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$vs3,
    sce$cluster,
    sce$group,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "vs3",
    "cluster",
    "group",
    "sample",
    "stage",
    "plate_number"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  annotation_row = data.frame(
    Sig = factor(
      ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
      levels = c("Yes", "No")),
    row.names = rownames(best_set)),
  main = paste0("Cluster-group: ", chosen, " ", x, " - ",
                y[1], "_top ", y[2], "_significance: ", y[3], " ; ",
                z[1], "_top ", z[2], "_significance: ", z[3]),  column_annotation_colors = list(
                  Sig = c("Yes" = "red", "No" = "lightgrey"),
                  vs3 = vs3_colours,
                  cluster = cluster_colours,
                  group = group_colours,
                  sample = sample_colours,
                  stage = stage_colours,
                  plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)










