
# NOTE: not suggest to narrow down into protein coding genes (pcg) as it remove all significant candidate in most of the comparison !!!
# TODO: is there any way to change the Sig colour using but still using plotHeatmap
# TODO: cluster 1,2 as a group, then compare cluster 1+2 (i.e. mostly Blood.S3) vs 3 (i.e. mostly Thymus.S1+.S2) vs 4 (i.e. mostly Thymus.S3)

library(here)
library(SingleCellExperiment)

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.cells_selected.SCE.rds"))

# dir.create(here("data", "marker_genes"), recursive = TRUE)
# dir.create(here("output", "marker_genes"), recursive = TRUE)

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

# gene-of-interest (suggested by Dan at meeting on 15 Jul 2021)
interest <- c("KLRB1", "CD4", # QC
             "NKG7", "GNLY", "CCL5", # cluster 1 selected unique markers
             "TCF7", "SOX4", "LEF1", "CD3D", "CCR9", "BCL11B", "NREP" # cluster 3 selected unique markers
             )

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



###################################
# cluster 1 (i.e. S3-mix more Blood) vs 2 (i.e. S3-mix more Thymus) vs 3 (i.e. mostly Thymus.S1 + .S2) vs 4 (i.e. mostly Thymus.S3)

# find unique DE ./. clusters
uniquely_up <- findMarkers(
  sce,
  groups = sce$cluster,
  block = sce$block,
  pval.type = "all",
  direction = "up")

##########################################
# look at cluster 1 (i.e. S3-mix more Blood)
chosen <- "1"
cluster1_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3-mix more Blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# cluster1_uniquely_up <- cluster1_uniquely_up[intersect(protein_coding_gene_set, rownames(cluster1_uniquely_up)), ]

# get rid of noise (i.e. pseudo, ribo, mito, sex) that collaborator not interested in
cluster1_uniquely_up_noiseR <- cluster1_uniquely_up[setdiff(rownames(cluster1_uniquely_up), c(pseudogene_set, mito_set, ribo_set, sex_set)), ]

# see if key marker, "CD4 and/or ""KLRB1/CD161"", contain in the DE list + if it is "significant (i.e FDR <0.05)
y <- c("CD4",
       which(rownames(cluster1_uniquely_up_noiseR) %in% "CD4"),
       cluster_uniquely_up_noiseR[which(rownames(cluster1_uniquely_up_noiseR) %in% "CD4"), ]$FDR < 0.05)
z <- c("KLRB1/CD161",
       which(rownames(cluster1_uniquely_up_noiseR) %in% "KLRB1"),
       cluster_uniquely_up_noiseR[which(rownames(cluster1_uniquely_up_noiseR) %in% "KLRB1"), ]$FDR < 0.05)

# top25 only + gene-of-interest
best_set <- cluster1_uniquely_up_noiseR[union(rownames(cluster1_uniquely_up_noiseR)[1:25], interest),]

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
# look at cluster 2 (i.e. S3-mix more Thymus)
chosen <- "2"
cluster2_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(S3-mix more Thymus)"

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

# top25 only + gene-of-interest
best_set <- cluster2_uniquely_up_noiseR[union(rownames(cluster2_uniquely_up_noiseR)[1:25], interest),]

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
  # annotation_row = data.frame(
  #   Sig = factor(
  #     ifelse(best_set[, "FDR"] < 0.05, "Yes", "No"),
  #     # TODO: temp trick to deal with the row-colouring problem
  #     # levels = c("Yes", "No")),
  #     levels = c("No")),
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
# look at cluster 3 (i.e. mostly Thymus.S1 + .S2)
chosen <- "3"
cluster3_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Thymus.S1 + .S2)"

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

# top25 only + gene-of-interest
best_set <- cluster3_uniquely_up_noiseR[union(rownames(cluster3_uniquely_up_noiseR)[1:25], interest),]

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
# look at cluster 4 (i.e. mostly Thymus.S3)
chosen <- "4"
cluster4_uniquely_up <- uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(mostly Thymus.S3)"

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

# top25 only + gene-of-interest
best_set <- cluster4_uniquely_up_noiseR[union(rownames(cluster4_uniquely_up_noiseR)[1:25], interest),]

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








#########
# A vs B
#########

##########################################################################
# cluster 1 + 2 + 4 (i.e. mostly S3) vs cluster 3 (i.e. mostly S1 and S2)

# checkpoint
cp <- sce

# classify cluster-group for comparison
cp$vs1 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 2 | cp$cluster == 4, "A", "B"))

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


##########################################################
# look at cluster-group A / cluster 1+2+4 (i.e. mostly S3)
chosen <- "A"
A_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1+2+4; mostly S3)"

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
best_set <- A_uniquely_up_noiseR[union(rownames(A_uniquely_up_noiseR)[1:25], interest),]

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

###########################################################
# look at cluster-group B / cluster 3 (i.e. mostly S1 + S2)
chosen <- "B"
B_uniquely_up <- vs1_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly S1 and S2)"

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
best_set <- B_uniquely_up_noiseR[union(rownames(B_uniquely_up_noiseR)[1:25], interest),]

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

##########################################################################################
# cluster 1 + 2 (i.e. mostly Blood.S3) vs cluster 3 (i.e. mostly Thymus.S1 and Thymus.S2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster != "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs2 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 2 , "C", "D"))

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

###############################################################
# look at cluster-group C / cluster 1+2 (i.e. mostly Blood.S3)
chosen <- "C"
C_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1+2; mostly Blood.S3)"

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
best_set <- C_uniquely_up_noiseR[union(rownames(C_uniquely_up_noiseR)[1:25], interest),]

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

###########################################################################
# look at cluster-group D / cluster 3 (i.e. mostly Thymus.S1 and Thymus.S2)
chosen <- "D"
D_uniquely_up <- vs2_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly Thymus.S1 and Thymus.S2)"

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
best_set <- D_uniquely_up_noiseR[union(rownames(D_uniquely_up_noiseR)[1:25], interest),]

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

##########################################################################
# cluster 1 + 2 (i.e. mostly Blood.S3) vs cluster 4 (i.e. mostly Thymus.S3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster != "3"]
colData(cp) <- droplevels(colData(cp))
# classify cluster-group for comparison

cp$vs3 <- factor(ifelse(cp$cluster == 1 | cp$cluster == 2 , "E", "F"))

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

###############################################################
# look at cluster-group E / cluster 1+2 (i.e. mostly Blood.S3)
chosen <- "E"
E_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1+2; mostly Blood.S3)"

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
best_set <- E_uniquely_up_noiseR[union(rownames(E_uniquely_up_noiseR)[1:25], interest),]

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

##############################################################
# look at cluster-group F / cluster 4 (i.e. mostly Thymus.S3)
chosen <- "F"
F_uniquely_up <- vs3_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly Thymus.S3)"

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
best_set <- F_uniquely_up_noiseR[union(rownames(F_uniquely_up_noiseR)[1:25], interest),]

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

##########################################################################
# cluster 1 (i.e. S3-mix more Blood) vs cluster 2 (i.e. S3-mix more Thymus)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "2"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs4 <- factor(ifelse(cp$cluster == 1 , "G", "H"))

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

###########################################################
# look at cluster-group G / cluster 1 (i.e. S3-mix more Blood)
chosen <- "G"
G_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix more Blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# G_uniquely_up_pcg <- G_uniquely_up[intersect(protein_coding_genG_set, rownames(G_uniquely_up)), ]

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
best_set <- G_uniquely_up_noiseR[union(rownames(G_uniquely_up_noiseR)[1:25], interest),]

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

############################################################
# look at cluster-group H / cluster 2 (i.e. S3-mix more Thymus)
chosen <- "H"
H_uniquely_up <- vs4_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# H_uniquely_up_pcg <- H_uniquely_up[intersect(protein_coding_genG_set, rownames(H_uniquely_up)), ]

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
best_set <- H_uniquely_up_noiseR[union(rownames(H_uniquely_up_noiseR)[1:25], interest),]

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
  #     levels = c("No")),
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

#################################################################################
# cluster 3 (i.e. mostly Thymus.S1 and .S2) vs cluster 4 (i.e. mostly Thymus.S3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "3" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs5 <- factor(ifelse(cp$cluster == 3 , "I", "J"))

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

#####################################################################
# look at cluster-group I / cluster 3 (i.e. mostly Thymus.S1 and .S2)
chosen <- "I"
I_uniquely_up <- vs5_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly Thymus.S1 and .S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# I_uniquely_up_pcg <- I_uniquely_up[intersect(protein_coding_genI_set, rownames(I_uniquely_up)), ]

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
best_set <- I_uniquely_up_noiseR[union(rownames(I_uniquely_up_noiseR)[1:25], interest),]

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

##############################################################
# look at cluster-group J / cluster 4  (i.e. mostly Thymus.S3)
chosen <- "J"
J_uniquely_up <- vs5_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# J_uniquely_up_pcg <- J_uniquely_up[intersect(protein_coding_genI_set, rownames(J_uniquely_up)), ]

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
best_set <- J_uniquely_up_noiseR[union(rownames(J_uniquely_up_noiseR)[1:25], interest),]

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



#########
# K vs L
#########

#################################################################################
# cluster 2 (i.e. S3-mix more Thymus) vs cluster 3 (i.e. mostly Thymus.S1 + .S2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs6 <- factor(ifelse(cp$cluster == 2 , "K", "L"))

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

#####################################################################
# look at cluster-group K / cluster 2 (i.e. S3-mix more Thymus)
chosen <- "K"
K_uniquely_up <- vs6_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix more Thymus)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# K_uniquely_up_pcg <- K_uniquely_up[intersect(protein_coding_genK_set, rownames(K_uniquely_up)), ]

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
best_set <- K_uniquely_up_noiseR[union(rownames(K_uniquely_up_noiseR)[1:25], interest),]

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

##############################################################
# look at cluster-group L / cluster 3  (i.e. mostly Thymus.S1 + .S2)
chosen <- "L"
L_uniquely_up <- vs6_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# L_uniquely_up_pcg <- L_uniquely_up[intersect(protein_coding_genK_set, rownames(L_uniquely_up)), ]

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
best_set <- L_uniquely_up_noiseR[union(rownames(L_uniquely_up_noiseR)[1:25], interest),]

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
# M vs N
#########

#################################################################################
# cluster 2 (i.e. S3-mix more Thymus) vs cluster 4 (i.e. mostly Thymus.S3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "2" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs7 <- factor(ifelse(cp$cluster == 2 , "M", "N"))

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

#####################################################################
# look at cluster-group M / cluster 2 (i.e. S3-mix more Thymus)
chosen <- "M"
M_uniquely_up <- vs7_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 2; S3-mix more Thymus)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# M_uniquely_up_pcg <- M_uniquely_up[intersect(protein_coding_genM_set, rownames(M_uniquely_up)), ]

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
best_set <- M_uniquely_up_noiseR[union(rownames(M_uniquely_up_noiseR)[1:25], interest),]

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
# look at cluster-group N / cluster 4  (i.e. mostly Thymus.S3)
chosen <- "N"
N_uniquely_up <- vs7_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# N_uniquely_up_pcg <- N_uniquely_up[intersect(protein_coding_genM_set, rownames(N_uniquely_up)), ]

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
best_set <- N_uniquely_up_noiseR[union(rownames(N_uniquely_up_noiseR)[1:25], interest),]

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



#########
# O vs P
#########

#################################################################################
# cluster 1 (i.e. S3-mix more Blood) vs cluster 3 (i.e. mostly Thymus.S1 and .S2)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "3"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs8 <- factor(ifelse(cp$cluster == 1 , "O", "P"))

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

#####################################################################
# look at cluster-group O / cluster 1 (i.e. S3-mix more Blood)
chosen <- "O"
O_uniquely_up <- vs8_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix more Blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# O_uniquely_up_pcg <- O_uniquely_up[intersect(protein_coding_genO_set, rownames(O_uniquely_up)), ]

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
best_set <- O_uniquely_up_noiseR[union(rownames(O_uniquely_up_noiseR)[1:25], interest),]

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
# look at cluster-group P / cluster 3  (i.e. mostly Thymus.S1 and .S2)
chosen <- "P"
P_uniquely_up <- vs8_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 3; mostly Thymus.S1 and .S2)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# P_uniquely_up_pcg <- P_uniquely_up[intersect(protein_coding_genO_set, rownames(P_uniquely_up)), ]

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
best_set <- P_uniquely_up_noiseR[union(rownames(P_uniquely_up_noiseR)[1:25], interest),]

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

#################################################################################
# cluster 1 (i.e. S3-mix more Blood) vs cluster 4 (i.e. mostly Thymus.S3)

# checkpoint
cp <- sce
# exclude cells of uninterested cluster from cp
cp <- cp[, cp$cluster == "1" | cp$cluster == "4"]
colData(cp) <- droplevels(colData(cp))

# classify cluster-group for comparison
cp$vs9 <- factor(ifelse(cp$cluster == 1 , "Q", "R"))

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

#####################################################################
# look at cluster-group Q / cluster 1 (i.e. S3-mix more Blood)
chosen <- "Q"
Q_uniquely_up <- vs9_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 1; S3-mix more Blood)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# Q_uniquely_up_pcg <- Q_uniquely_up[intersect(protein_coding_genQ_set, rownames(Q_uniquely_up)), ]

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
best_set <- Q_uniquely_up_noiseR[union(rownames(Q_uniquely_up_noiseR)[1:25], interest),]

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
# look at cluster-group R / cluster 4  (i.e. mostly Thymus.S3)
chosen <- "R"
R_uniquely_up <- vs9_uniquely_up[[chosen]]

# add description for the chosen cluster-group
x <- "(cluster 4; mostly Thymus.S3)"

# look only at protein coding gene (pcg)
# NOTE: not suggest to narrow down into pcg as it remove all significant candidates (FDR << 0.05) !
# R_uniquely_up_pcg <- R_uniquely_up[intersect(protein_coding_genQ_set, rownames(R_uniquely_up)), ]

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
best_set <- R_uniquely_up_noiseR[union(rownames(R_uniquely_up_noiseR)[1:25], interest),]

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
















############################################
# Sanity check (TRDV2 and TRGV9 expression)

assay(sce, "log2cpm") <- edgeR::cpm(
  counts(sce),
  log = TRUE,
  lib.size = colSums(counts(sce)))

plot_grid(

  plotExpression(
    sce,
    features = "TRDV2",
    x = "cluster",
    exprs_values = "log2cpm",
    colour_by = "cluster",
    other_fields = "stage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_manual(values = cluster_colours, name = "cluster") +
    facet_grid(~stage),
  plotExpression(
    sce,
    features = "TRGV9",
    x = "cluster",
    exprs_values = "log2cpm",
    colour_by = "cluster",
    other_fields = "stage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_manual(values = cluster_colours, name = "cluster") +
    facet_grid(~stage),
  ncol= 1,
  align ="h"
)

# NOTE: as suggested by Pete, we have no idea how well the CELSeq protocol could capture the expression of the T-cell-receptor component
# The suggested known markers , i.e. TRDV2, TRGV9, may not be a good candidate for defining whether a cluster contains cells relevant to gamma-delta T cells














############
# cluster 3
############

chosen <- "3"

# add description for the chosen cluster-group
x <- "(mostly Thymus.S1 and Thymus.S2)"

# retain only significant markers (FDR<0.05) + keep only required output columns
B_uniquely_up_noiseR_sig <- B_uniquely_up_noiseR[B_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
D_uniquely_up_noiseR_sig <- D_uniquely_up_noiseR[D_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
I_uniquely_up_noiseR_sig <- I_uniquely_up_noiseR[I_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
L_uniquely_up_noiseR_sig <- L_uniquely_up_noiseR[L_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
P_uniquely_up_noiseR_sig <- P_uniquely_up_noiseR[P_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
cluster3_uniquely_up_noiseR_sig <- cluster3_uniquely_up_noiseR[cluster3_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]

# add top column
B_uniquely_up_noiseR_sig$top <- 1:nrow(B_uniquely_up_noiseR_sig)
D_uniquely_up_noiseR_sig$top <- 1:nrow(D_uniquely_up_noiseR_sig)
I_uniquely_up_noiseR_sig$top <- 1:nrow(I_uniquely_up_noiseR_sig)
L_uniquely_up_noiseR_sig$top <- 1:nrow(L_uniquely_up_noiseR_sig)
P_uniquely_up_noiseR_sig$top <- 1:nrow(P_uniquely_up_noiseR_sig)
cluster3_uniquely_up_noiseR_sig$top <- 1:nrow(cluster3_uniquely_up_noiseR_sig)

# unify S4 objects, sort by top (ascending) then FDR (ascending), keep only first unique entry for each marker
BDILP3_uniquely_up_noiseR_sig <- rbind2(B_uniquely_up_noiseR_sig,
                                        D_uniquely_up_noiseR_sig,
                                        I_uniquely_up_noiseR_sig,
                                        L_uniquely_up_noiseR_sig,
                                        P_uniquely_up_noiseR_sig,
                                        cluster3_uniquely_up_noiseR_sig)

BDILP3_uniquely_up_noiseR_sig_sort <- BDILP3_uniquely_up_noiseR_sig[with(BDILP3_uniquely_up_noiseR_sig, order(top, FDR)), ]
BDILP3_uniquely_up_noiseR_sig_sort_uniq <- BDILP3_uniquely_up_noiseR_sig_sort[unique(rownames(BDILP3_uniquely_up_noiseR_sig_sort)), ]

# # de-select unannotated/ not well-characterised genes
# deselected <- c("NPIPB13", "NPIPB3", "NPIPB11", "NPIPB5", "NPIPB4", "EEF1A1", "ACTG1", "ACTB", "IFITM1")
deselected <- c()
BDILP3_uniquely_up_noiseR_sig_sort_uniq_selected <- BDILP3_uniquely_up_noiseR_sig_sort_uniq[!(rownames(BDILP3_uniquely_up_noiseR_sig_sort_uniq) %in% deselected), ]

# top only + gene-of-interest
best_set <- BDILP3_uniquely_up_noiseR_sig_sort_uniq_selected[1:50, ]

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
  annotation_row = data.frame(
    # TODO: temp trick to deal with the row-colouring problem
    cluster3.vs.1.2.4 = factor(ifelse(rownames(best_set) %in% rownames(B_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster3.vs.1.2 = factor(ifelse(rownames(best_set) %in% rownames(D_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster3.vs.4 = factor(ifelse(rownames(best_set) %in% rownames(I_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster3.vs.2 = factor(ifelse(rownames(best_set) %in% rownames(L_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster3.vs.1 = factor(ifelse(rownames(best_set) %in% rownames(P_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster3.vs.1.vs.2.vs.4 = factor(ifelse(rownames(best_set) %in% rownames(cluster3_uniquely_up_noiseR_sig_trim), "DE", "not DE"), levels = c("DE")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x),
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



############
# cluster 4
############

chosen <- "4"

# add description for the chosen cluster-group
x <- "(mostly Thymus.S3)"

# retain only significant markers (FDR<0.05) + keep only required output columns
F_uniquely_up_noiseR_sig <- F_uniquely_up_noiseR[F_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
J_uniquely_up_noiseR_sig <- J_uniquely_up_noiseR[J_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
N_uniquely_up_noiseR_sig <- N_uniquely_up_noiseR[N_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
R_uniquely_up_noiseR_sig <- R_uniquely_up_noiseR[R_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
cluster4_uniquely_up_noiseR_sig <- cluster4_uniquely_up_noiseR[cluster4_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]

# add top column
F_uniquely_up_noiseR_sig$top <- 1:nrow(F_uniquely_up_noiseR_sig)
J_uniquely_up_noiseR_sig$top <- 1:nrow(J_uniquely_up_noiseR_sig)
N_uniquely_up_noiseR_sig$top <- 1:nrow(N_uniquely_up_noiseR_sig)
R_uniquely_up_noiseR_sig$top <- 1:nrow(R_uniquely_up_noiseR_sig)
cluster4_uniquely_up_noiseR_sig$top <- 1:nrow(cluster4_uniquely_up_noiseR_sig)

# unify S4 objects, sort by top (ascending) then FDR (ascending), keep only first unique entry for each marker
FJNR4_uniquely_up_noiseR_sig <- rbind2(F_uniquely_up_noiseR_sig,
                                       J_uniquely_up_noiseR_sig,
                                       N_uniquely_up_noiseR_sig,
                                       R_uniquely_up_noiseR_sig,
                                       cluster4_uniquely_up_noiseR_sig)

FJNR4_uniquely_up_noiseR_sig_sort <- FJNR4_uniquely_up_noiseR_sig[with(FJNR4_uniquely_up_noiseR_sig, order(top, FDR)), ]
FJNR4_uniquely_up_noiseR_sig_sort_uniq <- FJNR4_uniquely_up_noiseR_sig_sort[unique(rownames(FJNR4_uniquely_up_noiseR_sig_sort)), ]

# # de-select unannotated/ not well-characterised genes
# deselected <- c("MTRNR1L12", "MTRNR2L8", "NPIPB12", "NPIPB4", "NPIPB3", "NPIPB13", "NPIPB11", "NPIPB15", "MALAT1")
deselected <- c()
FJNR_uniquely_up_noiseR_sig_sort_uniq_selected <- FJNR_uniquely_up_noiseR_sig_sort_uniq[!(rownames(FJNR_uniquely_up_noiseR_sig_sort_uniq) %in% deselected), ]

# top only + gene-of-interest
best_set <- FJNR_uniquely_up_noiseR_sig_sort_uniq_selected[1:18, ]

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
  annotation_row = data.frame(
    # TODO: temp trick to deal with the row-colouring problem
    cluster4.vs.1.2 = factor(ifelse(rownames(best_set) %in% rownames(F_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster4.vs.3 = factor(ifelse(rownames(best_set) %in% rownames(J_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster4.vs.2 = factor(ifelse(rownames(best_set) %in% rownames(N_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster4.vs.1 = factor(ifelse(rownames(best_set) %in% rownames(R_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster4.vs.1.vs.2.vs.3 = factor(ifelse(rownames(best_set) %in% rownames(cluster4_uniquely_up_noiseR_sig_trim), "DE", "not DE"), levels = c("DE")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x),
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



############
# cluster 2
############

chosen <- "2"

# add description for the chosen cluster-group
x <- "(S3-mix more Thymus)"

# retain only significant markers (FDR<0.05) + keep only required output columns
K_uniquely_up_noiseR_sig <- K_uniquely_up_noiseR[K_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
M_uniquely_up_noiseR_sig <- M_uniquely_up_noiseR[M_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
H_uniquely_up_noiseR_sig <- H_uniquely_up_noiseR[H_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
cluster2_uniquely_up_noiseR_sig <- cluster2_uniquely_up_noiseR[cluster2_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]

# add top column
K_uniquely_up_noiseR_sig$top <- 1:nrow(K_uniquely_up_noiseR_sig)
M_uniquely_up_noiseR_sig$top <- 1:nrow(M_uniquely_up_noiseR_sig)
# TODO: fix error - exclude list from rbind2 if empty
H_uniquely_up_noiseR_sig$top <- 1:nrow(H_uniquely_up_noiseR_sig)
cluster2_uniquely_up_noiseR_sig$top <- 1:nrow(cluster2_uniquely_up_noiseR_sig)

# unify S4 objects, sort by top (ascending) then FDR (ascending), keep only first unique entry for each marker
KMH2_uniquely_up_noiseR_sig <- rbind2(K_uniquely_up_noiseR_sig,
                                      M_uniquely_up_noiseR_sig,
                                      H_uniquely_up_noiseR_sig,
                                      cluster2_uniquely_up_noiseR_sig)

KMH2_uniquely_up_noiseR_sig_sort <- KMH2_uniquely_up_noiseR_sig[with(KMH2_uniquely_up_noiseR_sig, order(top, FDR)), ]
KMH2_uniquely_up_noiseR_sig_sort_uniq <- KMH2_uniquely_up_noiseR_sig_sort[unique(rownames(KMH2_uniquely_up_noiseR_sig_sort)), ]

# # de-select unannotated/ not well-characterised genes
deselected <- c()
KMH_uniquely_up_noiseR_sig_sort_uniq_selected <- KMH_uniquely_up_noiseR_sig_sort_uniq[!(rownames(KMH_uniquely_up_noiseR_sig_sort_uniq) %in% deselected), ]

# top only + gene-of-interest
best_set <- KMH_uniquely_up_noiseR_sig_sort_uniq_selected[1:50, ]

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
  annotation_row = data.frame(
    # TODO: temp trick to deal with the row-colouring problem
    cluster2.vs.3 = factor(ifelse(rownames(best_set) %in% rownames(K_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster2.vs.4 = factor(ifelse(rownames(best_set) %in% rownames(M_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    # TODO: work out how to remove `fill` from `not_DE`
    cluster2.vs.1 = factor(ifelse(rownames(best_set) %in% rownames(H_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE", "not DE")),
    cluster2.vs.1.vs.3.vs.4 = factor(ifelse(rownames(best_set) %in% rownames(cluster2_uniquely_up_noiseR_sig_trim), "DE", "not DE"), levels = c("DE", "not DE")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x),
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



############
# cluster 1
############

chosen <- "1"

# add description for the chosen cluster-group
x <- "(S3-mix more Blood)"

# retain only significant markers (FDR<0.05) + keep only required output columns
G_uniquely_up_noiseR_sig <- G_uniquely_up_noiseR[G_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
O_uniquely_up_noiseR_sig <- O_uniquely_up_noiseR[O_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
Q_uniquely_up_noiseR_sig <- Q_uniquely_up_noiseR[Q_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]
cluster1_uniquely_up_noiseR_sig <- cluster1_uniquely_up_noiseR[cluster1_uniquely_up_noiseR$FDR<0.05,][ ,c(1:3)]

# add top column
G_uniquely_up_noiseR_sig$top <- 1:nrow(G_uniquely_up_noiseR_sig)
O_uniquely_up_noiseR_sig$top <- 1:nrow(O_uniquely_up_noiseR_sig)
Q_uniquely_up_noiseR_sig$top <- 1:nrow(Q_uniquely_up_noiseR_sig)
cluster1_uniquely_up_noiseR_sig$top <- 1:nrow(cluster1_uniquely_up_noiseR_sig)

# unify S4 objects, sort by top (ascending) then FDR (ascending), keep only first unique entry for each marker
GOQ1_uniquely_up_noiseR_sig <- rbind2(G_uniquely_up_noiseR_sig,
                                      O_uniquely_up_noiseR_sig,
                                      Q_uniquely_up_noiseR_sig,
                                      cluster1_uniquely_up_noiseR_sig)

GOQ1_uniquely_up_noiseR_sig_sort <- GOQ1_uniquely_up_noiseR_sig[with(GOQ1_uniquely_up_noiseR_sig, order(top, FDR)), ]
GOQ1_uniquely_up_noiseR_sig_sort_uniq <- GOQ1_uniquely_up_noiseR_sig_sort[unique(rownames(GOQ1_uniquely_up_noiseR_sig_sort)), ]

# # de-select unannotated/ not well-characterised genes
deselected <- c()
GOQ_uniquely_up_noiseR_sig_sort_uniq_selected <- GOQ_uniquely_up_noiseR_sig_sort_uniq[!(rownames(GOQ_uniquely_up_noiseR_sig_sort_uniq) %in% deselected), ]

# top only + gene-of-interest
best_set <- GOQ_uniquely_up_noiseR_sig_sort_uniq_selected[1:50, ]

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
  annotation_row = data.frame(
    # TODO: temp trick to deal with the row-colouring problem
    cluster1.vs.2 = factor(ifelse(rownames(best_set) %in% rownames(G_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster1.vs.3 = factor(ifelse(rownames(best_set) %in% rownames(O_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster1.vs.4 = factor(ifelse(rownames(best_set) %in% rownames(Q_uniquely_up_noiseR_sig), "DE", "not DE"), levels = c("DE")),
    cluster1.vs.2.vs.3.vs.4 = factor(ifelse(rownames(best_set) %in% rownames(cluster1_uniquely_up_noiseR_sig_trim), "DE", "not DE"), levels = c("DE")),
    row.names = rownames(best_set)),
  main = paste0("Cluster: ", chosen, " ", x),
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











