
library(here)
library(SingleCellExperiment)

# read-in
sce <- readRDS(here("data", "SCEs", "C094_Pellicci_annotate.SCE.rds"))

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

# exclusion of cluster 3 and adjust factor levels
sce <- sce[,sce$cluster!="3"]
colData(sce) <- droplevels(colData(sce))


# Reprocessing
sce$batch <- sce$plate_number
var_fit <- modelGeneVarWithSpikes(sce, "ERCC", block = sce$batch)
hvg <- getTopHVGs(var_fit, var.threshold = 0)
hvg <- setdiff(hvg, c(ribo_set, mito_set, pseudogene_set))

set.seed(67726)
sce <- denoisePCA(
  sce,
  var_fit,
  subset.row = hvg,
  BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE))

set.seed(853)
sce <- runUMAP(sce, dimred = "PCA")

set.seed(4759)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA")
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
cluster_colours <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce$cluster))],
  levels(sce$cluster))
sce$colours$cluster_colours <- cluster_colours[sce$cluster]

# UMAP plot
p1 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = cluster), size = 0.25) +
  scale_colour_manual(values = cluster_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

p2 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = sample), size = 0.25) +
  scale_colour_manual(values = sample_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

p3 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = stage), size = 0.25) +
  scale_colour_manual(values = stage_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

p4 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = plate_number), size = 0.25) +
  scale_colour_manual(values = plate_number_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

p5 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = tissue), size = 0.25) +
  scale_colour_manual(values = tissue_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

p6 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = donor), size = 0.25) +
  scale_colour_manual(values = donor_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

(p1 | p2) / (p3 | p4) / (p5 | p6)

# stacked barplot
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

(p1 | p2) / (p3 | p4) / (p5 | p6)

# MNN correction
library(batchelor)
set.seed(1819)

mnn_out <- fastMNN(
  multiBatchNorm(sce, batch = sce$batch),
  batch = sce$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce, "PCA")),
  auto.merge = FALSE,
  merge.order = list(list("LCE513", "LCE514", "LCE509"), list("LCE508", "LCE511", "LCE512")),
  subset.row = hvg)

tab <- metadata(mnn_out)$merge.info$lost.var
knitr::kable(
  100 * tab,
  digits = 1,
  caption = "Percentage of estimated biological variation lost within each plate at each step of the merge (manual). Ideally, all these values should be small (e.g., < 5%).")

reducedDim(sce, "corrected") <- reducedDim(mnn_out, "corrected")

# generate UMAP
set.seed(1248)
sce <- runUMAP(sce, dimred = "corrected", name = "UMAP_corrected")

# re-clustering after each MNN correction
set.seed(4759)

snn_gr <- buildSNNGraph(sce, use.dimred = "corrected")
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
cluster_colours <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce$cluster))],
  levels(sce$cluster))
sce$colours$cluster_colours <- cluster_colours[sce$cluster]


# UMAP: before vs after
p1 <- plotReducedDim(sce, "UMAP", colour_by = "plate_number") +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p2 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "plate_number") +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")

p3 <- plotReducedDim(sce, "UMAP", colour_by = "sample") +
  scale_colour_manual(values = sample_colours, name = "sample")
p4 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "sample") +
  scale_colour_manual(values = sample_colours, name = "sample")

p5 <- plotReducedDim(sce, "UMAP", colour_by = "cluster") +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p6 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "cluster") +
  scale_colour_manual(values = cluster_colours, name = "cluster")

p7 <- plotReducedDim(sce, "UMAP", colour_by = "stage") +
  scale_colour_manual(values = stage_colours, name = "stage")
p8 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "stage") +
  scale_colour_manual(values = stage_colours, name = "stage")

p9 <- plotReducedDim(sce, "UMAP", colour_by = "tissue") +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p10 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "tissue") +
  scale_colour_manual(values = tissue_colours, name = "tissue")

p11 <- plotReducedDim(sce, "UMAP", colour_by = "donor") +
  scale_colour_manual(values = donor_colours, name = "donor")
p12 <- plotReducedDim(sce, "UMAP_corrected", colour_by = "donor") +
  scale_colour_manual(values = donor_colours, name = "donor")

p1 + p2 + p3 + p4 +
  p5 + p6 + p7 + p8 +
  p9 + p10 + p11 + p12 +
  plot_layout(ncol = 2, guides = "collect")


# UMAP - breakdown by plate
umap_df <- makePerCellDF(sce)
bg <- dplyr::select(umap_df, -plate_number)

plot_grid(

  ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
    geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
    geom_point(aes(colour = plate_number), alpha = 1, size = 0.5) +
    scale_fill_manual(values = plate_number_colours, name = "plate_number") +
    scale_colour_manual(values = plate_number_colours, name = "plate_number") +
    theme_cowplot(font_size = 10) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    facet_wrap(~plate_number, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    guides(colour = FALSE),

  ggplot(aes(x = UMAP_corrected.1, y = UMAP_corrected.2), data = umap_df) +
    geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
    geom_point(aes(colour = plate_number), alpha = 1, size = 0.5) +
    scale_fill_manual(values = plate_number_colours, name = "plate_number") +
    scale_colour_manual(values = plate_number_colours, name = "plate_number") +
    theme_cowplot(font_size = 10) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    facet_wrap(~plate_number, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    guides(colour = FALSE),

  ncol= 2,
  align ="h"
)


# UMAP breakdown by stage
umap_df <- makePerCellDF(sce)
bg <- dplyr::select(umap_df, -stage)

plot_grid(

  ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
    geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
    geom_point(aes(colour = stage), alpha = 1, size = 0.5) +
    scale_fill_manual(
      values = stage_colours,
      name = "stage") +
    scale_colour_manual(
      values = stage_colours,
      name = "stage") +
    theme_cowplot(font_size = 10) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    facet_wrap(~stage, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    guides(colour = FALSE),

  ggplot(aes(x = UMAP_corrected.1, y = UMAP_corrected.2), data = umap_df) +
    geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
    geom_point(aes(colour = stage), alpha = 1, size = 0.5) +
    scale_fill_manual(
      values = stage_colours,
      name = "stage") +
    scale_colour_manual(
      values = stage_colours,
      name = "stage") +
    theme_cowplot(font_size = 10) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    facet_wrap(~stage, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    guides(colour = FALSE),

  ncol= 2,
  align ="h"
)


# UMAP breakdown by tissue
umap_df <- makePerCellDF(sce)
bg <- dplyr::select(umap_df, -tissue)

plot_grid(

  ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
    geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
    geom_point(aes(colour = tissue), alpha = 1, size = 0.5) +
    scale_fill_manual(values = tissue_colours, name = "tissue") +
    scale_colour_manual(values = tissue_colours, name = "tissue") +
    theme_cowplot(font_size = 10) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    facet_wrap(~tissue, ncol = 2) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    guides(colour = FALSE),

  ggplot(aes(x = UMAP_corrected.1, y = UMAP_corrected.2), data = umap_df) +
    geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
    geom_point(aes(colour = tissue), alpha = 1, size = 0.5) +
    scale_fill_manual(values = tissue_colours, name = "tissue") +
    scale_colour_manual(values = tissue_colours, name = "tissue") +
    theme_cowplot(font_size = 10) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    facet_wrap(~tissue, ncol = 2) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    guides(colour = FALSE),

  ncol= 2,
  align ="h"
)


# integrated merge to SCE + clustering
reducedDim(sce, "corrected") <- reducedDim(mnn_out, "corrected")

set.seed(1248)
sce <- runUMAP(sce, dimred = "corrected", name = "UMAP_corrected")

set.seed(4759)
snn_gr <- buildSNNGraph(sce, use.dimred = "corrected")
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
cluster_colours <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce$cluster))],
  levels(sce$cluster))
sce$cluster_colours <- cluster_colours[sce$cluster]


# UMAP (after merge)
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

(p1 | p2) / (p3 | p4) / (p5 | p6)


# Stacked barplot (after merge)
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

(p1 | p2) / (p3 | p4) / (p5 | p6)



# Marker gene detection (unique)
sce$block <- paste0(sce$plate_number)

markers_uniquely_up <- findMarkers(
  sce,
  groups = sce$cluster,
  block = sce$block,
  pval.type = "all",
  direction = "up")

# cluster 1
chosen <- "1"
interesting_uniquely_up <- markers_uniquely_up[[chosen]]

best_set <- interesting_uniquely_up[1:25, ]

plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
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
  main = paste0("Cluster: ", chosen),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



# cluster 3
chosen <- "3"
interesting_uniquely_up <- markers_uniquely_up[[chosen]]

best_set <- interesting_uniquely_up[1:25, ]

plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
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
  main = paste0("Cluster: ", chosen),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)



# cluster 2
chosen <- "2"
interesting_uniquely_up <- markers_uniquely_up[[chosen]]

best_set <- interesting_uniquely_up[1:25, ]

plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$sample,
    sce$stage,
    sce$plate_number),
  colour_columns_by = c(
    "cluster",
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
  main = paste0("Cluster: ", chosen),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)
# NOTE: seems to be all insignificant marker (tried even top100, but still, all insignificant)

# cluster 2 (union any)

markers_union_any <- findMarkers(
  sce,
  groups = sce$cluster,
  block = sce$block,
  pval.type = "any",
  direction = "any")

chosen <- "2"
interesting <- markers_union_any[[chosen]]

.adf(head(interesting, 10)) %>%
  tibble::rownames_to_column("Gene") %>%
  knitr::kable(
    caption = paste0(
      "First 10 marker genes for cluster ",
      chosen,
      " across comparisons."))

best_set <- interesting[interesting$Top <= 10, ]

plotHeatmap(
  sce,
  features = rownames(best_set),
  columns = order(
    sce$cluster,
    sce$sample,
    sce$stage,
    sce$plate_number),
  order_columns_by = c(
    "cluster",
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
  main = paste0("Cluster: ", chosen),
  column_annotation_colors = list(
    Sig = c("Yes" = "red", "No" = "lightgrey"),
    cluster = cluster_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 6)
