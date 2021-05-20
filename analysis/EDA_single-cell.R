#' # Setup

library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(cowplot)
library(patchwork)

sce <- readRDS(here("data/SCEs/C094_Pellicci.preprocessed.SCE.rds"))

sample_name_colours <- setNames(
  unique(sce$colours$sample_name_colours),
  unique(names(sce$colours$sample_name_colours)))
plate_number_colours <- setNames(
  unique(sce$colours$plate_number_colours),
  unique(names(sce$colours$plate_number_colours)))

facs_markers <- rownames(altExp(sce, "FACS"))

#' # Summary of preprocessed data

plotColData(sce, "plate_number", "sample_name", colour_by = "sample_name") +
  scale_colour_manual(values = sample_name_colours, name = "sample_name")

umap_df <- makePerCellDF(sce)
bg <- dplyr::select(umap_df, -plate_number, -sample_name)
ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
  geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
  geom_point(aes(colour = sample_name), alpha = 1, size = 0.5) +
  scale_fill_manual(values = sample_name_colours, name = "sample_name") +
  scale_colour_manual(values = sample_name_colours, name = "sample_name") +
  theme_cowplot(font_size = 10) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  facet_wrap(~ plate_number, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

p <- lapply(facs_markers, function(f) {
  plotUMAP(
    sce,
    colour_by = I(assay(altExp(sce, "FACS"), "pseudolog")[f, ]),
    point_size = 0.5) +
    ggtitle(f) +
    theme_cowplot(font_size = 8)
})
#+ fig.asp = 4 / 3
wrap_plots(p, ncol = 3)

p <- lapply(facs_markers, function(f) {
  plotExpression(
    sce,
    f,
    x = "sample_name",
    colour_by = "sample_name",
    exprs_values = "pseudolog") +
    theme(axis.text.x = element_blank()) +
    scale_colour_manual(values = sample_name_colours, name = "sample_name") +
    theme_cowplot(font_size = 6)
})
#+ fig.asp = 4 / 3
wrap_plots(p, ncol = 3) + plot_layout(guides = "collect")

#' # Clustering

set.seed(4759)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA")
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
cluster_colours <- setNames(
  palette.colors(nlevels(sce$cluster), "Tableau 10"),
  levels(sce$cluster))
sce$colours$cluster_colours <- cluster_colours[sce$cluster]

plotUMAP(sce, colour_by = "cluster", text_by = "cluster", point_alpha = 1) +
  scale_colour_manual(values = cluster_colours, name = "cluster")

plotColData(sce, "cluster", "sample_name", colour_by = "sample_name") +
  scale_colour_manual(values = sample_name_colours, name = "sample_name")

plotColData(sce, "cluster", "plate_number", colour_by = "plate_number") +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")

p <- lapply(facs_markers, function(f) {
  plotExpression(
    sce,
    f,
    x = "cluster",
    colour_by = "cluster",
    exprs_values = "pseudolog") +
    theme(axis.text.x = element_blank()) +
    scale_colour_manual(values = cluster_colours, name = "cluster") +
    theme_cowplot(font_size = 6)
})
#+ fig.asp = 4 / 3
wrap_plots(p, ncol = 3) + plot_layout(guides = "collect")

#' # Cluster marker genes

out <- pairwiseTTests(
  sce,
  sce$cluster,
  direction = "up")
top_markers <- getTopMarkers(
  out$statistics,
  out$pairs,
  n = 10,
  pairwise = FALSE,
  pval.type = "some")
features <- unique(unlist(top_markers))
#+ fig.asp = 1.5
plotHeatmap(
  sce,
  features,
  order_columns_by = c("cluster", "sample_name"),
  cluster_rows = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  zlim = c(-3, 3),
  column_annotation_colors = list(
    cluster = cluster_colours,
    sample_name = sample_name_colours,
    plate_number = plate_number_colours),
  main = "Upregulated cluster marker genes",
  fontsize = 8)

#' # FACS data

tmp <- altExp(sce, "FACS")
colData(tmp) <- colData(sce)
assay(tmp, "pseudolog")[is.na(assay(tmp, "pseudolog"))] <- 0
plotHeatmap(
  tmp,
  facs_markers,
  order_columns_by = c("cluster", "sample_name"),
  cluster_rows = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  column_annotation_colors = list(
    cluster = cluster_colours,
    sample_name = sample_name_colours,
    plate_number = plate_number_colours),
  exprs_values = "pseudolog",
  main = "FACS markers",
  fontsize = 8)

plotHeatmap(
  tmp,
  facs_markers,
  colour_columns_by = c("cluster", "sample_name"),
  cluster_rows = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  center = TRUE,
  column_annotation_colors = list(
    cluster = cluster_colours,
    sample_name = sample_name_colours,
    plate_number = plate_number_colours),
  exprs_values = "pseudolog",
  main = "FACS markers",
  fontsize = 8)
