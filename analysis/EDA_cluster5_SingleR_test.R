# read in different human SingleR annotation reference
hpca <- HumanPrimaryCellAtlasData()
.adf(colData(hpca)) %>%
  dplyr::count(label.main, label.fine) %>%
  dplyr::arrange(label.main) %>%
  knitr::kable()

be <- BlueprintEncodeData()
.adf(colData(be)) %>%
  dplyr::count(label.main, label.fine) %>%
  dplyr::arrange(label.main) %>%
  knitr::kable()

dice <- DatabaseImmuneCellExpressionData()
.adf(colData(dice)) %>%
  dplyr::count(label.main, label.fine) %>%
  dplyr::arrange(label.main) %>%
  knitr::kable()

mi <- MonacoImmuneData()
.adf(colData(mi)) %>%
  dplyr::count(label.main, label.fine) %>%
  dplyr::arrange(label.main) %>%
  knitr::kable()

nh <- NovershternHematopoieticData()
.adf(colData(nh)) %>%
  dplyr::count(label.main, label.fine) %>%
  dplyr::arrange(label.main) %>%
  knitr::kable()

##########################################

# set reference
library(SingleR)
dice <- DatabaseImmuneCellExpressionData()
ref <- dice
labels_fine <- ref$label.fine

# define colour of collapsed labels
polychrome_colors <- c(
  Polychrome::glasbey.colors(),
  Polychrome::kelly.colors(),
  Polychrome::green.armytage.colors(),
  Polychrome::palette36.colors(),
  Polychrome::alphabet.colors(),
  Polychrome::light.colors(),
  Polychrome::dark.colors()
)
label_fine_collapsed_colours <- setNames(
  c(
  polychrome_colors[1:nlevels(factor(labels_fine))],
  "black"),
  c(levels(factor(labels_fine)), "other"))

# SingleR (cluster level)
pred_cluster_fine <- SingleR(
  test = sce,
  ref = ref[!grepl("^mt|^Rps|^Rpl", rownames(ref)), ],
  labels = labels_fine,
  cluster = sce$cluster,
  BPPARAM = bpparam())
sce$label_cluster_fine <- factor(pred_cluster_fine$pruned.labels[sce$cluster])
sce$label_cluster_fine_collapsed <- .collapseLabel(
  sce$label_cluster_fine,
  sce$batch)
sce$label_fine_collapsed_colours <- label_fine_collapsed_colours[
  as.character(sce$label_cluster_fine)]
umap_df <- makePerCellDF(sce)
umap_df$label_cluster_fine_collapsed <- sce$label_cluster_fine_collapsed
# tabulate sumup the annotation outcome
tabyl(
  data.frame(label.fine = sce$label_cluster_fine, cluster = sce$cluster),
  cluster,
  label.fine) %>%
  knitr::kable(
    caption = "Cluster-level assignments using the fine labels of the reference.")

# UMAP sumup
p1 <- ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
  geom_point(
    aes(colour = cluster),
    alpha = 1,
    size = 0.25) +
  scale_fill_manual(values = cluster_colours) +
  scale_colour_manual(values = cluster_colours) +
  theme_cowplot(font_size = 10) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  ggtitle("Clusters")
bg <- dplyr::select(umap_df, -label_cluster_fine_collapsed)
p2 <- ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
  geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
  geom_point(
    aes(colour = label_cluster_fine_collapsed),
    alpha = 1,
    size = 0.25) +
  scale_fill_manual(values = label_fine_collapsed_colours) +
  scale_colour_manual(values = label_fine_collapsed_colours) +
  theme_cowplot(font_size = 10) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  facet_wrap(~ label_cluster_fine_collapsed, ncol = 2) +
  guides(colour = FALSE) +
  ggtitle("'fine' cluster-level label")
p1 + p2 + plot_layout(widths = c(1, 2))







##########################################

all_markers <- metadata(pred_cluster_fine)$de.genes
# Get the top-50 marker genes
top_markers_1 <- head(all_markers[["T_cell:gamma-delta"]]$`T_cell:CD8+`, 50)
top_markers_2 <- head(all_markers[["T_cell:CD8+"]]$`T_cell:gamma-delta`, 50)
top_markers <- c(top_markers_1, top_markers_2)

# heatmap sumup
library(scater)
plotHeatmap(
  sce,
  features = top_markers,
  columns = order(
    sce$cluster,
    sce$label_cluster_fine_collapsed,
    sce$sample_name,
    sce$plate_number,
    sce$post_hoc_sample_gate,
    sce$tissue,
    sce$donor,
    sce$sample_type),
  colour_columns_by = c(
    "cluster",
    "label_cluster_fine_collapsed",
    "sample_name",
    "plate_number",
    "post_hoc_sample_gate",
    "tissue",
    "donor",
    "sample_type"),
  cluster_cols = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  show_colnames = FALSE,
  main = "Intermediate monocytes vs. gamma-delta T cells",
  column_annotation_colors = list(
     cluster = cluster_colours,
     label_cluster_fine_collapsed = label_fine_collapsed_colours,
     sample_name = sample_name_colours,
     plate_number = plate_number_colours,
     post_hoc_sample_gate = sample_gate_colours,
     tissue = tissue_colours,
     donor = donor_colours,
     sample_type = sample_type_colours),
  fontsize = 5)
