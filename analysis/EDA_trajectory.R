# EDA of scRNA-seq trajectory analysis
# Peter Hickey
# 2021-08-17

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(scater)
library(patchwork)

sce <- readRDS(
  here("data/SCEs/C094_Pellicci.single-cell.cell_selected.SCE.rds"))

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
protein_coding_gene_set <- rownames(sce)[
  any(grepl("protein_coding", rowData(sce)$ENSEMBL.GENEBIOTYPE))]

# Cluster-based minimum spanning tree ------------------------------------------

# NOTE: Need more than the 'default' number of clusters to do something
#       (potentially) interesting with TSCAN.
set.seed(4759)
snn_gr <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 5)
clusters <- igraph::cluster_louvain(snn_gr)
colLabels(sce) <- factor(clusters$membership)
plotUMAP(sce, colour_by = "label") +
  plotUMAP(sce, colour_by = "tissue") +
  plotUMAP(sce, colour_by = "stage")

by.cluster <- aggregateAcrossCells(
  sce,
  ids = colLabels(sce),
  use_altexps = FALSE)
centroids <- reducedDim(by.cluster, "PCA")

library(TSCAN)
mst <- createClusterMST(centroids, clusters = NULL)

line.data <- reportEdges(
  by.cluster,
  mst = mst,
  clusters = NULL,
  use.dimred = "UMAP")

plotUMAP(
  sce,
  colour_by = "label") +
  geom_line(data = line.data, mapping = aes(x = dim1, y = dim2, group = edge))

map.tscan <- mapCellsToEdges(sce, mst = mst, use.dimred = "PCA")
# NOTE: Manually set starting node using prior biological knowledge.
tscan.pseudo <- orderCells(map.tscan, mst, start = "5")
head(tscan.pseudo)

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have the same pseudo-time value for
# those paths anyway, so the rowMeans doesn't change anything.
common.pseudo <- rowMeans(tscan.pseudo, na.rm = TRUE)
plotUMAP(
  sce,
  colour_by = I(common.pseudo),
  text_by = "label",
  text_colour = "red") +
  geom_line(data = line.data, mapping = aes(x = dim1, y = dim2, group = edge))

# Summary: Forcing the pseudotime to respect the clusters seems suboptimal;
# prefer the slingshot results which don't require this but recapitulates a
# similar trajectory.

# Principal curves -------------------------------------------------------------

library(slingshot)
sce.sling <- slingshot(sce, reducedDim = "PCA")
head(sce.sling$slingPseudotime_1)

embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord, ])

plotUMAP(
  sce.sling,
  colour_by = "slingPseudotime_1",
  point_alpha = 1,
  point_size = 0.5) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2) +
  plotUMAP(sce.sling, colour_by = "stage", point_alpha = 1, point_size = 0.5) +
  plotUMAP(sce.sling, colour_by = "tissue", point_alpha = 1, point_size = 0.5) +
  plot_layout(ncol = 2)

# Changes along a trajectory ---------------------------------------------------

pseudo <- testPseudotime(sce, pseudotime = sce.sling$slingPseudotime_1)
pseudo[order(pseudo$p.value), ]

pseudo <- testPseudotime(
  sce.sling,
  # NOTE: Excluding mitochondrial genes from testing
  # sce2[!rownames(sce) %in% mito_set, ],
  pseudotime = sce.sling$slingPseudotime_1,
  block = sce.sling$plate_number)
sorted <- pseudo[order(pseudo$p.value), ]

up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)
best <- head(rownames(up.left), 10)
plotExpression(
  sce.sling,
  features = best,
  x = "slingPseudotime_1",
  # colour_by = "label",
  colour_by = "stage",
  show_smooth = TRUE)

up.right <- sorted[sorted$logFC > 0,]
head(up.right, 10)
best <- head(rownames(up.right), 10)
plotExpression(
  sce.sling,
  features = best,
  x = "slingPseudotime_1",
  # colour_by = "label",
  colour_by = "stage",
  show_smooth = TRUE)

plotHeatmap(
  sce.sling,
  order_columns_by = "slingPseudotime_1",
  # colour_columns_by = "label",
  colour_columns_by = c("stage", "tissue"),
  features = head(rownames(up.right), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"))

plotHeatmap(
  sce.sling,
  order_columns_by = "slingPseudotime_1",
  # colour_columns_by = "label",
  colour_columns_by = c("stage", "tissue"),
  features = head(rownames(up.left), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"))

# TODO: Might need some filter on aveLogCPM to avoid genes for which only a
#       handful of cells have non-zero expression.

# tradeSeq ---------------------------------------------------------------------

nonna.pseudo <- sce.sling$slingPseudotime_1
nonna.pseudo[is.na(nonna.pseudo)] <- 0
cell.weights <- !is.na(sce.sling$slingPseudotime_1)
storage.mode(cell.weights) <- "numeric"

library(tradeSeq)
set.seed(666)
# NOTE: Takes ~1 hour.
# fit <- fitGAM(
#   counts(sce),
#   pseudotime = nonna.pseudo,
#   cellWeights = cell.weights,
#   parallel = TRUE,
#   BPPARAM = BiocParallel::MulticoreParam(8))
# saveRDS(fit, here("tmp/tradeSeq_fit.rds"))
fit <- readRDS(here("tmp/tradeSeq_fit.rds"))

res <- associationTest(fit)
res <- res[order(res$pvalue), ]
res$FDR <- p.adjust(res$pvalue, method = "BH")

plotHeatmap(
  sce.sling,
  order_columns_by = "slingPseudotime_1",
  # colour_columns_by = "label",
  colour_columns_by = c("stage", "tissue"),
  features = head(setdiff(rownames(res[res$FDR < 0.05, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 8)

write.csv(
  res,
  here("tmp/C094_Pellicci.tradeSeq_results.csv"),
  row.names = TRUE,
  quote = FALSE)
