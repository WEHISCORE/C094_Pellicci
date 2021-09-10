# EDA of scRNA-seq trajectory analysis
# Peter Hickey
# 2021-09-09

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(scater)
library(patchwork)

# NOTE: Have to start from this object in order to have consistency with
#       William's analyses.
sce <- readRDS(
  here("data/SCEs/C094_Pellicci.single-cell.merged.whole_cell.SCE.rds"))

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
snn_gr <- scran::buildSNNGraph(sce, use.dimred = "corrected", k = 10)
clusters <- igraph::cluster_louvain(snn_gr)
colLabels(sce) <- factor(clusters$membership)
plotReducedDim(sce, "UMAP_corrected", colour_by = "cluster") +
  plotReducedDim(sce, "UMAP_corrected", colour_by = "label") +
  plotReducedDim(sce, "UMAP_corrected", colour_by = "tissue") +
  plotReducedDim(sce, "UMAP_corrected", colour_by = "stage")

by.cluster <- aggregateAcrossCells(
  sce,
  ids = colLabels(sce),
  use_altexps = FALSE)
centroids <- reducedDim(by.cluster, "corrected")

library(TSCAN)
mst <- createClusterMST(centroids, clusters = NULL)

line.data <- reportEdges(
  by.cluster,
  mst = mst,
  clusters = NULL,
  use.dimred = "UMAP_corrected")

plotReducedDim(sce, "UMAP_corrected", colour_by = "label") +
  geom_line(data = line.data, mapping = aes(x = dim1, y = dim2, group = edge))

map.tscan <- mapCellsToEdges(sce, mst = mst, use.dimred = "corrected")
# NOTE: Manually set starting node using prior biological knowledge.
tscan.pseudo <- orderCells(map.tscan, mst, start = "3")
head(tscan.pseudo)

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have the same pseudo-time value for
# those paths anyway, so the rowMeans doesn't change anything.
common.pseudo <- rowMeans(tscan.pseudo, na.rm = TRUE)
plotReducedDim(
  sce,
  "UMAP_corrected",
  colour_by = I(common.pseudo),
  text_by = "label",
  text_colour = "red") +
  geom_line(data = line.data, mapping = aes(x = dim1, y = dim2, group = edge))

# Summary: Forcing the pseudotime to respect the clusters seems suboptimal;
# prefer the slingshot results which don't require this but recapitulates a
# similar trajectory.

# Principal curves -------------------------------------------------------------

library(slingshot)
# NOTE: Can't directly use `slingshot(sce, reducedDim = "corrected")` because
#       it gives an error:
#       "Error in solve.default(s1 + s2) :
#         system is computationally singular: reciprocal condition
#         number = 2.91591e-17"
#       This is similar to previously reported issues
#       (https://github.com/kstreet13/slingshot/issues/87 and
#       https://github.com/kstreet13/slingshot/issues/35) for which a hacky
#       workaround is to drop the last dimension.
rd <- reducedDim(sce, "corrected")
reducedDim(sce, "corrected_1") <- rd[, seq_len(ncol(rd) - 1)]
sce <- slingshot(sce, reducedDim = "corrected_1")

embedded <- embedCurves(sce, "UMAP_corrected")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord, ])

plotReducedDim(
  sce,
  dimred = "UMAP_corrected",
  colour_by = "slingPseudotime_1",
  point_alpha = 1,
  point_size = 0.5) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2) +
  plotReducedDim(
    sce,
    dimred = "UMAP_corrected",
    colour_by = "stage",
    point_alpha = 1,
    point_size = 0.5) +
  plotReducedDim(
    sce,
    dimred = "UMAP_corrected",
    colour_by = "tissue",
    point_alpha = 1,
    point_size = 0.5) +
  plotReducedDim(
    sce,
    dimred = "UMAP_corrected",
    colour_by = "cluster",
    point_alpha = 1,
    point_size = 0.5) +
  plot_layout(ncol = 2)

# Changes along a trajectory ---------------------------------------------------

pseudo <- testPseudotime(sce, pseudotime = sce$slingPseudotime_1)
pseudo[order(pseudo$p.value), ]

pseudo <- testPseudotime(
  sce,
  # NOTE: Excluding mitochondrial genes from testing
  # sce2[!rownames(sce) %in% mito_set, ],
  pseudotime = sce$slingPseudotime_1,
  block = sce$plate_number)
sorted <- pseudo[order(pseudo$p.value), ]

up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)
best <- head(rownames(up.left), 10)
plotExpression(
  sce,
  features = best,
  x = "slingPseudotime_1",
  # colour_by = "label",
  colour_by = "stage",
  show_smooth = TRUE)

up.right <- sorted[sorted$logFC > 0,]
head(up.right, 10)
best <- head(rownames(up.right), 10)
plotExpression(
  sce,
  features = best,
  x = "slingPseudotime_1",
  # colour_by = "label",
  colour_by = "stage",
  show_smooth = TRUE)

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  # colour_columns_by = "label",
  colour_columns_by = c("stage", "tissue"),
  features = head(rownames(up.right), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"))

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  # colour_columns_by = "label",
  colour_columns_by = c("stage", "tissue"),
  features = head(rownames(up.left), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"))

# NOTE: Might need some filter on aveLogCPM to avoid genes for which only a
#       handful of cells have non-zero expression.

# tradeSeq ---------------------------------------------------------------------

nonna.pseudo <- sce$slingPseudotime_1
nonna.pseudo[is.na(nonna.pseudo)] <- 0
cell.weights <- !is.na(sce$slingPseudotime_1)
storage.mode(cell.weights) <- "numeric"

# Remove genes with zero counts
sce <- sce[rowSums(counts(sce)) > 0, ]

library(tradeSeq)
set.seed(666)
# NOTE: Takes ~1 hour and can crash due to memory constraints.
fit <- fitGAM(
  counts(sce),
  pseudotime = nonna.pseudo,
  cellWeights = cell.weights,
  parallel = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(4))
saveRDS(fit, here("tmp/tradeSeq_fit.rds"))
fit <- readRDS(here("tmp/tradeSeq_fit.rds"))

res <- associationTest(fit)
res <- res[order(res$pvalue), ]
res$FDR <- p.adjust(res$pvalue, method = "BH")

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
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

# Figures for paper ------------------------------------------------------------

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

sce$CD4 <- assay(altExp(sce, "FACS"), "pseudolog")["V525_50_A_CD4_BV510", ]
sce$CD161 <- assay(altExp(sce, "FACS"), "pseudolog")["B530_30_A_CD161_FITC", ]

# pdf()
plotReducedDim(
  sce,
  "UMAP_corrected",
  colour_by = "stage",
  point_alpha = 1) +
  scale_colour_manual(values = stage_colours, name = "stage") +
  ggtitle("All cells")

plotReducedDim(
  sce,
  "UMAP_corrected",
  colour_by = "tissue",
  point_alpha = 1) +
  scale_colour_manual(values = tissue_colours, name = "tissue") +
  ggtitle("All cells")

plotReducedDim(
  sce,
  dimred = "UMAP_corrected",
  colour_by = "slingPseudotime_1",
  point_alpha = 1) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2) +
  ggtitle("All cells")

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("slingPseudotime_1", "CD4", "CD161", "stage", "tissue"),
  features = head(setdiff(rownames(res[res$FDR < 0.05, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 8,
  main = "All cells: Selected pseudotime-associated genes",
  column_annotation_colors = list(
    stage = stage_colours,
    tissue = tissue_colours))

plotReducedDim(
  sce[, sce$stage != "Unknown"],
  "UMAP_corrected",
  colour_by = "stage",
  point_alpha = 1) +
  scale_colour_manual(values = stage_colours, name = "stage") +
  ggtitle("Excluding 'Unknown'-stage cells")

plotReducedDim(
  sce[, sce$stage != "Unknown"],
  "UMAP_corrected",
  colour_by = "tissue",
  point_alpha = 1) +
  scale_colour_manual(values = tissue_colours, name = "tissue") +
  ggtitle("Excluding 'Unknown'-stage cells")

plotReducedDim(
  sce[, sce$stage != "Unknown"],
  dimred = "UMAP_corrected",
  colour_by = "slingPseudotime_1",
  point_alpha = 1) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2) +
  ggtitle("Excluding 'Unknown'-stage cells")

plotHeatmap(
  sce[, sce$stage != "Unknown"],
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("slingPseudotime_1", "CD4", "CD161", "stage", "tissue"),
  features = head(setdiff(rownames(res[res$FDR < 0.05, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 8,
  main = "Excluding 'Unknown'-stage cells: Selected pseudotime-associated genes",
  column_annotation_colors = list(
    stage = stage_colours[levels(factor(sce[, sce$stage != "Unknown"]$stage))],
    tissue = tissue_colours))

# Heatmaps sorted by CD4 and CD161 FACS expression " showing the DEGs across the gamma delta T cell stages ordered by decreasing either CD4 or KLRB1."
# NOTE: This is an example but as I explained to dan I don't think these are very useful plots
tmp <- sce[, sce$group %in% c("Thymus.S1 (CD4+/CD161-)", "Thymus.S2 (CD4-/CD161-)", "Thymus.S3 (CD4-/CD161+)", "Blood.S3 (CD4-/CD161+)")]
tmp$group <- factor(
  tmp$group,
  c("Thymus.S1 (CD4+/CD161-)", "Thymus.S2 (CD4-/CD161-)", "Thymus.S3 (CD4-/CD161+)", "Blood.S3 (CD4-/CD161+)"))

x <- read.csv(
  here("output/DEGs/excluding_donors_1-3/Thymus.S3_vs_Thymus.S1.aggregated_tech_reps.DEGs.csv.gz"))

plotHeatmap(
  tmp,
  order_columns_by = c("group", "CD4"),
  colour_columns_by = c("CD4", "group"),
  features = head(x[x$FDR < 0.05, "ENSEMBL.SYMBOL"], 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 8,
  main = "Mini-bulk DEGs: Thymus.S3 vs. Thymus.S1",
  column_annotation_colors = list(group = group_colours[levels(tmp$group)]))

plotHeatmap(
  tmp,
  order_columns_by = c("group", "CD161"),
  colour_columns_by = c("CD161", "group"),
  features = head(x[x$FDR < 0.05, "ENSEMBL.SYMBOL"], 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 8,
  main = "Mini-bulk DEGs: Thymus.S3 vs. Thymus.S1",
  column_annotation_colors = list(group = group_colours[levels(tmp$group)]))

# Mock figure sent to Dan and co (2021-09-09).
p1 <- plotReducedDim(
  sce[, sce$stage != "Unknown"],
  "UMAP_corrected",
  colour_by = "stage",
  point_alpha = 1) +
  scale_colour_manual(values = stage_colours, name = "stage")
  # ggtitle("Excluding 'Unknown'-stage cells")

p2 <- plotReducedDim(
  sce[, sce$stage != "Unknown"],
  "UMAP_corrected",
  colour_by = "tissue",
  point_alpha = 1) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
  # ggtitle("Excluding 'Unknown'-stage cells")

p3 <- plotReducedDim(
  sce[, sce$stage != "Unknown"],
  dimred = "UMAP_corrected",
  colour_by = "slingPseudotime_1",
  point_alpha = 1) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2)
  # ggtitle("Excluding 'Unknown'-stage cells")

p4 <- plotHeatmap(
  sce[, sce$stage != "Unknown"],
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("slingPseudotime_1", "CD4", "CD161", "stage", "tissue"),
  features = head(setdiff(rownames(res[res$FDR < 0.05, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 5,
  # main = "Excluding 'Unknown'-stage cells: Selected pseudotime-associated genes",
  column_annotation_colors = list(
    stage = stage_colours[levels(factor(sce[, sce$stage != "Unknown"]$stage))],
    tissue = tissue_colours),
  silent = TRUE)

p1 + p2 + p3 + p4$gtable +
  plot_annotation(
    title = "Summary of scRNA-seq data",
    subtitle = "Excluding 'Unknown'-stage cells",
    tag_levels = "a")
