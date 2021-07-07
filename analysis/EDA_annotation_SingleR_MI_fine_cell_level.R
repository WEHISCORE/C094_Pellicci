library(SingleCellExperiment)
library(here)
library(ggplot2)
library(BiocParallel)
library(pheatmap)
library(dplyr)
library(tidyr)
source(here("code", "helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.cells_selected.SCE.rds"))

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

# clustering
set.seed(4759)
snn_gr <- buildSNNGraph(sce, use.dimred = "corrected")
clusters <- igraph::cluster_louvain(snn_gr)
sce$cluster <- factor(clusters$membership)
cluster_colours <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce$cluster))],
  levels(sce$cluster))
sce$colours$cluster_colours <- cluster_colours[sce$cluster]

# SingleR annotation (ref: mi, fine label, cell level annotaion)
library(SingleR)

mi <- MonacoImmuneData()

pred_mi_cell_fine <- SingleR(
test = sce,
ref = mi[!grepl("^mt|^Rps|^Rpl", rownames(mi)), ],
labels = mi$label.fine,
BPPARAM = SerialParam())

tabyl(data.frame(label.fine = pred_mi_cell_fine$labels), label.fine) %>%
adorn_pct_formatting(digits = 1) %>%
dplyr::arrange(desc(n)) %>%
knitr::kable(
caption = "Cell label assignments using the fine labels of the `MonacoImmuneData` reference data.")

# heatmap (group score)
stopifnot(identical(rownames(pred_mi_cell_fine), colnames(sce)))
plotScoreHeatmap(
  pred_mi_cell_fine,
  annotation_col = data.frame(
    cluster = sce$cluster,
    sample = sce$sample,
    stage = sce$stage,
    plate_number = sce$plate_number,
    tissue = sce$tissue,
    donor = sce$donor,
    row.names = rownames(pred_mi_cell_fine)),
  annotation_colors = list(
    cluster = cluster_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours,
    tissue = tissue_colours,
    donor = donor_colours),
  labels.use = unique(pred_mi_cell_fine$labels),
  max.labels = Inf,
  normalize = TRUE,
  show.labels = FALSE,
  fontsize = 6)

# heatmap (group by cluster)
stopifnot(identical(rownames(pred_mi_cell_fine), colnames(sce)))
plotScoreHeatmap(
  pred_mi_cell_fine,
  annotation_col = data.frame(
    cluster = sce$cluster,
    sample = sce$sample,
    stage = sce$stage,
    plate_number = sce$plate_number,
    tissue = sce$tissue,
    donor = sce$donor,
    row.names = rownames(pred_mi_cell_fine)),
  annotation_colors = list(
    cluster = cluster_colours,
    sample = sample_colours,
    stage = stage_colours,
    plate_number = plate_number_colours,
    tissue = tissue_colours,
    donor = donor_colours),
  cells.order = sce$cluster,
  labels.use = unique(pred_mi_cell_fine$labels),
  max.labels = Inf,
  normalize = TRUE,
  show.labels = FALSE,
  fontsize = 6)

# overlays the normalized scores for the 10 labels with the largest maximum scores on the UMAP plot
plotScoreReducedDim(
  pred_mi_cell_fine, sce,
  max.labels = 10,
  ncol = 2,
  dimred = "UMAP_corrected")
