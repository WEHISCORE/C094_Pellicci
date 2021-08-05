


# setup data
library(SingleCellExperiment)
library(here)
library(scater)
library(scran)
library(ggplot2)
library(cowplot)
library(edgeR)
library(Glimma)
library(BiocParallel)
library(patchwork)
library(pheatmap)
library(janitor)
library(distill)
source(here("code", "helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.single-cell.cell_selection.SCE.rds"))

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

# also define group (i.e. tissue.stage)
sce$group <- factor(paste0(sce$tissue, ".", sce$stage))
group_colours <- setNames(
  Polychrome::kelly.colors(nlevels(sce$group) + 1)[-1],
  levels(sce$group))
sce$colours$group_colours <- group_colours[sce$group]

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

### Re-processing (whole cell)
sce$batch <- sce$plate_number
var_fit <- modelGeneVarWithSpikes(sce, "ERCC", block = sce$batch)
hvg <- getTopHVGs(var_fit, var.threshold = 0)
hvg <- intersect(hvg, protein_coding_gene_set)

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

# UMAP
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

p7 <- ggcells(sce, aes(x = UMAP.1, y = UMAP.2)) +
  geom_point(aes(colour = group), size = 0.25) +
  scale_colour_manual(values = group_colours) +
  theme_cowplot(font_size = 8) +
  xlab("Dimension 1") +
  ylab("Dimension 2")

(p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | plot_spacer())

#stacked bar
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

# Data integration - test























```{r}
# backup before test
sce0 <- sce
```

```{r}
sce0$batch <- sce0$plate_number
var_fit <- modelGeneVarWithSpikes(sce0, "ERCC", block = sce0$batch)
hvg <- getTopHVGs(var_fit, var.threshold = 0)
hvg <- setdiff(hvg, c(ribo_set, mito_set, pseudogene_set))
```

```{r, results = "hide"}
library(batchelor)
set.seed(1819)
# auto merge
# NOTE: auto merge does not make sense; this seems like merging plate with adjacent number first (i.e. LCE513 + LCE514, + LCE509, LCE508, LCE511, LCE512); remark: all plates were in the same sequencing run, i.e. NN215
mnn_out_1 <- fastMNN(
  multiBatchNorm(sce0, batch = sce0$batch),
  batch = sce0$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce0, "PCA")),
  auto.merge = TRUE,
  subset.row = hvg)
# manual merge 1
# NOTE: order defined based on the similarities in cell number between plates for each`sample_name`
mnn_out_2 <- fastMNN(
  multiBatchNorm(sce0, batch = sce0$batch),
  batch = sce0$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce0, "PCA")),
  auto.merge = FALSE,
  merge.order = list(list("LCE508", "LCE511", "LCE512", "LCE513", "LCE514"), "LCE509"),
  subset.row = hvg)
# manual merge 2
# NOTE: to reduced var loss of LCE509, which contained different proportion of cells from each sample, I decide to merge LCE509 after the merge of LCE503 + LCE504 (as indicated in the auto-merge)
mnn_out_3 <- fastMNN(
  multiBatchNorm(sce0, batch = sce0$batch),
  batch = sce0$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce0, "PCA")),
  auto.merge = FALSE,
  merge.order = list(list("LCE513", "LCE514", "LCE509"), list("LCE508", "LCE511", "LCE512")),
  subset.row = hvg)
```

One useful diagnostic of the MNN algorithm is the proportion of variance within each batch that is lost during MNN correction^[Specifically, this refers to the within-batch variance that is removed during orthogonalization with respect to the average correction vector at each merge step.].
Large proportions of lost variance ($>10 \%$) suggest that correction is removing genuine biological heterogeneity.
This would occur due to violations of the assumption of orthogonality between the batch effect and the biological subspace [@haghverdi2018batch].

Among the three merge orders, the *manual merge 2* (Table \@ref(tab:tab3)) lead to the least loss of biological variance from the dataset.

```{r tab1}
tab1 <- metadata(mnn_out_1)$merge.info$lost.var
knitr::kable(
  100 * tab1,
  digits = 1,
  caption = "Percentage of estimated biological variation lost within each plate at each step of the merge (auto). Ideally, all these values should be small (e.g., < 5%).")
```

```{r tab2}
tab2 <- metadata(mnn_out_2)$merge.info$lost.var
knitr::kable(
  100 * tab2,
  digits = 1,
  caption = "Percentage of estimated biological variation lost within each plate at each step of the merge (manual 1). Ideally, all these values should be small (e.g., < 5%).")
```

```{r tab3}
tab3 <- metadata(mnn_out_3)$merge.info$lost.var
knitr::kable(
  100 * tab3,
  digits = 1,
  caption = "Percentage of estimated biological variation lost within each plate at each step of the merge (manual 2). Ideally, all these values should be small (e.g., < 5%).")
```

```{r}
reducedDim(sce0, "corrected_1") <- reducedDim(mnn_out_1, "corrected")
reducedDim(sce0, "corrected_2") <- reducedDim(mnn_out_2, "corrected")
reducedDim(sce0, "corrected_3") <- reducedDim(mnn_out_3, "corrected")
# generate UMAP
set.seed(18191)
sce0 <- runUMAP(sce0, dimred = "PCA", name = "UMAP")
uncorrected_umap <- cbind(
  data.frame(
    x = reducedDim(sce0, "UMAP")[, 1],
    y = reducedDim(sce0, "UMAP")[, 2]),
  as.data.frame(colData(sce0)))
set.seed(11901)
sce0 <- runUMAP(sce0, dimred = "corrected_1", name = "UMAP_corrected_1")
corrected_umap_1 <- cbind(
  data.frame(
    x = reducedDim(sce0, "UMAP_corrected_1")[, 1],
    y = reducedDim(sce0, "UMAP_corrected_1")[, 2]),
  as.data.frame(colData(sce0)))
set.seed(11901)
sce0 <- runUMAP(sce0, dimred = "corrected_2", name = "UMAP_corrected_2")
corrected_umap_2 <- cbind(
  data.frame(
    x = reducedDim(sce0, "UMAP_corrected_2")[, 1],
    y = reducedDim(sce0, "UMAP_corrected_2")[, 2]),
  as.data.frame(colData(sce0)))
set.seed(11901)
sce0 <- runUMAP(sce0, dimred = "corrected_3", name = "UMAP_corrected_3")
corrected_umap_3 <- cbind(
  data.frame(
    x = reducedDim(sce0, "UMAP_corrected_3")[, 1],
    y = reducedDim(sce0, "UMAP_corrected_3")[, 2]),
  as.data.frame(colData(sce0)))
# re-clustering after each MNN correction
set.seed(4759)
snn_gr_1 <- buildSNNGraph(sce0, use.dimred = "corrected_1")
clusters_1 <- igraph::cluster_louvain(snn_gr_1)
sce0$cluster_1 <- factor(clusters_1$membership)
cluster_colours_1 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_1))],
  levels(sce0$cluster_1))
sce0$colours$cluster_colours_1 <- cluster_colours_1[sce0$cluster_1]
set.seed(4759)
snn_gr_2 <- buildSNNGraph(sce0, use.dimred = "corrected_2")
clusters_2 <- igraph::cluster_louvain(snn_gr_2)
sce0$cluster_2 <- factor(clusters_2$membership)
cluster_colours_2 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_2))],
  levels(sce0$cluster_2))
sce0$colours$cluster_colours_2 <- cluster_colours_2[sce0$cluster_2]
set.seed(4759)
snn_gr_3 <- buildSNNGraph(sce0, use.dimred = "corrected_3")
clusters_3 <- igraph::cluster_louvain(snn_gr_3)
sce0$cluster_3 <- factor(clusters_3$membership)
cluster_colours_3 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce0$cluster_3))],
  levels(sce0$cluster_3))
sce0$colours$cluster_colours_3 <- cluster_colours_3[sce0$cluster_3]
```

Figure \@ref(fig:batch-correction-test) shows an overview of comparisons between the unmerged and merged (left to right: unmerged, auto, manual merge 1, and manual merge 2) data.

```{r batch-correction-test, fig.asp = 1.5, fig.cap = "Comparison between batch-uncorrected data (leftmost column) and -corrected data by auto merge (2nd column) and two different manual merge orders (3rd and rightmost columns, respectively).", layout="l-page"}
p1 <- plotReducedDim(sce0, "UMAP", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p2 <- plotReducedDim(sce0, "UMAP_corrected_1", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p3 <- plotReducedDim(sce0, "UMAP_corrected_2", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p4 <- plotReducedDim(sce0, "UMAP_corrected_3", colour_by = "plate_number", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = plate_number_colours, name = "plate_number")
p5 <- plotReducedDim(sce0, "UMAP", colour_by = "sample_name", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_name_colours, name = "sample_name")
p6 <- plotReducedDim(sce0, "UMAP_corrected_1", colour_by = "sample_name", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_name_colours, name = "sample_name")
p7 <- plotReducedDim(sce0, "UMAP_corrected_2", colour_by = "sample_name", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_name_colours, name = "sample_name")
p8 <- plotReducedDim(sce0, "UMAP_corrected_3", colour_by = "sample_name", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_name_colours, name = "sample_name")
p9 <- plotReducedDim(sce0, "UMAP", colour_by = "cluster", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours, name = "cluster")
p10 <- plotReducedDim(sce0, "UMAP_corrected_1", colour_by = "cluster_1", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_1, name = "cluster_1")
p11 <- plotReducedDim(sce0, "UMAP_corrected_2", colour_by = "cluster_2", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_2, name = "cluster_2")
p12 <- plotReducedDim(sce0, "UMAP_corrected_3", colour_by = "cluster_3", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = cluster_colours_3, name = "cluster_3")
p13 <- plotReducedDim(sce0, "UMAP", colour_by = "post_hoc_sample_gate", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_gate_colours, name = "post_hoc_sample_gate")
p14 <- plotReducedDim(sce0, "UMAP_corrected_1", colour_by = "post_hoc_sample_gate", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_gate_colours, name = "post_hoc_sample_gate")
p15 <- plotReducedDim(sce0, "UMAP_corrected_2", colour_by = "post_hoc_sample_gate", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_gate_colours, name = "post_hoc_sample_gate")
p16 <- plotReducedDim(sce0, "UMAP_corrected_3", colour_by = "post_hoc_sample_gate", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = sample_gate_colours, name = "post_hoc_sample_gate")
p17 <- plotReducedDim(sce0, "UMAP", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p18 <- plotReducedDim(sce0, "UMAP_corrected_1", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p19 <- plotReducedDim(sce0, "UMAP_corrected_2", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p20 <- plotReducedDim(sce0, "UMAP_corrected_3", colour_by = "tissue", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = tissue_colours, name = "tissue")
p21 <- plotReducedDim(sce0, "UMAP", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p22 <- plotReducedDim(sce0, "UMAP_corrected_1", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p23 <- plotReducedDim(sce0, "UMAP_corrected_2", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p24 <- plotReducedDim(sce0, "UMAP_corrected_3", colour_by = "donor", theme_size = 7, point_size = 0.2) +
  scale_colour_manual(values = donor_colours, name = "donor")
p1 + p2 + p3 + p4 +
  p5 + p6 + p7 + p8 +
  p9 + p10 + p11 + p12 +
  p13 + p14 + p15 + p16 +
  p17 + p18 + p19 + p20 +
  p21 + p22 + p23 + p24 +
  plot_layout(ncol = 4, guides = "collect")
```

To get an insight, here we breakdown the UMAP plot above further by `plate_number` (Figure \@ref(fig:umap-plate-number)), `post_hoc_sample_gate` (Figure \@ref(fig:umap-post-hoc-sample-gate)), and `tissue` (Figure \@ref(fig:umap-tissue)).

From the perspective of correcting the batch effect, the MNN correction can effectively alleviate the plate-specific grouping of cells (Figure \@ref(fig:umap-plate-number)).

From the perspective of preserving the biological-feature-of-interest, after the batch correction, we can still observe: (i) the "shifting" of the transcriptomic stages of cells P6 to P7, then P8 (Figure \@ref(fig:umap-post-hoc-sample-gate)); (ii) Blood contains cells mostly from P8 only, whilst thymus contains cells from all 3 phases, i.e. P6, P7 and P8 (Figure \@ref(fig:umap-tissue)).

```{r umap-plate-number, fig.cap = "UMAP plot of the dataset. Each point represents a cell and each panel highlights cells from a particular `plate_number` when data is unmerged (left) and merged by manual merge 2 (right).", fig.asp = 3/4, layout="l-body"}
umap_df <- makePerCellDF(sce0)
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
ggplot(aes(x = UMAP_corrected_3.1, y = UMAP_corrected_3.2), data = umap_df) +
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
```

```{r umap-post-hoc-sample-gate, fig.cap = "UMAP plot of the dataset. Each point represents a cell and each panel highlights cells from a particular `post_hoc_sample_gate` when data is unmerged (left) and merged by manual merge 2 (right).", fig.asp = 3/4, layout="l-body"}
umap_df <- makePerCellDF(sce0)
bg <- dplyr::select(umap_df, -post_hoc_sample_gate)
plot_grid(
ggplot(aes(x = UMAP.1, y = UMAP.2), data = umap_df) +
  geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
  geom_point(aes(colour = post_hoc_sample_gate), alpha = 1, size = 0.5) +
  scale_fill_manual(
    values = sample_gate_colours,
    name = "post_hoc_sample_gate") +
  scale_colour_manual(
    values = sample_gate_colours,
    name = "post_hoc_sample_gate") +
  theme_cowplot(font_size = 10) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  facet_wrap(~post_hoc_sample_gate, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  guides(colour = FALSE),
ggplot(aes(x = UMAP_corrected_3.1, y = UMAP_corrected_3.2), data = umap_df) +
  geom_point(data = bg, colour = scales::alpha("grey", 0.5), size = 0.125) +
  geom_point(aes(colour = post_hoc_sample_gate), alpha = 1, size = 0.5) +
  scale_fill_manual(
    values = sample_gate_colours,
    name = "post_hoc_sample_gate") +
  scale_colour_manual(
    values = sample_gate_colours,
    name = "post_hoc_sample_gate") +
  theme_cowplot(font_size = 10) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  facet_wrap(~post_hoc_sample_gate, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  guides(colour = FALSE),
ncol= 2,
align ="h"
)
```

```{r umap-tissue, fig.cap = "UMAP plot of the dataset. Each point represents a cell and each panel highlights cells from a particular `tissue`when data is unmerged (left) and merged by manual merge 2 (right).", fig.asp = 3/4, layout="l-body"}
umap_df <- makePerCellDF(sce0)
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
ggplot(aes(x = UMAP_corrected_3.1, y = UMAP_corrected_3.2), data = umap_df) +
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
```

Further, we also investigate the usefulness of the control cell line included for batch correction.
The cell line alone, unfortunately, seems not be able to help with the batch correction, but either exacerbated or leave the batch impact unresolved in this scenario (Figure \@ref(fig:umap-plate-number-cell-line-guided-MNN)).

```{r, results = "hide"}
sce1 <- sce
sce1$batch <- sce1$plate_number
var_fit <- modelGeneVarWithSpikes(sce1, "ERCC", block = sce1$batch)
hvg <- getTopHVGs(var_fit, var.threshold = 0)
hvg <- setdiff(hvg, c(ribo_set, mito_set, pseudogene_set))
library(batchelor)
set.seed(1819)
# auto merge
mnn_out_4 <- fastMNN(
  multiBatchNorm(sce1, batch = sce1$batch),
  batch = sce1$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce1, "PCA")),
  auto.merge = TRUE,
  restrict = list(sce1$sample_name == "Cell line"),
  subset.row = hvg)
# manual merge 1
# NOTE: order defined based on the similarities in cell number between plates for each`sample_name`
mnn_out_5 <- fastMNN(
  multiBatchNorm(sce1, batch = sce1$batch),
  batch = sce1$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce1, "PCA")),
  auto.merge = FALSE,
  merge.order = list(list("LCE508", "LCE511", "LCE512", "LCE513", "LCE514"), "LCE509"),
  restrict = list(sce1$sample_name == "Cell line"),
  subset.row = hvg)
# manual merge 2
# NOTE: to reduced var loss of LCE509, which contained different proportion of cells from each sample, I decide to merge LCE509 after the merge of LCE503 + LCE504 (as indicated in the auto-merge)
mnn_out_6 <- fastMNN(
  multiBatchNorm(sce1, batch = sce1$batch),
  batch = sce1$batch,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce1, "PCA")),
  auto.merge = FALSE,
  merge.order = list(list("LCE513", "LCE514", "LCE509"), list("LCE508", "LCE511", "LCE512")),
  restrict = list(sce1$sample_name == "Cell line"),
  subset.row = hvg)
reducedDim(sce1, "corrected_1") <- reducedDim(mnn_out_4, "corrected")
reducedDim(sce1, "corrected_2") <- reducedDim(mnn_out_5, "corrected")
reducedDim(sce1, "corrected_3") <- reducedDim(mnn_out_6, "corrected")
# generate UMAP
set.seed(18191)
sce1 <- runUMAP(sce1, dimred = "PCA", name = "UMAP")
uncorrected_umap <- cbind(
  data.frame(
    x = reducedDim(sce1, "UMAP")[, 1],
    y = reducedDim(sce1, "UMAP")[, 2]),
  as.data.frame(colData(sce1)))
set.seed(11901)
sce1 <- runUMAP(sce1, dimred = "corrected_1", name = "UMAP_corrected_1")
corrected_umap_1 <- cbind(
  data.frame(
    x = reducedDim(sce1, "UMAP_corrected_1")[, 1],
    y = reducedDim(sce1, "UMAP_corrected_1")[, 2]),
  as.data.frame(colData(sce1)))
set.seed(11901)
sce1 <- runUMAP(sce1, dimred = "corrected_2", name = "UMAP_corrected_2")
corrected_umap_2 <- cbind(
  data.frame(
    x = reducedDim(sce1, "UMAP_corrected_2")[, 1],
    y = reducedDim(sce1, "UMAP_corrected_2")[, 2]),
  as.data.frame(colData(sce1)))
set.seed(11901)
sce1 <- runUMAP(sce1, dimred = "corrected_3", name = "UMAP_corrected_3")
corrected_umap_3 <- cbind(
  data.frame(
    x = reducedDim(sce1, "UMAP_corrected_3")[, 1],
    y = reducedDim(sce1, "UMAP_corrected_3")[, 2]),
  as.data.frame(colData(sce1)))
# re-clustering after each MNN correction
set.seed(4759)
snn_gr_1 <- buildSNNGraph(sce1, use.dimred = "corrected_1")
clusters_1 <- igraph::cluster_louvain(snn_gr_1)
sce1$cluster_1 <- factor(clusters_1$membership)
cluster_colours_1 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce1$cluster_1))],
  levels(sce1$cluster_1))
sce1$colours$cluster_colours_1 <- cluster_colours_1[sce1$cluster_1]
set.seed(4759)
snn_gr_2 <- buildSNNGraph(sce1, use.dimred = "corrected_2")
clusters_2 <- igraph::cluster_louvain(snn_gr_2)
sce1$cluster_2 <- factor(clusters_2$membership)
cluster_colours_2 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce1$cluster_2))],
  levels(sce1$cluster_2))
sce1$colours$cluster_colours_2 <- cluster_colours_2[sce1$cluster_2]
set.seed(4759)
snn_gr_3 <- buildSNNGraph(sce1, use.dimred = "corrected_3")
clusters_3 <- igraph::cluster_louvain(snn_gr_3)
sce1$cluster_3 <- factor(clusters_3$membership)
cluster_colours_3 <- setNames(
  scater:::.get_palette("tableau10medium")[seq_len(nlevels(sce1$cluster_3))],
  levels(sce1$cluster_3))
sce1$colours$cluster_colours_3 <- cluster_colours_3[sce1$cluster_3]
```

```{r umap-plate-number-cell-line-guided-MNN, fig.cap = "UMAP plot of the dataset. Each point represents a cell and each panel highlights cells from a particular `plate_number` when data is unmerged (left) and merged by manual merge 2 (right).", fig.asp = 3/4, layout="l-body"}
umap_df <- makePerCellDF(sce1)
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
ggplot(aes(x = UMAP_corrected_3.1, y = UMAP_corrected_3.2), data = umap_df) +
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
```






