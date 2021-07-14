
# NOTE: no fixed workflow for "multi-sample comparison"; highly customized and vary from project to project
# NOTE: for this project, DE analysis may not be necessary as ncell per aggregate could be too little
# client will have to determine DE from minibulk, then apply same set of DE to singlecell to confirm only

library(SingleCellExperiment)
library(here)
library(scater)
library(scran)
library(pheatmap)
# source(here("code", "helper_functions.R"))

# read in SCE
sce <- readRDS(here("data", "SCEs", "C094_Pellicci.annotate.SCE.rds"))

# remove stage `Unknown` (as Dan suggested) to facilitate DE analysis
sce <- sce[,sce$stage !="Unknown"]
colData(sce) <- droplevels(colData(sce))

# simplify `stage` labels
sce$stage <- factor(
  dplyr::case_when(
    sce$stage == "S1 (CD4+/CD161-)" ~ "S1",
    sce$stage == "S2 (CD4-/CD161-)" ~ "S2",
    sce$stage == "S3 (CD4-/CD161+)" ~ "S3"))

# define group and replicates
# NOTE: I: `group` = combination of  factor-of-interest
# NOTE: I: `group` should not be factorized, or error at `pseudoBulkDGE`
sce$group <- paste0(sce$tissue, ".", sce$stage)
# NOTE: I: `rep` = finest possible sample
# NOTE: as ncells per aggregate would be too little (<<100); have to skip `plate_number` when defining `rep` for this case only
sce$rep <- paste0(sce$group, ".", sce$donor)

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
# need to reset `stage` colour after simplification above
stage_colours <- setNames(
  palette.colors(nlevels(sce$stage), "Dark 2"),
  levels(sce$stage))
sce$colours$stage_colours <- stage_colours[sce$stage]
# define group colours (without factorizing `group`)
group_colours <- setNames(
  Polychrome::kelly.colors(nlevels(factor(sce$group))+1)[-1],
  levels(factor(sce$group)))
# cluster -> label for including "cluster_" only
cluster_colours <- setNames(
  unique(sce$colours$cluster_colours),
  unique(names(sce$colours$cluster_colours)))
cluster_colours <- cluster_colours[levels(sce$cluster)]
sce$label <- factor(
  paste0("cluster_", sce$cluster),
  levels = paste0("cluster_", 1:nlevels(sce$cluster)))
label_colours <- setNames(
  cluster_colours,
  paste0("cluster_", names(cluster_colours)))

# Re-level to set the baseline/ref
sce$stage <- relevel(sce$stage, "S1")
sce$tissue <- relevel(sce$tissue, "Thymus")
sce$donor <- relevel(sce$donor, "1")

# define cutoff
min_ncells <- 10
fdr <- 0.05


############################
# `Thymus.S2` vs `Thymus.S1`

# checkpoint
tmp <- sce
# focus on subset
# NOTE: can either subset before/after aggregate; but must subset, or aggregate out-of-comparison-scope cells
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S2" |
                tmp$tissue == "Thymus" & tmp$stage == "S1")]

# abundance
table(tmp$tissue, tmp$stage)

# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  # NOTE: line determine ncell per aggregate; be careful + include essential factor only
  # NOTE: to secure statistical power (optimal >100 cells per aggregate); check by: colData(summed)
  # NOTE: I: usu must include `cluster` here, or per aggregate contain multiple clusters - too much variance by itself already and limit DE
  # NOTE: also, if one ignore `cluster` then it will be DE of minibulk rather than singlecell
  # NOTE: here, `tissue` and `stage` is enveloped by `group`; `donor` is enveloped by `rep`
  # thus number of aggregate is the same for: c("group", "rep", "cluster"), c("group", "tissue", "stage", "donor", "rep", "label")
  # NOTE: `donor/rep` has to be included too or duplicated colnames will lead to error
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "rep", "label")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
# NOTE: colnames of `summed` must needed to be defined or no labels in MDS plot; usu: `replicate` (+ `cluster/label`)
# NOTE: colLabels could also be defined, but optional; usu: `cluster/label`
colData(summed)
colnames(summed) <- paste0(summed$rep, ".", summed$label)

# # Construct values for plotting in heatmaps (new)
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)
# Construct values for plotting in heatmaps (old)
assay(summed, "logCPM") <- edgeR::cpm(counts(summed), log = TRUE)
assay(summed, "correctedLogCPM") <- limma::removeBatchEffect(
  assay(summed, "logCPM"),
  batch = summed$plate_number,
  design = model.matrix(~group, colData(summed)))

# MDS (difference in `tissue`)
library(edgeR)
plotMDS(summed, col = as.integer(factor(summed$tissue)))
legend("topleft", legend = levels(factor(summed$tissue)), col = 1:nlevels(factor(summed$tissue)), pch = 16)
# MDS (difference in `stage`)
# NOTE: use MDS to trace for "outlier"; usu; be aggregate with too library ncells (and thus library size)
# NOTE: remember to colnames `summed` to help trace for outlier (by naming each sample)
plotMDS(summed, col = as.integer(factor(summed$stage)))
legend("topleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)

# focus on large enough aggregates
summed_filt <- summed[, summed$ncells >= min_ncells]

# DE
groups <- c("groupThymus.S2", "groupThymus.S1")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$label,
  # NOTE: should block on `plate_number` due to batch effect; but no choice or ncells per aggregate <<100
  # NOTE: "design = ~ group + donor"  means interest in variation between groups and difference between donors
  design = ~ group + donor,
  coef = groups[1],
  condition = summed_filt$group)

all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
# Need at least 2 features to make a heatmap.
stopifnot(length(features) < 2)

is_de <- decideTestsPerLabel(de_results, threshold = fdr)
summarizeTestsPerLabel(is_de) %>%
  knitr::kable(
    caption = "Number of DEGs per label. Each row corresponds to a label and each column corresponds to the number of downregulated genes (`-1`), the number of non-differentially expressed genes (`0`), the number of upregulated genes (`1`), and the number of genes not tested (`NA`).")



############################
# `Thymus.S3` vs `Thymus.S1`

# checkpoint
tmp <- sce
# focus on subset
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S3" |
                tmp$tissue == "Thymus" & tmp$stage == "S1")]

# abundance
table(tmp$tissue, tmp$stage)

# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "rep", "label")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
colData(summed)
colnames(summed) <- paste0(summed$rep, ".", summed$label)

# # Construct values for plotting in heatmaps (new)
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)
# Construct values for plotting in heatmaps (old)
assay(summed, "logCPM") <- edgeR::cpm(counts(summed), log = TRUE)
assay(summed, "correctedLogCPM") <- limma::removeBatchEffect(
  assay(summed, "logCPM"),
  batch = summed$plate_number,
  design = model.matrix(~group, colData(summed)))

# MDS (difference in `tissue`)
library(edgeR)
plotMDS(summed, col = as.integer(factor(summed$tissue)))
legend("topleft", legend = levels(factor(summed$tissue)), col = 1:nlevels(factor(summed$tissue)), pch = 16)
# MDS (difference in `stage`)
plotMDS(summed, col = as.integer(factor(summed$stage)))
legend("topleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)

# focus on large enough aggregates
summed_filt <- summed[, summed$ncells >= min_ncells]

# DE
groups <- c("groupThymus.S3", "groupThymus.S1")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$label,
  design = ~ group + donor,
  coef = groups[1],
  condition = summed_filt$group)

all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
# Need at least 2 features to make a heatmap.
stopifnot(length(features) < 2)

# heatmap
# NOTE: I: order by `group` then `cluster/label`
# then one can see if there are different cell type/status under each group
# e.g. Thymus.S3 has cluster.1,2,3, whilst Blood.S3 only got cluster.3
# where both cluster.3 should have similar expression/transcriptome
all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
plotHeatmap(
  summed_filt,
  features = features,
  columns = order(
    factor(summed_filt$group),
    summed_filt$label,
    summed_filt$tissue,
    summed_filt$stage,
    summed_filt$donor),
  colour_columns_by = c(
    "group",
    "label",
    "tissue",
    "stage",
    "donor"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  annotation_row = data.frame(
    cluster_4 = ifelse(features %in% all_features$cluster_4, "DE", "not DE"),
    cluster_3 = ifelse(features %in% all_features$cluster_3, "DE", "not DE"),
    cluster_2 = ifelse(features %in% all_features$cluster_2, "DE", "not DE"),
    cluster_1 = ifelse(features %in% all_features$cluster_1, "DE", "not DE"),
    row.names = features),
  column_annotation_colors = list(
    cluster_1 = c("not DE" = "white", "DE" = label_colours[["cluster_1"]]),
    cluster_2 = c("not DE" = "white", "DE" = label_colours[["cluster_2"]]),
    cluster_3 = c("not DE" = "white", "DE" = label_colours[["cluster_3"]]),
    cluster_4 = c("not DE" = "white", "DE" = label_colours[["cluster_4"]]),
    group = group_colours[levels(factor(summed_filt$group))],
    label = label_colours,
    tissue = tissue_colours,
    stage = stage_colours,
    donor = donor_colours),
  main = paste0(groups[1], "_vs_", groups[2]),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

is_de <- decideTestsPerLabel(de_results, threshold = fdr)
summarizeTestsPerLabel(is_de) %>%
  knitr::kable(
    caption = "Number of DEGs per label. Each row corresponds to a label and each column corresponds to the number of downregulated genes (`-1`), the number of non-differentially expressed genes (`0`), the number of upregulated genes (`1`), and the number of genes not tested (`NA`).")



############################
# `Thymus.S3` vs `Thymus.S2`

# checkpoint
tmp <- sce
# focus on subset
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S3" |
                tmp$tissue == "Thymus" & tmp$stage == "S2")]

# abundance
table(tmp$tissue, tmp$stage)

# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "rep", "label")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
colData(summed)
colnames(summed) <- paste0(summed$rep, ".", summed$label)

# # Construct values for plotting in heatmaps (new)
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)
# Construct values for plotting in heatmaps (old)
assay(summed, "logCPM") <- edgeR::cpm(counts(summed), log = TRUE)
assay(summed, "correctedLogCPM") <- limma::removeBatchEffect(
  assay(summed, "logCPM"),
  batch = summed$plate_number,
  design = model.matrix(~group, colData(summed)))

# MDS (difference in `tissue`)
library(edgeR)
plotMDS(summed, col = as.integer(factor(summed$tissue)))
legend("topleft", legend = levels(factor(summed$tissue)), col = 1:nlevels(factor(summed$tissue)), pch = 16)
# MDS (difference in `stage`)
plotMDS(summed, col = as.integer(factor(summed$stage)))
legend("topleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)

# focus on large enough aggregates
summed_filt <- summed[, summed$ncells >= min_ncells]

# DE
groups <- c("groupThymus.S3", "groupThymus.S2")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$label,
  design = ~ group + donor,
  coef = groups[1],
  condition = summed_filt$group)

all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
# Need at least 2 features to make a heatmap.
stopifnot(length(features) < 2)

# heatmap
all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
plotHeatmap(
  summed_filt,
  features = features,
  columns = order(
    factor(summed_filt$group),
    summed_filt$label,
    summed_filt$tissue,
    summed_filt$stage,
    summed_filt$donor),
  colour_columns_by = c(
    "group",
    "label",
    "tissue",
    "stage",
    "donor"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  annotation_row = data.frame(
    cluster_4 = ifelse(features %in% all_features$cluster_4, "DE", "not DE"),
    cluster_3 = ifelse(features %in% all_features$cluster_3, "DE", "not DE"),
    cluster_2 = ifelse(features %in% all_features$cluster_2, "DE", "not DE"),
    cluster_1 = ifelse(features %in% all_features$cluster_1, "DE", "not DE"),
    row.names = features),
  column_annotation_colors = list(
    cluster_1 = c("not DE" = "white", "DE" = label_colours[["cluster_1"]]),
    cluster_2 = c("not DE" = "white", "DE" = label_colours[["cluster_2"]]),
    cluster_3 = c("not DE" = "white", "DE" = label_colours[["cluster_3"]]),
    cluster_4 = c("not DE" = "white", "DE" = label_colours[["cluster_4"]]),
    group = group_colours[levels(factor(summed_filt$group))],
    label = label_colours,
    tissue = tissue_colours,
    stage = stage_colours,
    donor = donor_colours),
  main = paste0(groups[1], "_vs_", groups[2]),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

is_de <- decideTestsPerLabel(de_results, threshold = fdr)
summarizeTestsPerLabel(is_de) %>%
  knitr::kable(
    caption = "Number of DEGs per label. Each row corresponds to a label and each column corresponds to the number of downregulated genes (`-1`), the number of non-differentially expressed genes (`0`), the number of upregulated genes (`1`), and the number of genes not tested (`NA`).")



############################
# `Thymus.S3` vs `Blood.S3`

# checkpoint
tmp <- sce
# focus on subset
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S3" |
                tmp$tissue == "Blood" & tmp$stage == "S3")]

# abundance
table(tmp$tissue, tmp$stage)

# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "rep", "label")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
colData(summed)
colnames(summed) <- paste0(summed$rep, ".", summed$label)

# # Construct values for plotting in heatmaps (new)
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)
# Construct values for plotting in heatmaps (old)
assay(summed, "logCPM") <- edgeR::cpm(counts(summed), log = TRUE)
assay(summed, "correctedLogCPM") <- limma::removeBatchEffect(
  assay(summed, "logCPM"),
  batch = summed$plate_number,
  design = model.matrix(~group, colData(summed)))

# MDS (difference in `tissue`)
library(edgeR)
plotMDS(summed, col = as.integer(factor(summed$tissue)))
legend("topleft", legend = levels(factor(summed$tissue)), col = 1:nlevels(factor(summed$tissue)), pch = 16)
# MDS (difference in `stage`)
plotMDS(summed, col = as.integer(factor(summed$stage)))
legend("topleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)

# focus on large enough aggregates
summed_filt <- summed[, summed$ncells >= min_ncells]

# DE
groups <- c("groupThymus.S3", "groupBlood.S3")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$label,
  design = ~ group + donor,
  coef = groups[1],
  condition = summed_filt$group)

all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
# Need at least 2 features to make a heatmap.
stopifnot(length(features) < 2)

# heatmap
all_features <- lapply(
  de_results,
  function(x) rownames(x[which(x$FDR < fdr), ]))
features <- unique(unlist(all_features))
plotHeatmap(
  summed_filt,
  features = features,
  columns = order(
    factor(summed_filt$group),
    summed_filt$label,
    summed_filt$tissue,
    summed_filt$stage,
    summed_filt$donor),
  colour_columns_by = c(
    "group",
    "label",
    "tissue",
    "stage",
    "donor"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  annotation_row = data.frame(
    cluster_4 = ifelse(features %in% all_features$cluster_4, "DE", "not DE"),
    cluster_3 = ifelse(features %in% all_features$cluster_3, "DE", "not DE"),
    cluster_2 = ifelse(features %in% all_features$cluster_2, "DE", "not DE"),
    cluster_1 = ifelse(features %in% all_features$cluster_1, "DE", "not DE"),
    row.names = features),
  column_annotation_colors = list(
    cluster_1 = c("not DE" = "white", "DE" = label_colours[["cluster_1"]]),
    cluster_2 = c("not DE" = "white", "DE" = label_colours[["cluster_2"]]),
    cluster_3 = c("not DE" = "white", "DE" = label_colours[["cluster_3"]]),
    cluster_4 = c("not DE" = "white", "DE" = label_colours[["cluster_4"]]),
    group = group_colours[levels(factor(summed_filt$group))],
    label = label_colours,
    tissue = tissue_colours,
    stage = stage_colours,
    donor = donor_colours),
  main = paste0(groups[1], "_vs_", groups[2]),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 7)

is_de <- decideTestsPerLabel(de_results, threshold = fdr)
summarizeTestsPerLabel(is_de) %>%
  knitr::kable(
    caption = "Number of DEGs per label. Each row corresponds to a label and each column corresponds to the number of downregulated genes (`-1`), the number of non-differentially expressed genes (`0`), the number of upregulated genes (`1`), and the number of genes not tested (`NA`).")
