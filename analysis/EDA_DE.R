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
# NOTE: `group` should not be factorized, or error at `pseudoBulkDGE`
sce$group <- paste0(sce$tissue, ".", sce$stage)
sce$rep <- paste0(sce$group, ".", sce$donor, ".", sce$plate_number)

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
# need to reset stage colour after simplification above
stage_colours <- setNames(
  palette.colors(nlevels(sce$stage), "Dark 2"),
  levels(sce$stage))
sce$colours$stage_colours <- stage_colours[sce$stage]
cluster_colours <- setNames(
  unique(sce$colours$cluster_colours),
  unique(sce$stage))
cluster_colours <- cluster_colours[levels(sce$cluster)]
# define group colours (without factorizing `group`)
group_colours <- setNames(
  Polychrome::kelly.colors(nlevels(factor(sce$group))+1)[-1],
  levels(factor(sce$group)))

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
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S2" |
                tmp$tissue == "Thymus" & tmp$stage == "S1")]

# abundance
table(tmp$tissue, tmp$stage)

# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "plate_number", "rep")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
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

groups <- c("groupThymus.S2", "groupThymus.S1")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$group,
  design = ~plate_number + group + donor,
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
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "plate_number", "rep")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
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


groups <- c("groupThymus.S3", "groupThymus.S1")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$group,
  design = ~plate_number + group + donor,
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
    summed_filt$plate_number,
    summed_filt$tissue,
    summed_filt$stage,
    summed_filt$donor),
  colour_columns_by = c(
    "group",
    "plate_number",
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
    Thymus.S1 = ifelse(features %in% all_features$Thymus.S1, "DE", "not DE"),
    Thymus.S3 = ifelse(features %in% all_features$Thymus.S3, "DE", "not DE"),
    row.names = features),
  main = paste0(names(all_features)[1], "_vs_", names(all_features)[2]),
  column_annotation_colors = list(
    Thymus.S1 = c("not DE" = "white", "DE" = group_colours[["Thymus.S1"]]),
    Thymus.S3 = c("not DE" = "white", "DE" = group_colours[["Thymus.S3"]]),
    group = group_colours[levels(factor(summed_filt$group))],
    plate_number = plate_number_colours,
    tissue = tissue_colours,
    stage = stage_colours,
    donor = donor_colours),
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
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "plate_number", "rep")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
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


groups <- c("groupThymus.S3", "groupThymus.S2")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$group,
  design = ~plate_number + group + donor,
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
# `Thymus.S3` vs `Blood.S3`

# checkpoint
tmp <- sce
# focus on subset
tmp <- tmp[, (tmp$tissue == "Blood" & tmp$stage == "S3" |
                tmp$tissue == "Thymus" & tmp$stage == "S3")]

# abundance
table(tmp$tissue, tmp$stage)

# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "plate_number", "rep")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
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


groups <- c("groupThymus.S3", "groupBlood.S3")
de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$group,
  design = ~plate_number + group + donor,
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











