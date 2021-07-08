
#### `Thymus.S3` vs `Thymus.S2`


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

# Re-level to set the baseline/ref
sce$stage <- relevel(sce$stage, "S1")
sce$tissue <- relevel(sce$tissue, "Thymus")
sce$donor <- relevel(sce$donor, "1")

# define group and replicates
sce$group <- paste0(sce$tissue, ".", sce$stage)
sce$rep <- paste0(sce$group, ".", sce$donor, ".", sce$plate_number)




# checkpoint
tmp <- sce
# focus on subset
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S2") |
             (tmp$tissue == "Thymus" & tmp$stage == "S3")]
# abundance
table(tmp$tissue, tmp$stage)




# aggregate replicates
summed <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "plate_number", "rep")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
# logNormCounts
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)
# vvv not sure
colnames(summed) <- summed$rep
# colLabels(summed) <- summed$group




# # MDS (difference in `tissue`)
# library(edgeR)
# plotMDS(summed, col = as.integer(factor(summed$tissue)))
# legend("topleft", legend = levels(factor(summed$tissue)), col = 1:nlevels(factor(summed$tissue)), pch = 16)
# # MDS (difference in `stage`)
# plotMDS(summed, col = as.integer(factor(summed$stage)))
# legend("topleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)



# focus on large enough aggregates (>10 cells)
summed_filt <- summed[, summed$ncells >= 10]








de_results <- pseudoBulkDGE(
  tmp,
  label = summed_filt$group,
  design = ~plate_number + group + donor,
  coef = "groupThymus.S3",
  condition = summed_filt$group)

metadata(de_results)

