
#### `Thymus.S2` vs `Blood.S3`


# read in SCE
sce <- readRDS(here("data", "SCEs", "C094_Pellicci.annotate.SCE.rds"))

# remove stage `Unknown` (as Dan suggested) to facilitate DE analysis
sce <- sce[,sce$stage !="Unknown"]
colData(sce) <- droplevels(colData(sce))

# simplify `stage` labels
abbre <- factor(
  dplyr::case_when(
    sce$stage == "S1 (CD4+/CD161-)" ~ "S1",
    sce$stage == "S2 (CD4-/CD161-)" ~ "S2",
    sce$stage == "S3 (CD4-/CD161+)" ~ "S3"))
sce$stage <- abbre

# Re-level to set the baseline.
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



# aggregate replicates
tmp <- aggregateAcrossCells(
  tmp,
  id = colData(tmp)[, c("group", "tissue", "stage", "donor", "plate_number", "rep")],
  coldata_merge = FALSE,
  use.dimred = FALSE,
  use.altexps = FALSE)
colnames(tmp) <- paste0(tmp$rep)

# logNormCounts
sizeFactors(tmp) <- NULL
tmp <- logNormCounts(tmp)
colLabels(tmp) <- tmp$group
colnames(tmp) <- tmp$rep




# # MDS (difference in `tissue`)
# library(edgeR)
# plotMDS(tmp, col = as.integer(factor(tmp$tissue)))
# legend("topleft", legend = levels(factor(tmp$tissue)), col = 1:nlevels(factor(tmp$tissue)), pch = 16)
# # MDS (difference in `stage`)
# plotMDS(tmp, col = as.integer(factor(tmp$stage)))
# legend("topleft", legend = levels(factor(tmp$stage)), col = 1:nlevels(factor(tmp$stage)), pch = 16)



# focus on large enough aggregates (>10 cells)
tmp <- tmp[, tmp$ncells >= 10]





group <- c("Thymus.S2", "Thymus.S3")



de_results <- pseudoBulkDGE(
  tmp,
  label = tmp$group,
  design = ~plate_number + group + donor,
  coef = group[1],
  condition = tmp$group)

metadata(de_results)

