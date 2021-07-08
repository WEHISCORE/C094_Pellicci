
library(SingleCellExperiment)
library(here)
library(scater)
library(scran)



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




### DE between stages in thymus

# checkpoint + subset
tmp <- sce
tmp <- tmp[, tmp$tissue == "Thymus"]

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

# MDS plot
plotMDS(summed, col = as.integer(factor(summed$stage)))
legend("bottomleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)

# design
# TODO: need to block by `plate_number` (need to workaround: all ref levels hide if block)
design <- model.matrix(~0 + group + donor, colData(summed))

# TODO: (need to workaround: unmatch ./. nrow(design) and ncol(summed_filt))
# # focus on large enough aggregates (>10 cells)
# summed_filt <- summed[, summed$ncells >= 10]

# model fit
summed <- calcNormFactors(summed)
summed <- estimateDisp(summed, design)
fit <- glmQLFit(summed, design)

# define contrast
contr <- makeContrasts(
  Thymus.S3_vs_Thymus.S1 = groupThymus.S3 - groupThymus.S1,
  Thymus.S2_vs_Thymus.S1 = groupThymus.S2 - groupThymus.S1,
  Thymus.S3_vs_Thymus.S2 = groupThymus.S3 - groupThymus.S2,
  levels = design)

# sumup table
for (j in colnames(contr)) {
  qlf <- glmQLFTest(fit, contrast = contr[, j])
  print(summary(decideTests(qlf)))
  glMDPlot(
    qlf,
    counts = summed$counts,
    groups = summed$group,
    status = decideTests(qlf),
    transform = TRUE,
    path = here("output"),
    html = j,
    main = j,
    launch = FALSE)
}

# -1*groupThymus.S1 1*groupThymus.S3
# Down                                   48
# NotSig                              26993
# Up                                    958
# -1*groupThymus.S1 1*groupThymus.S2
# Down                                    1
# NotSig                              27981
# Up                                     17
# -1*groupThymus.S2 1*groupThymus.S3
# Down                                   31
# NotSig                              27659
# Up                                    309





### DE between blood (S3) and thymus (S3)

# checkpoint + subset
tmp <- sce
tmp <- tmp[, (tmp$tissue == "Thymus" & tmp$stage == "S3") |
             (tmp$tissue == "Blood" & tmp$stage == "S3")]
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

# MDS plots
plotMDS(summed, col = as.integer(factor(summed$tissue)))
legend("topleft", legend = levels(factor(summed$tissue)), col = 1:nlevels(factor(summed$tissue)), pch = 16)
plotMDS(summed, col = as.integer(factor(summed$stage)))
legend("topleft", legend = levels(factor(summed$stage)), col = 1:nlevels(factor(summed$stage)), pch = 16)

# design
# TODO: need to block by `plate_number` (need to workaround: all ref levels hide if block)
design <- model.matrix(~0 + group + donor, colData(summed))

# TODO: (need to workaround: unmatch ./. nrow(design) and ncol(summed_filt))
# # focus on large enough aggregates (>10 cells)
# summed_filt <- summed[, summed$ncells >= 10]

# model fit
summed <- calcNormFactors(summed)
summed <- estimateDisp(summed, design)
fit <- glmQLFit(summed, design)

# define contrast
contr <- makeContrasts(
  Thymus.S3_vs_Blood.S3 = groupThymus.S3 - groupBlood.S3,
  levels = design)

# sumup table
for (j in colnames(contr)) {
  qlf <- glmQLFTest(fit, contrast = contr[, j])
  print(summary(decideTests(qlf)))
  glMDPlot(
    qlf,
    counts = summed$counts,
    groups = summed$group,
    status = decideTests(qlf),
    transform = TRUE,
    path = here("output"),
    html = j,
    main = j,
    launch = FALSE)
}

# -1*groupBlood.S3 1*groupThymus.S3
# Down                                  14
# NotSig                             27845
# Up                                   140
