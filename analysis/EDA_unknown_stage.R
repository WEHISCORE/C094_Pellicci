# EDA to look into cells of 'Unknown' stage.
# Peter Hickey
# 2021-08-04

library(here)
library(scater)

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.cells_selected.SCE.rds"))

# Hmm, removing 20% of cells if removing 'Unknown' cells.
table(sce$stage)
proportions(table(sce$stage))

# An association between donor and stage.
table(sce$donor, sce$stage)
round(proportions(table(sce$donor, sce$stage), 1), 2)
plotColData(sce, "donor", "stage")
chisq.test(table(sce$donor, sce$stage))

# An association between tissue and stage.
table(sce$tissue, sce$stage)
round(proportions(table(sce$tissue, sce$stage), 1), 2)
plotColData(sce, "tissue", "stage")
chisq.test(table(sce$tissue, sce$stage))

# Most, but not all, 'Unknown' cells cluster with S1 and S2 cells in cluster 3.
table(sce$cluster, sce$stage)
plotColData(sce, "cluster", "stage")

plotExpression(
  sce,
  features = c("V525_50_A_CD4_BV510", "B530_30_A_CD161_FITC"),
  x = "stage",
  exprs_values = "pseudolog")

plotExpression(
  sce,
  features = "V525_50_A_CD4_BV510",
  x = "cluster",
  other_fields = "stage",
  exprs_values = "pseudolog") +
  facet_wrap(~stage, ncol = 2)
plotExpression(
  sce,
  features = "B530_30_A_CD161_FITC",
  x = "cluster",
  other_fields = "stage",
  exprs_values = "pseudolog") +
  facet_wrap(~stage, ncol = 2)

plotReducedDim(sce, "UMAP_corrected", colour_by = "stage", point_alpha = 1)
