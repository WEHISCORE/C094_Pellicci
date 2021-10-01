# Start vs. End test

# Using tradeSeq
# # set.seed(666)
# # i <- sort(sample(nrow(sce), 1000))
# i <- head(rownames(res), 1000)
# system.time(svet <- startVsEndTest(fit[i, ]))
# system.time(svet <- startVsEndTest(fit))
# svet <- svet[order(svet$pvalue), ]
# svet$FDR <- p.adjust(svet$pvalue, method = "BH")
# NOTE: Using the signed Wald statistic as the statistic, following by
#       analogy what Gordon recommends
#       (https://support.bioconductor.org/p/112937/#112938).
# svet$statistic <- svet$waldStat * sign(svet$logFClineage1)
# svet$ENTREZID <- sapply(
#   rowData(sce[rownames(svet), ])$NCBI.ENTREZID,
#   "[[",
#   1)
# svet$logFC <- svet$logFClineage1

# Using TSCAN
pseudo <- testPseudotime(
  # TODO: What sort of gene filtering?
  sce[rowSums(counts(sce) > 0) >= 10, ],
  pseudotime = sce$slingPseudotime_1,
  block = sce$plate_number)
svet <- pseudo[order(pseudo$p.value), ]
# TODO: Anything better? Perhaps convert P-value to signed z-value?
# svet$statistic <- svet$logFC
# NOTE: Convert 2-sided P-value to a signed z-score
svet$statistic <- qnorm(svet$p.value / 2, lower.tail = FALSE) * sign(svet$logFC)
svet$ENTREZID <- sapply(
  rowData(sce[rownames(svet), ])$NCBI.ENTREZID,
  "[[",
  1)

up_left <- rownames(svet[order(svet$statistic, decreasing = FALSE), ])[1]
# plotSmoothers(fit, counts(sce), gene = up_left)
plotExpression(sce, up_left, x = "slingPseudotime_1", show_smooth = TRUE)

up_right <- rownames(svet[order(svet$statistic, decreasing = TRUE), ])[1]
# plotSmoothers(fit, counts(sce), gene = up_right)
plotExpression(sce, up_right, x = "slingPseudotime_1", show_smooth = TRUE)

plotReducedDim(
  sce,
  dimred = "UMAP_corrected",
  colour_by = up_left,
  point_alpha = 1,
  point_size = 0.5) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2)

plotReducedDim(
  sce,
  dimred = "UMAP_corrected",
  colour_by = up_right,
  point_alpha = 1,
  point_size = 0.5) +
  geom_path(data = embedded, aes(x = Dim.1, y = Dim.2), size = 1.2)

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("stage", "tissue"),
  features = head(setdiff(rownames(svet[svet$FDR < 0.05, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 6)

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("stage", "tissue"),
  features = head(setdiff(rownames(svet[svet$FDR < 0.05 & svet$logFC < 0, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 6)

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("stage", "tissue"),
  features = head(setdiff(rownames(svet[svet$FDR < 0.05 & svet$logFC > 0, ]), mito_set), 50),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 6)

pseudotime_degs <- list(
  up_left =
    sapply(
      rowData(
        sce[rownames(svet[svet$FDR < 0.05 & svet$logFC < 0, ]), ])$NCBI.ENTREZID,
      "[[",
      1),
  up_right = sapply(
    rowData(
      sce[rownames(svet[svet$FDR < 0.05 & svet$logFC > 0, ]), ])$NCBI.ENTREZID,
    "[[",
    1))
pseudotime_degs <- sapply(pseudotime_degs, function(x) x[!is.na(x)])

go <- goana(de = pseudotime_degs, species = "Hs")
topGO(go, sort = "up_left")
topGO(go, sort = "up_right")

kegg <- kegga(de = pseudotime_degs, species = "Hs")
topKEGG(kegg, sort = "up_left")
topKEGG(kegg, sort = "up_right")

# NOTE: Using BiocFileCache to avoid re-downloading these gene sets everytime
#       the report is rendered.
library(BiocFileCache)
bfc <- BiocFileCache()
# NOTE: Creating list of gene sets in this slightly convoluted way so that each
#       gene set name is prepended by its origin (e.g. H, C2, or C7).
msigdb <- do.call(
  c,
  list(
    H = readRDS(
      bfcrpath(
        bfc,
        "http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds")),
    C2 = readRDS(
      bfcrpath(
        bfc,
        "http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.all.v7.1.entrez.rds")),
    C7 = readRDS(
      bfcrpath(
        bfc,
        "http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c7.all.v7.1.entrez.rds"))))

y_idx <- ids2indices(msigdb, id = svet$ENTREZID)
cam <- cameraPR(
  statistic = setNames(svet$statistic, svet$ENTREZID),
  index = y_idx,
  # TODO: use.ranks if don't trust parametric assumptions.
  use.ranks = FALSE)
head(cam)

term <- "C2.REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"
library(org.Hs.eg.db)
key <- msigdb[[term]]
tmp <- AnnotationDbi::select(
  org.Hs.eg.db,
  key = key,
  columns = "SYMBOL",
  keytype = "ENTREZID")
# y_index <- ids2indices(unique(tmp$SYMBOL), rownames(svet))
y_index <- ids2indices(msigdb[[term]], svet$ENTREZID)

plotHeatmap(
  sce,
  order_columns_by = "slingPseudotime_1",
  colour_columns_by = c("stage", "tissue"),
  features = intersect(tmp$SYMBOL, rownames(svet)),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-3, 3),
  color = hcl.colors(101, "Blue-Red 3"),
  fontsize = 6)

barcodeplot(statistics = svet$statistic, index = y_index[[1]])
