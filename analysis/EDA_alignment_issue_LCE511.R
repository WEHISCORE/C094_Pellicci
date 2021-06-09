library(SingleCellExperiment)
sce <- readRDS("data/SCEs/C094_Pellicci.scPipe.SCE.rds")
sce$sum <- colSums(counts(sce))
sce$row <- substr(sce$well_position, 1, 1)
sce$col <- as.integer(substr(sce$well_position, 2, nchar(as.character(sce$well_position))))

library(scater)
z <- sce[, sce$plate_number == "LCE511"]
z <- z[, order(z$row, z$col)]
plotColData(z, "sum", x = "col", other_fields = "row", colour_by = I(sapply(strsplit(z$sample_name, " "), "[[", 1))) +
  scale_y_log10() +
  facet_wrap(~row, ncol = 4) +
  geom_vline(xintercept = c(3, 22), lty = 2, col = "blue")

z <- sce[, sce$plate_number == "LCE511" & substr(sce$well_position, 1, 1) %in% c("L", "M", "N", "O")]
z <- z[, order(z$row, z$col)]
plotColData(z, "sum", x = "col", other_fields = "row", colour_by = I(sapply(strsplit(z$sample_name, " "), "[[", 1))) +
  scale_y_log10() +
  facet_grid(~row) +
  geom_vline(xintercept = c(3, 22), lty = 2, col = "blue")

# Order wells and make a matrix of libsizes. then row normalize and heatmap
zz <- matrix(
  log10(z$sum),
  ncol = 24,
  nrow = 4,
  dimnames = list(unique(z$row), as.character(unique((z$col)))))
library(pheatmap)
pheatmap(zz, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(zz, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(zz, scale = "column", cluster_rows = FALSE, cluster_cols = FALSE)

zzz <- t(apply(zz, 1, function(x) x - mean(x)))
zzz1 <- colMedians(zzz[c("L", "M", "N"), ])
zzz2 <- zzz[c("O"), ]
plot(
  1:24,
  zzz2 - zzz1,
  ylim = c(-1, 1),
  xlab = "col",
  ylab = "log(O / median(L, M, N)",
  main = "Library size ratios")
abline(h = 0, lty = 2, col = "red")
abline(v = c(3, 22), lty = 2, col = "blue")

z1 <- z[, z$row %in% "N"]
z2 <- z[, z$row %in% "O"]

a <- tapply(z1$sum, z1$col, median)
b <- tapply(z2$sum, z2$col, max)
plot(
  z2$col,
  log10(a / b),
  ylim = c(-1, 1),
  xlab = "col",
  ylab = "log(O / N)",
  main = "Library size ratios")
abline(h = 0, lty = 2, col = "red")
abline(v = c(3, 22), lty = 2, col = "blue")
