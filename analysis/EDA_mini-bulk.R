# EDA of mini-bulk samples for Dan's grant application.
# Peter Hickey
# 2021-04-16

#' # Setup
#' ## Load mini-bulk data

library(SingleCellExperiment)
library(here)
sce <- readRDS(here("data/SCEs/C094_Pellicci.scPipe.SCE.rds"))
sce <- sce[, sce$sample_type != "Single cell"]
sce <- sce[, sce$sample_name != "Cell line"]
sce$sample_name <- sub(" P[0-9]$", "", sce$sample_name)
sce$tissue <- sapply(strsplit(sce$sample_name, " "), "[[", 1)
sce$replicate <- sapply(strsplit(sce$sample_name, " "), "[[", 2)

#' ## Helper functions

# Take a DataFrame with AtomicList columns and return a DataFrame where these
# columns have been flattened by paste-ing together the elements separated by
# `sep`.
flattenDF <- function(x, sep = "; ") {
  DataFrame(
    endoapply(x, function(xx) {
      if (!is(xx, "AtomicList")) {
        return(xx)
      }
      unstrsplit(as(xx, "CharacterList"), sep = sep)
    }),
    row.names = rownames(x))
}

#' ## Add gene annotations

# Extract rownames (Ensembl IDs) to use as key in database lookups.
ensembl <- rownames(sce)

# Pull out useful gene-based annotations from the Ensembl-based database.
library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
ensdb <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]
# NOTE: These columns were customised for this project.
ensdb_columns <- c(
  "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
stopifnot(all(ensdb_columns %in% columns(ensdb)))
ensdb_df <- DataFrame(
  lapply(ensdb_columns, function(column) {
    mapIds(
      x = ensdb,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "GENEID",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)
# NOTE: Can't look up GENEID column with GENEID key, so have to add manually.
ensdb_df$ENSEMBL.GENEID <- ensembl
# NOTE: Homo.sapiens combines org.Hs.eg.db and
#       TxDb.Hsapiens.UCSC.hg38.knownGene (as well as others) and therefore
#       uses entrez gene and RefSeq based data.
library(Homo.sapiens)
# NOTE: These columns were customised for this project.
ncbi_columns <- c(
  # From TxDB: None required
  # From OrgDB
  "ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
stopifnot(all(ncbi_columns %in% columns(Homo.sapiens)))
ncbi_df <- DataFrame(
  lapply(ncbi_columns, function(column) {
    mapIds(
      x = Homo.sapiens,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "ENSEMBL",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)
rowData(sce) <- cbind(ensdb_df, ncbi_df)

# Replace the row names of the SCE by the gene symbols (where available).
library(scater)
rownames(sce) <- uniquifyFeatureNames(
  ID = rownames(sce),
  # NOTE: An Ensembl ID may map to 0, 1, 2, 3, ... gene symbols.
  #       When there are multiple matches only the 1st match is used.
  names = vapply(rowData(sce)$ENSEMBL.SYMBOL, function(x) {
    if (length(x)) {
      x[[1]]
    } else {
      NA_character_
    }
  },
  character(1)))

#' # MDS plots

library(edgeR)
x <- SE2DGEList(sce)

#' ## Individual technical replicates

#+ fig.asp = 1, fig.cap = "Coloured by tissue"
plotMDS(x, col = as.integer(factor(x$samples$tissue)))
legend("bottomleft", legend = levels(factor(x$samples$tissue)), col = 1:nlevels(factor(x$samples$tissue)), pch = 16)
#+ fig.asp = 1, fig.cap = "Coloured by sample_gate"
plotMDS(x, col = as.integer(factor(x$samples$sample_gate)))
legend("bottomleft", legend = levels(factor(x$samples$sample_gate)), col = 1:nlevels(factor(x$samples$sample_gate)), pch = 16)
#+ fig.asp = 1, fig.cap = "Coloured by replicate"
plotMDS(x, col = as.integer(factor(x$samples$replicate)))
legend("bottomleft", legend = levels(factor(x$samples$replicate)), col = 1:nlevels(factor(x$samples$replicate)), pch = 16)

#' ## Aggregated technical replicates

y <- sumTechReps(
  x,
  paste0(x$samples$tissue, ".", x$sample$sample_gate, ".", x$sample$replicate))
#+ fig.asp = 1, fig.cap = "Coloured by tissue"
plotMDS(y, col = as.integer(factor(y$samples$tissue)))
legend("topright", legend = levels(factor(x$samples$tissue)), col = 1:nlevels(factor(x$samples$tissue)), pch = 16)
#+ fig.asp = 1, fig.cap = "Coloured by sample_gate"
plotMDS(y, col = as.integer(factor(y$samples$sample_gate)))
legend("topright", legend = levels(factor(x$samples$sample_gate)), col = 1:nlevels(factor(x$samples$sample_gate)), pch = 16)
#+ fig.asp = 1, fig.cap = "Coloured by replicate"
plotMDS(y, col = as.integer(factor(y$samples$replicate)))
legend("topright", legend = levels(factor(x$samples$replicate)), col = 1:nlevels(factor(x$samples$replicate)), pch = 16)

#' # DE analyses

#' ## Dan's marker proteins
marker_proteins <- read.csv(
  here("data/marker_proteins/C094_Pellicci.marker_proteins.csv"))
knitr::kable(marker_proteins)

#' ## DE between stages in thymus

library(Glimma)

x <- SE2DGEList(sce)
x <- x[, x$samples$tissue == "Thymus"]
x$genes <- x$genes[, c("ENSEMBL.GENEID", "ENSEMBL.SYMBOL", "ENSEMBL.GENEBIOTYPE")]
x$samples$group <- x$samples$sample_gate
x <- sumTechReps(
  x,
  paste0(x$samples$tissue, ".", x$sample$sample_gate, ".", x$sample$replicate))

#+ fig.asp = 1, fig.cap = "Coloured by sample_gate"
plotMDS(x, col = as.integer(factor(x$samples$sample_gate)))
legend("bottomleft", legend = levels(factor(x$samples$sample_gate)), col = 1:nlevels(factor(x$samples$sample_gate)), pch = 16)

design <- model.matrix(~0 + sample_gate + replicate, x$samples)
# NOTE: Decreasing `min.count` in order to test a few more genes.
keep <- filterByExpr(x, min.count = 2)
x <- x[keep, , keep.lib.sizes = FALSE]
x <- calcNormFactors(x)
x <- estimateDisp(x, design)
fit <- glmQLFit(x, design)
contr <- makeContrasts(
  Thymus.P8_vs_Thymus.P6 = sample_gateP8 - sample_gateP6,
  Thymus.P7_vs_Thymus.P6 = sample_gateP7 - sample_gateP6,
  Thymus.P8_vs_Thymus.P7 = sample_gateP8 - sample_gateP7,
  levels = design)

idx <- ids2indices(list(marker_proteins = marker_proteins$Gene), rownames(x))

for (j in colnames(contr)) {
  qlf <- glmQLFTest(fit, contrast = contr[, j])
  print(summary(decideTests(qlf)))
  glMDPlot(
    qlf,
    counts = x$counts,
    groups = x$samples$group,
    status = decideTests(qlf),
    transform = TRUE,
    path = here("output"),
    html = j,
    main = j,
    launch = FALSE)
  print(camera(x, idx, design, contrast = contr[, j]))
  barcodeplot(qlf$table$logFC, index = idx$marker_proteins, main = j)
}

#' ## DE between blood (P5) and thymus (P8)

x <- SE2DGEList(sce)
x <- x[, (x$samples$tissue == "Thymus" & x$samples$sample_gate == "P8") |
         (x$samples$tissue == "Blood" & x$samples$sample_gate == "P5")]
x$genes <- x$genes[, c("ENSEMBL.GENEID", "ENSEMBL.SYMBOL", "ENSEMBL.GENEBIOTYPE")]
x$samples$group <- paste0(x$samples$tissue, ".", x$samples$sample_gate)
x <- sumTechReps(
  x,
  paste0(x$samples$tissue, ".", x$sample$sample_gate, ".", x$sample$replicate))

#+ fig.asp = 1, fig.cap = "Coloured by tissue"
plotMDS(x, col = as.integer(factor(x$samples$tissue)))
legend("topleft", legend = levels(factor(x$samples$tissue)), col = 1:nlevels(factor(x$samples$tissue)), pch = 16)
#+ fig.asp = 1, fig.cap = "Coloured by sample_gate"
plotMDS(x, col = as.integer(factor(x$samples$sample_gate)))
legend("topleft", legend = levels(factor(x$samples$sample_gate)), col = 1:nlevels(factor(x$samples$sample_gate)), pch = 16)

design <- model.matrix(~0 + group + replicate, x$samples)
# NOTE: Decreasing `min.count` in order to test a few more genes.
keep <- filterByExpr(x, min.count = 2)
x <- x[keep, , keep.lib.sizes = FALSE]
x <- calcNormFactors(x)
x <- estimateDisp(x, design)
fit <- glmQLFit(x, design)
contr <- makeContrasts(
  Thymus.P8_vs_Blood.P6 = groupThymus.P8 - groupBlood.P5,
  levels = design)

for (j in colnames(contr)) {
  qlf <- glmQLFTest(fit, contrast = contr[, j])
  print(summary(decideTests(qlf)))
  glMDPlot(
    qlf,
    counts = x$counts,
    groups = x$samples$group,
    status = decideTests(qlf),
    transform = TRUE,
    path = here("output"),
    html = j,
    main = j,
    launch = FALSE)
  print(camera(x, idx, design, contrast = contr[, j]))
  barcodeplot(qlf$table$logFC, index = idx$marker_proteins, main = j)
}
