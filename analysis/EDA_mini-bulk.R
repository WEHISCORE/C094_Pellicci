# EDA of mini-bulk samples for Dan's grant application.
# Peter Hickey
# 2021-04-16

# Setup mini-bulk data ---------------------------------------------------------

library(SingleCellExperiment)
library(here)
sce <- readRDS("data/SCEs/C094_Pellicci.scPipe.SCE.rds")
sce <- sce[, sce$sample_type != "Single cell"]
sce <- sce[, sce$sample_name != "Cell line"]
sce$sample_name <- sub(" P[0-9]$", "", sce$sample_name)

# Gene annotations -------------------------------------------------------------

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

# Helper functions -------------------------------------------------------------

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


library(scater)
createClusterMarkerOutputs <- function(sce, outdir, markers, k, ...) {

  # Create CSVs
  message("Creating CSVs")
  for (label in names(markers)) {
    message("\t", label)
    gzout <- gzfile(
      description = file.path(
        outdir,
        sprintf("%s.csv.gz", label)),
      open = "wb")
    write.csv(
      as.data.frame(markers[[label]]),
      gzout,
      # NOTE: quote = TRUE needed because some fields contain commas.
      quote = TRUE,
      row.names = FALSE)
    close(gzout)
  }

  # Create PDFs
  message("Creating PDFs")
  for (label in names(markers)) {
    message("\t", label)
    label_markers <- markers[[label]]
    # Filter genes
    # NOTE: No FDR cut-off since we're using the underpowered t-test
    features <- rownames(label_markers)
    # TODO: Useful?
    features <- setdiff(features, c(ribo_set, mito_set))
    # Select top-k
    features <- head(features, k)
    if (length(features) >= 2) {

      plotHeatmap(
        sce,
        features,
        color = hcl.colors(101, "Blue-Red 3"),
        label_cols = TRUE,
        show_colnames = FALSE,
        order_columns_by = c("tissue", "sample_name", "sample_gate"),
        center = TRUE,
        symmetric = TRUE,
        zlim = c(-3, 3),
        fontsize = 6,
        main = label,
        # NOTE: Leave it to pheatmap() to create the PDFs because multipage
        #       PDFs containing pheatmap() outputs are often broken (I
        #       suspect it's something to do with not closing the PDF
        #       graphics device, but I wasted too much time trying to fix
        #       this).
        filename = file.path(
          outdir,
          sprintf("%s.cell-level_heatmap.pdf", label)),
        ...)
    }
  }
}

# Marker genes -----------------------------------------------------------------

library(scuttle)
library(scran)
sce <- logNormCounts(sce)
sce$tissue <- ifelse(grepl("Blood", sce$sample_name), "Blood", "Thymus")
tissue_markers <- findMarkers(
  sce,
  sce$tissue,
  direction = "up",
  pval.type = "all")
sample_name_markers <- findMarkers(
  sce,
  sce$sample_name,
  direction = "up",
  pval.type = "all")
blood_sample_name_markers <- findMarkers(
  sce[sce$tissue == "Blood"],
  sce$sample_name[sce$tissue == "Blood"],
  direction = "up",
  pval.type = "all")
thymus_sample_name_markers <- findMarkers(
  sce[sce$tissue == "Thymus"],
  sce$sample_name[sce$tissue == "Thymus"],
  direction = "up",
  pval.type = "all")

outdir <- here("tmp")
dir.create(file.path(outdir, "tissue"), recursive = TRUE)
createClusterMarkerOutputs(sce, file.path(outdir, "tissue"), tissue_markers, k = 50)
dir.create(file.path(outdir, "sample_name"), recursive = TRUE)
createClusterMarkerOutputs(sce, file.path(outdir, "sample_name"), sample_name_markers, k = 50)
dir.create(file.path(outdir, "blood"), recursive = TRUE)
createClusterMarkerOutputs(sce, file.path(outdir, "blood"), blood_sample_name_markers, k = 50)
dir.create(file.path(outdir, "thymus"), recursive = TRUE)
createClusterMarkerOutputs(sce, file.path(outdir, "thymus"), thymus_sample_name_markers, k = 50)
