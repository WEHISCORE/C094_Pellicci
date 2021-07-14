# Helper function to coerce a DataFrame to a data.frame while preserving column
# names as is.
.adf <- function(x) {
  setNames(as.data.frame(x), colnames(x))
}

# Collapse labels only found in fewer than `cutoff` proportion of cells in all
# patients as 'other'
.collapseLabel <- function(labels, patient, cutoff = 0.01) {
  tmp <- table(labels, patient)
  tmp2 <- apply(tmp, 2, function(x) (x / sum(x)) > cutoff)
  tmp3 <- rownames(tmp2)[rowAnys(tmp2)]
  tmp4 <- ifelse(
    labels %in% tmp3,
    as.character(labels),
    "other")
  factor(tmp4, names(sort(table(tmp4), decreasing = TRUE)))
}

# Helper function to Combine data from 2 SCEs using gene names.
# NOTE: This assumes more than I'd like about the rowData and doesn't do much
#       checking of these assumptions.
.combine <- function(x, y, rowData_by = c("ENSEMBL", "SYMBOL", "CHR")) {
  if (is.null(rowData_by)) {
    rowData <- dplyr::full_join(
      as.data.frame(rowData(x)) %>%
        tibble::rownames_to_column(var = "gene"),
      as.data.frame(rowData(y)) %>%
        tibble::rownames_to_column(var = "gene")) %>%
      tibble::column_to_rownames("gene") %>%
      DataFrame(., row.names = rownames(.))
  } else {
    rowData <- dplyr::full_join(
      as.data.frame(rowData(x)[, rowData_by, drop = FALSE]),
      as.data.frame(rowData(y)[, rowData_by, drop = FALSE]),
      by = rowData_by) %>%
      DataFrame(row.names = scater::uniquifyFeatureNames(
        .$ENSEMBL,
        .$SYMBOL))
    rownames(x) <- rownames(rowData)[match(rowData(x)$ENSEMBL, rowData$ENSEMBL)]
    rownames(y) <- rownames(rowData)[match(rowData(y)$ENSEMBL, rowData$ENSEMBL)]
  }

  colData <- rbind(colData(x), colData(y))

  counts <- matrix(
    data = 0L,
    nrow = nrow(rowData), ncol = nrow(colData),
    dimnames = list(rownames(rowData), rownames(colData)))
  counts[rownames(x), colnames(x)] <- counts(
    x,
    withDimnames = FALSE)
  counts[rownames(y), colnames(y)] <- counts(
    y,
    withDimnames = FALSE)

  stopifnot(
    identical(
      metadata(x)$scPipe$version,
      metadata(y)$scPipe$version))
  stopifnot(
    identical(
      metadata(x)$scPipe$QC_cols,
      metadata(y)$scPipe$QC_cols))
  stopifnot(
    identical(
      metadata(x)$scPipe$demultiplex_info$status,
      metadata(y)$scPipe$demultiplex_info$status))
  stopifnot(
    identical(
      metadata(x)$scPipe$UMI_dup_info$duplication.number,
      metadata(y)$scPipe$UMI_dup_info$duplication.number))
  stopifnot(identical(metadata(x)$Biomart, metadata(y)$Biomart))
  metadata <- list(
    scPipe = list(
      version = metadata(x)$scPipe$version,
      QC_cols = metadata(x)$scPipe$QC_cols,
      demultiplex_info = data.frame(
        status = metadata(x)$scPipe$demultiplex_info$status,
        count = metadata(x)$scPipe$demultiplex_info$count +
          metadata(y)$scPipe$demultiplex_info$count),
      UMI_dup_info = data.frame(
        duplication.number = metadata(
          x)$scPipe$UMI_dup_info$duplication.number,
        count = metadata(x)$scPipe$UMI_dup_info$count +
          metadata(y)$scPipe$UMI_dup_info$count)),
    Biomart = metadata(x)$Biomart)

  sce <- SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    assays = list(counts = counts),
    metadata = metadata)

  stopifnot(identical(int_metadata(x), int_metadata(y)))
  int_metadata(sce) <- int_metadata(x)

  # NOTE: Not trying to combine int_elementMetadata of objects. Each is a
  #       DataFrame with a zero-column DataFrame as a `rowPairs` column. This
  #       is effectively no data and the SCE constructor makes one, anyway.

  stopifnot(validObject(sce))
  sce
}

.cbindSCEs <- function(list_of_sce, rowData_by = 1:6) {
  do.call(
    cbind,
    lapply(list_of_sce, function(sce) {
      # NOTE: Some fudging to combine only the necessary bits of each SCE
      #       (basically, don't include any QC metrics).
      rowData(sce) <- rowData(sce)[, rowData_by]
      sce
    }))
}

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

# Plot SingleR scores on a reduced dimension plot from a SCE.
plotScoreReducedDim <- function(results, sce, dimred = "TSNE",
                                max.labels = 20, normalize = TRUE, ncol = 5,
                                ...) {
  scores <- results$scores
  rownames(scores) <- rownames(results)
  m <- rowMaxs(scale(t(scores)))
  to.keep <- head(order(m, decreasing = TRUE), max.labels)
  if (normalize) {
    mmax <- rowMaxs(scores)
    mmin <- rowMins(scores)
    scores <- (scores - mmin) / (mmax - mmin)
    scores <- scores ^ 3
  }
  scores <- scores[, to.keep, drop = FALSE]
  cns <- colnames(scores)
  p <- lapply(cns, function(cn) {
    scater::plotReducedDim(
      sce,
      dimred = dimred,
      colour_by = data.frame(Score = scores[, cn]),
      ...) +
      ggtitle(cn) +
      scale_fill_viridis_c(limits = force(if(normalize) c(0, 1) else NULL)) +
      guides(fill = guide_colourbar(title = "Score"))
  })
  cowplot::plot_grid(plotlist = p, ncol = ncol)
}

# NOTE: Need to use my own gene counting function because not using UMI
#       deduplication.
geneCountingNoUMIDedup <- function(outdir, bc_anno) {
  files <- list.files(file.path(outdir, "count"), full.names = TRUE)
  names(files) <- sub("\\.csv", "", basename(files))
  counts <- lapply(files, function(file) {
    message(basename(file))
    data.table::fread(file, select = 1)[, table(gene_id)]
  })
  genes <- Reduce(union, lapply(counts, names))
  x <- matrix(
    0L,
    nrow = length(genes),
    ncol = length(files),
    dimnames = list(genes, names(counts)))
  for (j in names(counts)) {
    xx <- counts[[j]]
    x[names(xx), j] <- xx
  }
  z <- cbind(
    data.frame(gene_id = rownames(x)),
    as.data.frame(x))
  data.table::fwrite(
    x = z,
    file = file.path(paste0(outdir, "_no_dedup"), "gene_count.csv"),
    row.names = FALSE,
    nThread = 1)
}
