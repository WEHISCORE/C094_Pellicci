# Prepare C094_Pellicci data for GEO submission
# Peter Hickey
# 2021-09-23

library(here)

outdir <- here("GEO")
dir.create(outdir, recursive = TRUE)
dir.create(file.path(outdir, "FASTQ"))
dir.create(file.path(outdir, "SCE"))

# FASTQs -----------------------------------------------------------------------

# NOTE: These RPI-level FASTQ files are created by code/scPipe.R
for (rpi in c("1", "2", "3", "10", "12", "15", "16")) {
  message(rpi)
  rpi <- paste0("RPI-", rpi)
  file.copy(
    from = here("extdata/scPipe/NN215", rpi, paste0(rpi, ".R1.fastq.gz")),
    to = file.path(outdir, "FASTQ"),
    recursive = FALSE,
    overwrite = FALSE)
  file.copy(
    from = here("extdata/scPipe/NN215", rpi, paste0(rpi, ".R2.fastq.gz")),
    to = file.path(outdir, "FASTQ"),
    recursive = FALSE,
    overwrite = FALSE)
}
for (rpi in c("5", "9")) {
  message(rpi)
  rpi <- paste0("RPI-", rpi)
  file.copy(
    from = here("extdata/scPipe/NN227", rpi, paste0(rpi, ".R1.fastq.gz")),
    to = file.path(outdir, "FASTQ"),
    recursive = FALSE,
    overwrite = FALSE)
  file.copy(
    from = here("extdata/scPipe/NN227", rpi, paste0(rpi, ".R2.fastq.gz")),
    to = file.path(outdir, "FASTQ"),
    recursive = FALSE,
    overwrite = FALSE)
}

# scRNA-seq SCE ----------------------------------------------------------------

# NOTE: Update SCE similar to that done in
#       analysis/C094_Pellicci.single-cell.preprocess.Rmd, but then undo some
#       of it (specifically, revert changes to the rowData and some cosmetic
#       changes to colData).
sce <- readRDS(here("data", "SCEs", "C094_Pellicci.scPipe.SCE.rds"))
source(here("analysis", "C094_Pellicci.preprocess.R"))
sce <- sce[, sce$sample_type == "Single cell"]
colData(sce) <- droplevels(colData(sce))
rownames(sce) <- rowData(sce)$ENSEMBL.GENEID
rowData(sce) <- S4Vectors::make_zero_col_DFrame(nrow(sce))
sce$colours <- NULL

# Gene counts
write.csv(
  x = as.data.frame(as.matrix(counts(sce))),
  file = gzfile(file.path(outdir, "SCE", "scRNA-seq.gene_counts.csv.gz")),
  row.names = TRUE)

# ERCC counts
write.csv(
  x = as.data.frame(as.matrix(counts(altExp(sce, "ERCC")))),
  file = gzfile(file.path(outdir, "SCE", "scRNA-seq.ERCC_counts.csv.gz")),
  row.names = TRUE)

# FACS data
write.csv(
  x = as.data.frame(as.matrix(assay(altExp(sce, "FACS"), "raw"))),
  file = gzfile(file.path(outdir, "SCE", "scRNA-seq.FACS.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(file.path(outdir, "SCE", "scRNA-seq.sample_sheet.csv.gz")),
  row.names = TRUE)

# Mini-bulk SCE ----------------------------------------------------------------

# NOTE: Update SCE similar to that done in
#       analysis/C094_Pellicci.mini-bulk.preprocess.Rmd, but then undo some of
#       it (specifically, revert changes to the rowData and some cosmetic
#       changes to colData).
sce_deduped <- readRDS(here("data", "SCEs", "C094_Pellicci.scPipe.SCE.rds"))
sce_not_deduped <- readRDS(
  here("data", "SCEs", "C094_Pellicci.not_UMI_deduped.scPipe.SCE.rds"))
colnames(sce_not_deduped) <- sub(
  "\\.not_UMI_deduped",
  "",
  colnames(sce_not_deduped))
stopifnot(
  identical(rownames(sce_deduped), rownames(sce_not_deduped)),
  identical(colnames(sce_deduped), colnames(sce_not_deduped)))
# Combine UMI and read counts a single SCE.
sce <- sce_deduped
assay(sce, "UMI_counts") <- assay(sce_deduped, "counts")
assay(sce, "read_counts") <- assay(sce_not_deduped, "counts")
# NOTE: Have to reorder ERCCs
assay(altExp(sce, "ERCC"), "UMI_counts") <- assay(
  altExp(sce_deduped, "ERCC"), "counts")[rownames(altExp(sce, "ERCC")), ]
assay(altExp(sce, "ERCC"), "read_counts") <- assay(
  altExp(sce_not_deduped, "ERCC"), "counts")[rownames(altExp(sce, "ERCC")), ]
# Nullify some now-unrequired data.
assay(sce, "counts") <- NULL
assay(altExp(sce, "ERCC"), "counts") <- NULL
sce$UMI_deduped <- NULL
source(here("analysis", "C094_Pellicci.preprocess.R"))
sce <- sce[, sce$sample_type != "Single cell"]
colData(sce) <- droplevels(colData(sce))
rownames(sce) <- rowData(sce)$ENSEMBL.GENEID
rowData(sce) <- S4Vectors::make_zero_col_DFrame(nrow(sce))
sce$colours <- NULL

# Gene counts
# NOTE: Exporting read counts.
write.csv(
  x = as.data.frame(as.matrix(assay(sce, "read_counts"))),
  file = gzfile(file.path(outdir, "SCE", "mini-bulk.gene_counts.csv.gz")),
  row.names = TRUE)

# ERCC counts
# NOTE: Exporting read counts.
write.csv(
  x = as.data.frame(as.matrix(assay(altExp(sce, "ERCC"), "read_counts"))),
  file = gzfile(file.path(outdir, "SCE", "mini-bulk.ERCC_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(file.path(outdir, "SCE", "mini-bulk.sample_sheet.csv.gz")),
  row.names = TRUE)
