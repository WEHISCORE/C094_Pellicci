# Process C118 (NN215 and NN227) with scPipe including trimming R2.
# Peter Hickey
# 2022-04-07

# Setup ------------------------------------------------------------------------

library(scPipe)
library(Rsubread)
library(here)
library(readxl)
library(dplyr)
library(janitor)

source(here("code", "helper_functions.R"))

options("mc.cores" = 1L)

# Construct NN215 sample sheet -------------------------------------------------

file_nn215 <- here(
  "data",
  "sample_sheets",
  "C094_DanPellici_SC_MB_NN215_Seqprimer19Mar21.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn215,
  sheet = "Samples and indexing",
  skip = 2,
  n_max = 1)

# NOTE: FACS data in columns >= "K"
facs_data_idx <- seq(which(LETTERS == "K"), ncol(header_row))
header_row <- c(
  paste0(colnames(header_row[, -facs_data_idx]), header_row[1, -facs_data_idx]),
  unlist(header_row[1, facs_data_idx], use.names = FALSE))
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn215 <- read_excel(
  path = file_nn215,
  sheet = "Samples and indexing",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)
# Tidy up names and empty rows/columns.
sample_sheet_nn215 <- bind_cols(
  clean_names(sample_sheet_nn215[, -facs_data_idx]),
  clean_names(sample_sheet_nn215[, facs_data_idx], case = "parsed"))
sample_sheet_nn215 <- remove_empty(
  sample_sheet_nn215,
  which = c("rows", "cols"))
sample_sheet_nn215 <- filter(
  sample_sheet_nn215,
  !is.na(plate_number),
  illumina_index_index_number_separate_index_read != "removed")
# Ensure FACS columns are stored as numeric (readxl sometimes fails, presumably
# to weird pattern of empty cells).
sample_sheet_nn215 <- sample_sheet_nn215 %>%
  mutate_at(facs_data_idx, as.numeric)
facs_markers <- colnames(sample_sheet_nn215)[facs_data_idx]

# Some final tidying.
sample_sheet_nn215 <- sample_sheet_nn215 %>%
  mutate(
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    sequencing_run = "NN215",
    # NOTE: Need to fix up a typo.
    illumina_index_index_number_separate_index_read =
      sub(
        "(PRI|RPI) ",
        "RPI-",
        illumina_index_index_number_separate_index_read)) %>%
  arrange(well_position)

# Construct NN227 sample sheet -------------------------------------------------

file_nn227 <- here(
  "data",
  "sample_sheets",
  "C094_DanPellici_MB NN227_Seqprimer10Jun21.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn227,
  sheet = "Samples",
  skip = 2,
  n_max = 1)

# NOTE: No FACS data
header_row <- paste0(colnames(header_row), header_row[1, ])
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn227 <- read_excel(
  path = file_nn227,
  sheet = "Samples",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)
# Tidy up names and empty rows/columns.
sample_sheet_nn227 <- clean_names(sample_sheet_nn227)
sample_sheet_nn227 <- remove_empty(
  sample_sheet_nn227,
  which = c("rows", "cols"))

# Filter out those wells without a cell index sequence, with no cell, or that
# were otherwise removed.
# NOTE: The 200-cell samples were not included in the pool sequenced in NN227
#       (the `sample_type` column says "200 cell (not included in NN226)" (sic)
#       for these samples).
sample_sheet_nn227 <- filter(
  sample_sheet_nn227,
  !is.na(plate_number),
  illumina_index_index_number_separate_index_read != "removed",
  sample_type != "200 cell (not included in NN226)")

# Some final tidying.
sample_sheet_nn227 <- sample_sheet_nn227 %>%
  mutate(
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    sequencing_run = "NN227",
    # NOTE: Need to fix up a typo.
    illumina_index_index_number_separate_index_read =
      sub("PRI ", "RPI-", illumina_index_index_number_separate_index_read)) %>%
  arrange(well_position)

# Construct final sample sheet -------------------------------------------------

sample_sheet <- bind_rows(sample_sheet_nn215, sample_sheet_nn227) %>%
  mutate(rowname = paste0(
    plate_number,
    "_",
    well_position)) %>%
  tibble::column_to_rownames("rowname") %>%
  DataFrame(., check.names = FALSE)

# NOTE: Check that there aren't any malformed well_positions (e.g., 'I19=A1').
stopifnot(!anyNA(sample_sheet$well_position))

# Key variables ----------------------------------------------------------------

rpis <- unique(sample_sheet$illumina_index_index_number_separate_index_read)
names(rpis) <- rpis
sequencing_runs <- tapply(
  sample_sheet$sequencing_run,
  sample_sheet$illumina_index_index_number_separate_index_read,
  unique)
outdir <- here("data", "SCEs")
dir.create(outdir, recursive = TRUE)
extdir <- here("extdata", "scPipe_trimmed_R2", sequencing_runs, rpis)
names(extdir) <- rpis
sapply(extdir, dir.create, recursive = TRUE)
# NOTE: Only using first 7 nt of barcode.
read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
# NOTE: Must be an element of biomaRt::listDatasets(), e.g.,
#       biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))[["dataset"]]
organism <- "hsapiens_gene_ensembl"
# NOTE: Must be an element of biomaRt::listAttributes(), e.g.,
#       biomaRt::listAttributes(biomaRt::useEnsembl("ensembl", organism))[["name"]]
gene_id_type <- "ensembl_gene_id"

# Input files ------------------------------------------------------------------

# FASTQ files
r1_fq <- c(
  grep(
    pattern = "Undetermined",
    x = list.files(
      path = here("extdata", "NN215", "merged"),
      full.names = TRUE,
      pattern = glob2rx("*R1*.fastq.gz")),
    invert = TRUE,
    value = TRUE),
  list.files(
    path = here("extdata", "NN227", "merged"),
    full.names = TRUE,
    pattern = glob2rx("LCE516-LCE531*R1*.fastq.gz")))
r2_fq <- gsub("R1", "R2", r1_fq)
stopifnot(all(file.exists(r2_fq)))
tx_fq <- file.path(extdir, paste0(rpis, ".R2.fastq.gz"))
names(tx_fq) <- rpis
barcode_fq <- gsub("\\.R2", "\\.R1", tx_fq)

mclapply(rpis, function(rpi) {
  message(rpi)
  if (rpi == "RPI-15") {
    cmd <- paste0(
      "cat ",
      grep("LCE516-LCE531-50", r1_fq, value = TRUE),
      " > ",
      barcode_fq[[rpi]],
      "\n",
      "cat ",
      grep("LCE516-LCE531-50", r2_fq, value = TRUE),
      " > ",
      tx_fq[[rpi]])
  } else if (rpi == "RPI-16") {
    cmd <- paste0(
      "cat ",
      grep("LCE516-LCE531-100", r1_fq, value = TRUE),
      " > ",
      barcode_fq[[rpi]],
      "\n",
      "cat ",
      grep("LCE516-LCE531-100", r2_fq, value = TRUE),
      " > ",
      tx_fq[[rpi]])
  } else if (rpi %in% paste0("RPI-", c(1, 2, 3, 5, 9, 10, 12))) {
    plate <- unique(
      sample_sheet[
        sample_sheet$illumina_index_index_number_separate_index_read == rpi,
        "plate_number"])
    cmd <- paste0(
      "cat ",
      grep(plate, r1_fq, value = TRUE),
      " > ",
      barcode_fq[[rpi]],
      "\n",
      "cat ",
      grep(plate, r2_fq, value = TRUE),
      " > ",
      tx_fq[[rpi]])
  } else {
    stop("Unknown RPI")
  }
  system(cmd)
})

# Genome index
genome_index <- here("extdata", "GRCh38.p13", "GRCh38.p13_with_ERCC")

# Genome annotation(s)
annofn <- c(
  here("extdata", "GRCh38.p13", "gencode.v35.primary_assembly.annotation.gff3"),
  system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))

# Cell barcodes
bc_anno <- file.path(extdir, paste0(rpis, ".barcode_annotation.csv"))
names(bc_anno) <- rpis

for (rpi in rpis) {
  message(rpi)
  tmp <- sample_sheet[
    sample_sheet$illumina_index_index_number_separate_index_read == rpi, ]
  barcode_df <- data.frame(
    cell_id = row.names(tmp),
    # NOTE: Only using first 7 nt of barcode.
    barcode = strtrim(
      tmp$rd1_index_cell_index_index_sequence_as_in_c_rt1_primer,
      7))
  stopifnot(!anyDuplicated(barcode_df$barcode))
  write.csv(
    x = barcode_df,
    file = bc_anno[[rpi]],
    quote = FALSE,
    row.names = FALSE)
}

# Output files -----------------------------------------------------------------

combined_fq <- file.path(extdir, gsub("R[12]", "combined", basename(tx_fq)))
names(combined_fq) <- names(tx_fq)
subread_bam <- gsub("fastq.gz", "subread.bam", combined_fq, fixed = TRUE)
exon_bam <- gsub("subread", "exon", subread_bam)

# FASTQ reformatting -----------------------------------------------------------

filter_settings <- list(rmlow = TRUE, rmN = FALSE, minq = 20, numbq = 2)
# NOTE: Have to loop over files because sc_trim_barcode() is not vectorised.
mclapply(seq_along(tx_fq), function(i) {
  message(combined_fq[i])
  sc_trim_barcode(
    outfq = combined_fq[i],
    r1 = tx_fq[i],
    r2 = barcode_fq[i],
    read_structure = read_structure,
    filter_settings = filter_settings)
})

# Aligning reads to a reference genome -----------------------------------------

# NOTE: Trimming R2 to simulate sequencing with the shorter read length
#       (nominally the 50-cycle kit but Stephen squeezes out much more than
#       that, apparently); current R2 is 71 cycles, new kit would be 88 in
#       total for I1 + R1 + R2 where I1=6 and R1=26.
align(
  index = genome_index,
  readfile1 = combined_fq,
  output_file = subread_bam,
  nTrim3 = 71 - (88 - 6 - 26),
  nthreads = 16)

# Assigning reads to annotated exons -------------------------------------------

bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
mclapply(seq_along(subread_bam), function(i) {
  message(i)
  sc_exon_mapping(
    inbam = subread_bam[i],
    outbam = exon_bam[i],
    annofn = annofn,
    bam_tags = bam_tags,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = stnd,
    fix_chr = fix_chr)
})

# De-multiplexing data ---------------------------------------------------------

max_mis <- 1
has_UMI <- TRUE
mito <- "chrM"
mclapply(seq_along(exon_bam), function(i) {
  message(i)
  sc_demultiplex(
    inbam = exon_bam[i],
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    max_mis = max_mis,
    bam_tags = bam_tags,
    mito = mito,
    has_UMI = has_UMI)
})

# Gene counting deduped data ---------------------------------------------------

UMI_cor <- 1
gene_fl <- FALSE
mclapply(seq_along(bc_anno), function(i) {
  message(i)
  sc_gene_counting(
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
})

# Create and save deduped SingleCellExperiment ---------------------------------

list_of_sce <- lapply(rpis, function(rpi) {
  create_sce_by_dir(
    datadir = extdir[[rpi]],
    organism = organism,
    gene_id_type = gene_id_type,
    pheno_data = sample_sheet[
      sample_sheet$illumina_index_index_number_separate_index_read == rpi, ],
    # NOTE: Create the report separately for more fine-grained control.
    report = FALSE)
})
sce <- Reduce(function(x, y) .combine(x, y, rowData_by = NULL), list_of_sce)
assay(sce, withDimnames = FALSE) <- as(
  assay(sce, withDimnames = FALSE),
  "dgCMatrix")
sce <- splitAltExps(
  sce,
  ifelse(grepl("^ERCC", rownames(sce)), "ERCC", "Endogenous"))
facs_data <- as.matrix(endoapply(colData(sce)[, facs_markers], as.numeric))
colData(sce)[, facs_markers] <- NULL
altExp(sce, "FACS") <- SummarizedExperiment(assays = list(raw = t(facs_data)))
saveRDS(
  sce,
  file.path(outdir, "C094_Pellicci.scPipe_trimmed_R2.SCE.rds"),
  compress = "xz")

# Create QC report of deduped data----------------------------------------------

library(readr)
library(plotly)
library(DT)
library(scran)
library(Rtsne)
# NOTE: Needs a fix for https://github.com/LuyiTian/scPipe/issues/100.
dir.create(here("output", "scPipe_trimmed_R2"), recursive = TRUE)
# NOTE: Tends to crap itself if using mclapply().
lapply(rpis, function(rpi) {
  try(create_report(
    sample_name = rpi,
    outdir = extdir[[rpi]],
    r1 = tx_fq[[rpi]],
    r2 = barcode_fq[[rpi]],
    outfq = combined_fq[[rpi]],
    read_structure = read_structure,
    filter_settings = filter_settings,
    align_bam = subread_bam[[rpi]],
    genome_index = genome_index,
    map_bam = exon_bam[[rpi]],
    exon_anno = annofn,
    stnd = stnd,
    fix_chr = fix_chr,
    barcode_anno = bc_anno[[rpi]],
    max_mis = max_mis,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl,
    organism = organism,
    gene_id_type = gene_id_type))

  # NOTE: Workaround bug in create_report() and stop output after 'Data summary'
  #       section.
  tmp <- readLines(file.path(extdir[[rpi]], "report.Rmd"))
  tmp <- c(tmp[1:161], "knitr::knit_exit()", tmp[162:length(tmp)])
  writeLines(tmp, file.path(extdir[[rpi]], "report.Rmd"))
  knitr::wrap_rmd(
    file = file.path(extdir[[rpi]], "report.Rmd"),
    width = 120,
    backup = NULL)
  rmarkdown::render(
    input = file.path(extdir[[rpi]], "report.Rmd"),
    output_file = file.path(extdir[[rpi]], "report.html"),
    knit_root_dir = ".")

  # NOTE: Copy the QC report to the repository.
  file.copy(
    from = file.path(extdir[[rpi]], "report.nb.html"),
    to = here(
      "output",
      "scPipe_trimmed_R2",
      paste0(rpi, ".scPipe_QC_report.nb.html")),
    overwrite = TRUE)
})

# Gene counting non-deduped data -----------------------------------------------

sapply(paste0(extdir, "_no_dedup"), dir.create)
lapply(names(bc_anno), function(n) {
  message(n)
  geneCountingNoUMIDedup(
    outdir = extdir[[n]],
    bc_anno = bc_anno[[n]])
})

# Create and save non-deduped SingleCellExperiment -----------------------------

list_of_no_dedup_sce <- lapply(rpis, function(rpi) {
  x <- data.table::fread(
    file.path(paste0(extdir[[rpi]], "_no_dedup"), "gene_count.csv"))
  counts <- as.matrix(x[, -1])
  sce <- SingleCellExperiment(
    list(counts = counts),
    colData = sample_sheet[colnames(counts), ])
  rownames(sce) <- x$gene_id
  sce
})

no_dedup_sce <- Reduce(
  function(x, y) .combine(x, y, rowData_by = NULL),
  list_of_no_dedup_sce)
no_dedup_sce$UMI_deduped <- FALSE
colnames(no_dedup_sce) <- paste0(colnames(no_dedup_sce), ".not_UMI_deduped")
assay(no_dedup_sce, withDimnames = FALSE) <- unname(
  as(assay(no_dedup_sce, withDimnames = FALSE), "dgCMatrix"))
no_dedup_sce <- splitAltExps(
  no_dedup_sce,
  ifelse(grepl("^ERCC", rownames(no_dedup_sce)), "ERCC", "Endogenous"))
# NOTE: Order genes and samples as in `sce`.
no_dedup_sce <- no_dedup_sce[rownames(sce),
                             paste0(colnames(sce), ".not_UMI_deduped")]
saveRDS(
  no_dedup_sce,
  file.path(outdir, "C094_Pellicci.not_UMI_deduped.scPipe_trimmed_R2.SCE.rds"),
  compress = "xz")
