# Pre-processing that is common to the single-cell and mini-bulk samples
# Peter Hickey
# 2021-06-30

# This script is complicated and ad hoc because the sample sheets are an
# inconsistent mess that require a lot of normalization and post-processing.

library(here)
library(SingleCellExperiment)
library(dplyr)
library(readxl)
library(janitor)

# Setting up the data ----------------------------------------------------------

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.scPipe.SCE.rds"))

# Incorporating cell-based annotation ------------------------------------------

# NOTE: Some `sample_name` include `sample_gate`; remove this redundancy.
sce$sample_name <- case_when(
  grepl("Cell line", sce$sample_name) ~ "Cell line",
  TRUE ~ sub(" P[0-9]$", "", sce$sample_name))

# NOTE: Some `sample_gate` have a `(stage X)` suffix repeated 1-2 times; remove
#       this unnecessary information.
sce$sample_gate[is.na(sce$sample_gate)] <- "NA"
sce$sample_gate <- strtrim(sce$sample_gate, 2)
# NOTE: P2 and P3 are gates used for the SKW3 and are meaningless;
sce$sample_gate[sce$sample_gate %in% c("P2", "P3")] <- "NA"

# NOTE: `sample_name` = <tissue>_<donor>
sce$tissue <- factor(
  dplyr::case_when(
    grepl("Blood", sce$sample_name) ~ "Blood",
    grepl("Thymus", sce$sample_name) ~ "Thymus",
    sce$sample_name == "Cell line" ~ "SKW3"),
  levels = c("Blood", "Thymus", "SKW3"))
# NOTE: Although there are donors 1-3 in both NN225 and NN227, these are not
#       the same donor. There are in fact 8 distinct donors across the 2
#       batches.
sce$donor <- factor(
  dplyr::case_when(
    grepl("1", sce$sample_name) & sce$sequencing_run == "NN215" ~ "1",
    grepl("2", sce$sample_name) & sce$sequencing_run == "NN215" ~ "2",
    grepl("3", sce$sample_name) & sce$sequencing_run == "NN215" ~ "3",
    grepl("1", sce$sample_name) & sce$sequencing_run == "NN227" ~ "4",
    grepl("2", sce$sample_name) & sce$sequencing_run == "NN227" ~ "5",
    grepl("3", sce$sample_name) & sce$sequencing_run == "NN227" ~ "6",
    grepl("4", sce$sample_name) & sce$sequencing_run == "NN227" ~ "7",
    grepl("5", sce$sample_name) & sce$sequencing_run == "NN227" ~ "8",
    sce$sample_name == "Cell line" ~ "NA"),
  levels = c("1", "2", "3", "4", "5", "6", "7", "8", "NA"))

# Create factor columns, which are easier to work with.
# NOTE: Manually order some factors' levels, others we don't care about the
#       level order.
sce$sample_type <- factor(
  case_when(
    sce$sample_type %in% c("50 cell", "50 cells") ~ "50 cells",
    sce$sample_type %in% c("100 cell", "100 cells") ~ "100 cells",
    sce$sample_type == "99 cell" ~ "99 cells",
    TRUE ~ sce$sample_type),
  c("Single cell", "50 cells", "99 cells", "100 cells"))
sce$sample_name <- factor(
  sce$sample_name,
  levels = c(paste("Blood", 1:5), paste("Thymus", 1:5), "Cell line"))
colData(sce) <- DataFrame(
  endoapply(colData(sce), function(x) {
    if (is.character(x)) {
      as.factor(x)
    } else {
      x
    }
  }),
  row.names = colnames(sce))

# Incorporating FACS data ------------------------------------------------------

pseudoLog <- scales::pseudo_log_trans(sigma = 150 / 2)$transform
assay(altExp(sce, "FACS"), "pseudolog") <- pseudoLog(
  assay(altExp(sce, "FACS"), "raw"))

# Post-hoc gating of single cells into stages ----------------------------------

file_nn215 <- here(
  "data",
  "sample_sheets",
  "C094_DanPellici_SC_MB_NN215_AdjustedGates.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn215,
  sheet = "Samples and indexing",
  skip = 2,
  n_max = 1)

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
sample_sheet_nn215 <- dplyr::filter(
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
    plate_number = factor(plate_number, levels = levels(sce$plate_number)),
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
    post_hoc_sample_gate = ifelse(is.na(sample_gate), "NA", sample_gate)) %>%
  arrange(well_position)

# NOTE: P5 is the general gate such that (P6, P7, P8) are non-intersecting
#       subsets of P5.
gate_to_stage_df <- data.frame(
  gate = c("NA",
           "P5",
           "P6", "P7", "P8"),
  stage = c("NA",
            "Unknown",
            "S1 (CD4+/CD161-)", "S2 (CD4-/CD161-)", "S3 (CD4-/CD161+)"))

# Add `stage` to colData(sce).
colData(sce) <- left_join(
  as.data.frame(colData(sce)),
  select(sample_sheet_nn215, plate_number, well_position, post_hoc_sample_gate),
  by = c("plate_number", "well_position")) %>%
  mutate(
    sample_gate = factor(sample_gate, levels = c("P5", "P6", "P7", "P8", "NA")),
    post_hoc_sample_gate = ifelse(
      is.na(post_hoc_sample_gate),
      as.character(sample_gate),
      post_hoc_sample_gate)) %>%
  inner_join(
    gate_to_stage_df,
    by = c("post_hoc_sample_gate" = "gate")) %>%
  select(
    plate_number, well_position, sample_type, sample_gate, sequencing_run, tissue, donor, stage,
    unaligned, aligned_unmapped, mapped_to_exon, mapped_to_intron, ambiguous_mapping, mapped_to_ERCC, mapped_to_MT) %>%
  DataFrame(row.names = colnames(sce))
sce$stage <- factor(
  sce$stage,
  levels = c("S1 (CD4+/CD161-)", "S2 (CD4-/CD161-)", "S3 (CD4-/CD161+)", "Unknown", "NA"))

# Some useful colours ----------------------------------------------------------

sce$colours <- S4Vectors::make_zero_col_DFrame(ncol(sce))
plate_number_colours <- setNames(
  palette.colors(nlevels(sce$plate_number), "Okabe-Ito"),
  levels(sce$plate_number))
sce$colours$plate_number_colours <- plate_number_colours[sce$plate_number]
sample_type_colours <- setNames(
  palette.colors(nlevels(sce$sample_type), "R4"),
  levels(sce$sample_type))
sce$colours$sample_type_colours <- sample_type_colours[sce$sample_type]
sample_gate_colours <- setNames(
  palette.colors(nlevels(sce$sample_gate), "Dark 2"),
  levels(sce$sample_gate))
sce$colours$sample_gate_colours <- sample_gate_colours[sce$sample_gate]
sequencing_run_colours <- c("NN215" = "orange", "NN227" = "dodgerBlue")
sce$colours$sequencing_run_colours <- sequencing_run_colours[sce$sequencing_run]
tissue_colours <- setNames(
  palette.colors(nlevels(sce$tissue), "ggplot2"),
  levels(sce$tissue))
sce$colours$tissue_colours <- tissue_colours[sce$tissue]
donor_colours <- setNames(
  palette.colors(nlevels(sce$donor), "Tableau 10"),
  levels(sce$donor))
sce$colours$donor_colours <- donor_colours[sce$donor]
stage_colours <- c(
  "S1 (CD4+/CD161-)" = unname(sample_gate_colours["P6"]),
  "S2 (CD4-/CD161-)" = unname(sample_gate_colours["P7"]),
  "S3 (CD4-/CD161+)" = unname(sample_gate_colours["P8"]),
  "Unknown" = "black",
  "NA" = unname(sample_gate_colours["NA"]))
sce$colours$stage_colours <- stage_colours[sce$stage]

# Incorporating gene-based annotation ------------------------------------------

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
