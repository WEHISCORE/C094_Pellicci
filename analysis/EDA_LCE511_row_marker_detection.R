library(here)
library(BiocStyle)
library(dplyr)
library(janitor)
library(scater)
library(patchwork)
library(scales)
source(here("code", "helper_functions.R"))

library(SingleCellExperiment)
sce <- readRDS("data/SCEs/C094_Pellicci.scPipe.SCE.rds")

library(dplyr)
# NOTE: Some `sample_name` include `sample_gate`; remove this redundancy.
sce$sample_name <- case_when(
  grepl("Cell line", sce$sample_name) ~ "Cell line",
  TRUE ~ sub(" P[0-9]$", "", sce$sample_name))
sce$sample_gate[is.na(sce$sample_gate)] <- "NA"
# NOTE: Some `sample_gate` have a `(stage X)` suffix repeated 1-2 times; remove
#       this unnecessary information.
sce$sample_gate <- strtrim(sce$sample_gate, 2)
# NOTE: `sample_name` = <tissue>_<donor>
sce$tissue <- factor(
  dplyr::case_when(
    grepl("Blood", sce$sample_name) ~ "Blood",
    grepl("Thymus", sce$sample_name) ~ "Thymus",
    sce$sample_name == "Cell line" ~ "SKW3"),
  levels = c("Blood", "Thymus", "SKW3"))
sce$donor <- factor(
  dplyr::case_when(
    grepl("1", sce$sample_name) ~ "1",
    grepl("2", sce$sample_name) ~ "2",
    grepl("3", sce$sample_name) ~ "3",
    grepl("4", sce$sample_name) ~ "4",
    grepl("5", sce$sample_name) ~ "5",
    sce$sample_name == "Cell line" ~ "NA"),
  levels = c("1", "2", "3", "4", "5", "NA"))
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

# Some useful colours
sce$colours <- S4Vectors::make_zero_col_DFrame(ncol(sce))
sample_type_colours <- setNames(
  palette.colors(nlevels(sce$sample_type), "R4"),
  levels(sce$sample_type))
sce$colours$sample_type_colours <- sample_type_colours[sce$sample_type]
sample_name_colours <- setNames(
  palette.colors(nlevels(sce$sample_name), "Alphabet"),
  levels(sce$sample_name))
sce$colours$sample_name_colours <- sample_name_colours[sce$sample_name]
sample_gate_colours <- setNames(
  palette.colors(nlevels(sce$sample_gate), "Tableau 10"),
  levels(sce$sample_gate))
sce$colours$sample_gate_colours <- sample_gate_colours[sce$sample_gate]
plate_number_colours <- setNames(
  palette.colors(nlevels(sce$plate_number), "Okabe-Ito"),
  levels(sce$plate_number))
sce$colours$plate_number_colours <- plate_number_colours[sce$plate_number]
tissue_colours <- setNames(
  palette.colors(nlevels(sce$tissue), "ggplot2"),
  levels(sce$tissue))
sce$colours$tissue_colours <- tissue_colours[sce$tissue]
donor_colours <- setNames(
  palette.colors(nlevels(sce$donor), "Set 1"),
  levels(sce$donor))
sce$colours$donor_colours <- donor_colours[sce$donor]

pseudoLog <- scales::pseudo_log_trans(sigma = 150 / 2)$transform
assay(altExp(sce, "FACS"), "pseudolog") <- pseudoLog(
  assay(altExp(sce, "FACS"), "raw"))

library(here)
library(readxl)
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
    post_hoc_sample_gate = sample_gate) %>%
  arrange(well_position)

# NOTE: A poor-man's inner_join() to avoid the lossy coercion of colData(sce)
#       to a data.frame.
j <- match(
  colnames(sce),
  paste0(sample_sheet_nn215$plate_number, "_", sample_sheet_nn215$well_position))
colData(sce) <- cbind(
  colData(sce),
  sample_sheet_nn215[j, "post_hoc_sample_gate", drop = FALSE])
# Tidy up mini-bulk samples (which were supplied with the gate info).
sce$post_hoc_sample_gate[sce$sample_type != "Single cell"] <-
  as.character(sce$sample_gate[sce$sample_type != "Single cell"])
sce$post_hoc_sample_gate <- factor(sce$post_hoc_sample_gate)

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

sce <- logNormCounts(sce)

sizeFactors(altExp(sce, "ERCC")) <- librarySizeFactors(altExp(sce, "ERCC"))





### these are the `Thymus 2` and `Blood 2`  unique markers that defined based only on plates other than `LCE511` (and not the 2 new plates, i.e. LCE515 and 516)

# sample_name-specific upregulated genes
# select all plate except LCE511, and NOT LCE515 and LCE516 (mini-bulk only)
sce1 <- sce[,!(sce$plate_number=="LCE511" | sce$plate_number=="LCE515" | sce$plate_number=="LCE516")]
sce1$plate_number <- droplevels(sce1$plate_number)

# limit to only Blood 2 and Thymus 2
sce1 <- sce1[,(sce1$sample_name == "Blood 2" | sce1$sample_name == "Thymus 2")]
sce1$sample_name <- droplevels(sce1$sample_name)

# find unique markers for thymus2 and blood2
sample_name_markers <- findMarkers(
  sce1,
  groups = sce1$sample_name,
  direction = "up",
  pval.type = "all",
  row.data = rowData(sce1))

# sumup table
print(sapply(sample_name_markers, function(x) sum(x$FDR < 0.05)))

# plot
collected <- list()
for (lab in names(sample_name_markers)) {
  lab_markers <- sample_name_markers[[lab]]
  # NOTE: Select top-50 markers for plotting.
  m <- head(rownames(lab_markers), 50)
  collected[[lab]] <- plotHeatmap(
    object = sce1,
    features = m,
    color = hcl.colors(101, "Blue-Red 3"),
    center = TRUE,
    zlim = c(-3, 3),
    order_columns_by = c("sample_name", "plate_number", "post_hoc_sample_gate", "tissue", "donor"),
    cluster_rows = TRUE,
    fontsize = 5,
    column_annotation_colors = list(
      sample_name = sample_name_colours,
      post_hoc_sample_gate = sample_gate_colours,
      plate_number = plate_number_colours,
      tissue = tissue_colours,
      donor = donor_colours),
    main = lab,
    # annotation_legend = FALSE,
    silent = TRUE)[[4]]
}
wrap_plots(collected, ncol = 1)

# plot (if grouping the columns)
collected <- list()
for (lab in names(sample_name_markers)) {
  lab_markers <- sample_name_markers[[lab]]
  # NOTE: Select top-50 markers for plotting.
  m <- head(rownames(lab_markers), 50)
  collected[[lab]] <- plotHeatmap(
    object = sce1,
    features = m,
    color = hcl.colors(101, "Blue-Red 3"),
    center = TRUE,
    zlim = c(-3, 3),
    colour_columns_by = c("sample_name", "plate_number", "post_hoc_sample_gate", "tissue", "donor"),
    cluster_rows = TRUE,
    fontsize = 5,
    column_annotation_colors = list(
      sample_name = sample_name_colours,
      post_hoc_sample_gate = sample_gate_colours,
      plate_number = plate_number_colours,
      tissue = tissue_colours,
      donor = donor_colours),
    main = lab,
    # annotation_legend = FALSE,
    silent = TRUE)[[4]]
}
wrap_plots(collected, ncol = 1)




### then if we make use of this `Thymus 2` and `Blood 2`-specific markers and look into the last 5 rows of the plate LCE511

# subset to plate- and row-of-interest
sce2 <- sce[, sce$plate_number == "LCE511" & substr(sce$well_position, 1, 1) %in% c("L", "M", "N", "O", "P")]
sce2$row <- factor(sub("^([[:alpha:]]*).*", "\\1", sce2$well_position))

# set row colours
library(Polychrome)
row_colours <- setNames(
  Polychrome::kelly.colors(nlevels(sce2$row)+1)[-1],
  levels(sce2$row))
sce2$row_colours <- row_colours[sce2$row]

# plot
collected <- list()
for (lab in names(sample_name_markers)) {
  lab_markers <- sample_name_markers[[lab]]
  # NOTE: Select top-50 markers for plotting.
  m <- head(rownames(lab_markers), 50)
  collected[[lab]] <- plotHeatmap(
    object = sce2,
    features = m,
    color = hcl.colors(101, "Blue-Red 3"),
    center = TRUE,
    zlim = c(-3, 3),
    order_columns_by = c("row", "sample_name", "plate_number", "post_hoc_sample_gate", "tissue", "donor"),
    cluster_rows = TRUE,
    fontsize = 5,
    column_annotation_colors = list(
      row = row_colours,
      sample_name = sample_name_colours,
      post_hoc_sample_gate = sample_gate_colours,
      plate_number = plate_number_colours,
      tissue = tissue_colours,
      donor = donor_colours),
    main = lab,
    # annotation_legend = FALSE,
    silent = TRUE)[[4]]
}
wrap_plots(collected, ncol = 1)

# plot (if grouping the columns)
collected <- list()
for (lab in names(sample_name_markers)) {
  lab_markers <- sample_name_markers[[lab]]
  # NOTE: Select top-50 markers for plotting.
  m <- head(rownames(lab_markers), 50)
  collected[[lab]] <- plotHeatmap(
    object = sce2,
    features = m,
    color = hcl.colors(101, "Blue-Red 3"),
    center = TRUE,
    zlim = c(-3, 3),
    colour_columns_by = c("row", "sample_name", "plate_number", "post_hoc_sample_gate", "tissue", "donor"),
    cluster_rows = TRUE,
    fontsize = 5,
    column_annotation_colors = list(
      row = row_colours,
      sample_name = sample_name_colours,
      post_hoc_sample_gate = sample_gate_colours,
      plate_number = plate_number_colours,
      tissue = tissue_colours,
      donor = donor_colours),
    main = lab,
    # annotation_legend = FALSE,
    silent = TRUE)[[4]]
}
wrap_plots(collected, ncol = 1)

# comments:
# so our hypothesis is, row P (is Blood 2) sample accidentally being sorted onto row O (Thymus 2),
# it means we are supposed to see, on row O, both Blood 2 and Thymus 2  markers if this mixing really occurs;
# if not, for row O, we should only see Thymus 2 markers only
#
# if based on the heatmap without grouping of cells by colour,
# it seems to be true that row O have expression of both Thymus 2 and Blood 2 markers, but same is true for row L to row N
# Unless the contamination with Blood 2 occur far "upward" that cover all row L  to row N  samples,
# I am afraid based on this, it still hard to tell if row O was contaminated by row P  based on these markers
#
# but if based heatmap grouped cells by colour
# the thymus 2  markers seems not seeing consistency of marker expression between cells for each row (which may not be of useful);
# but for blood 2  one , it seems like cells on row O (orange) are showing up-regulation of blood 2 markers,
# indicating that blood 2 sample on row P goes up to row O  (???)







