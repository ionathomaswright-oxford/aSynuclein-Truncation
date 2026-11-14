# File description --------------------------------------------------------

### Goal
# make_feature_metadata.R produces gene and transcript metadata from GENECODE annotation.

### Execution steps
# ===>>> rwm/make_feature_metadata.R
# - rwm/salmon.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Iona Thomas-Wright (iona.thomaswright@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Thomas-Wright_2026/data/rwm/make_feature_metadata"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Define annotation parameters --------------------------------------------

species <- "human"

release <- 46

# Update paths ------------------------------------------------------------

sr_dir <- ifelse(species == "mouse",
                 paste0("vM", release),
                 paste0("v", release))

out_dir_sub_path <- file.path(out_dir_path, sr_dir)

# Download GENCODE annotation ---------------------------------------------

download_gencode_annotation(species = species,
                            release = release,
                            out_dir_path = out_dir_sub_path)

# Import GENCODE annotation -----------------------------------------------

gencode_annotation <- import_gencode_annotation(file = list.files(path = out_dir_sub_path,
                                                                  pattern = "^.*\\.gtf\\.gz$",
                                                                  full.names = TRUE))

# Add transgene annotation ------------------------------------------------

transgenes <- c("mTagBFP2", "puromycinR")

transgene_annotation <- tibble(gene_id = transgenes,
                               transcript_id = transgenes,
                               gene_name = transgenes,
                               gene_type = "protein_coding")

gencode_annotation <- bind_rows(gencode_annotation, transgene_annotation)

write_csv(gencode_annotation, file.path(out_dir_sub_path, "gencode_annotation.csv"))

# Make feature metadata ---------------------------------------------------

gene_metadata <- make_gene_metadata(input = gencode_annotation)

transcript_metadata <- make_transcript_metadata(input = gencode_annotation)

# Write feature metadata to file ------------------------------------------

write_csv(gene_metadata, file.path(out_dir_sub_path, "gene_metadata.csv"))

write_csv(transcript_metadata, file.path(out_dir_sub_path, "transcript_metadata.csv"))

# Make SNCA metadata ------------------------------------------------------

# Construct SNCA metadata from GENCODE annotation
snca_metadata <- gencode_annotation %>%
  filter(gene_name == "SNCA") |> 
  unique() |>
  select(chr_name = seqnames,
         gene_symbol = gene_name,
         gene_type,
         ensembl_gene_id_version = gene_id,
         transcript_name,
         transcript_type,
         ensembl_transcript_id_version = transcript_id,
         ensembl_protein_id_version = protein_id,
         type, exon_id, exon_number, start, end, width, strand) |>
  nest(data = c(type, exon_id, exon_number, start, end, width)) |>
  mutate(transcript_data = map(.x = data,
                               .f = ~.x[.x$type == "transcript", c("start", "end", "width")] |>
                                 rename_with(.fn = function(x) paste0(x, "_transcript"),
                                             .cols = everything())),
         exon_data = map2(.x = data, .y = strand,
                          .f = ~ {
                            dat <- .x[.x$type != "transcript", ]
                            cds_start <- suppressWarnings(min(dat$start[dat$type == "CDS"], na.rm = TRUE))
                            cds_end <- suppressWarnings(max(dat$end[dat$type == "CDS"], na.rm = TRUE))
                            
                            dat <- dat |>
                              mutate(type = case_when(
                                type != "UTR" ~ as.character(type),
                                .y == "+" & end < cds_start ~ "UTR5",
                                .y == "+" & start > cds_end ~ "UTR3",
                                .y == "-" & start > cds_end ~ "UTR5",
                                .y == "-" & end < cds_start ~ "UTR3",
                                TRUE ~ as.character(type)
                              ))
                            
                            dat |>
                              nest(data = -c(exon_id, exon_number)) |>
                              mutate(data = map(.x = data,
                                                .f = ~ pivot_wider(data = .x, names_from = "type", names_sep = "_",
                                                                   values_from = c("start", "end", "width"),
                                                                   names_vary = "slowest"))) |>
                              unnest(cols = data) |>
                              rename_with(.fn = tolower, .cols = everything())
                          }),
         n_exon = map_int(.x = exon_data,
                          .f = ~ dim(.x)[[1]]),
         width_cds = map_int(.x = exon_data,
                             .f = ~ ifelse("width_cds" %in% colnames(.x),
                                           sum(.x[["width_cds"]], na.rm = T),
                                           NA_integer_)),
         width_peptide = ceiling(width_cds / 3),
         cds_complete = ifelse(width_cds %% 3 == 0, T, F),
         transcript_id = ifelse(is.na(width_peptide),
                                paste0(transcript_name, " CDS unknown"),
                                paste0(transcript_name, " ", width_peptide, "aa")),
         transcript_id_split = str_replace(transcript_id, " ", "\n")) |>
  unnest(cols = transcript_data) |> 
  select(transcript_id, everything(), -c(data, exon_data)) |> 
  arrange(ensembl_transcript_id_version)

# Retrieve nucleic acid and protein sequences using biomaRt
biomaRt::listMarts()

biomart_db <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")

biomaRt::listDatasets(mart = biomart_db) %>% dplyr::filter(str_detect(dataset, "sapiens"))

biomart_ds <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filters <- biomaRt::searchFilters(mart = biomart_ds, pattern = "ensembl")

attributes <- biomaRt::searchAttributes(mart = biomart_ds, pattern = "sequence")

attributes <- attributes %>%
  dplyr::filter(page == "sequences",
                !str_detect(name, "flank")) %>%
  select(name, description) %>%
  unique() %>%
  group_by(name) %>%
  mutate(description = paste0(description, collapse = " / ")) %>%
  ungroup() %>%
  unique()

snca_sequences <- lapply(X = c("5utr", "3utr", "cdna", "coding", "peptide"),
                         FUN = function(x) {
                           biomaRt::getBM(attributes = c("ensembl_transcript_id_version", x),
                                          filters = "ensembl_transcript_id_version",
                                          values = snca_metadata$ensembl_transcript_id_version,
                                          mart = biomart_ds)
                         })

snca_sequences <- Reduce(f = function(...) { merge(..., all = T) },
                         x = snca_sequences)

names(snca_sequences) <- c("ensembl_transcript_id_version",
                           "5utr_seq", "3utr_seq", "cdna_seq", "coding_seq", "peptide_seq")

# Update SNCA metadata with sequences
snca_metadata <- snca_metadata %>%
  left_join(snca_sequences,
            by = join_by(ensembl_transcript_id_version))

write_csv(snca_metadata, file.path(out_dir_sub_path, "snca_metadata.csv"))

