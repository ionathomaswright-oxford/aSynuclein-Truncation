# File description --------------------------------------------------------

### Goal
# make_sample_metadata.R constructs the sample metadata for the FOUNDIN-PD bulk RNA-seq dataset.

### Execution steps
# ===>>> foundinpd/make_sample_metadata.R
# - foundinpd/make_feature_metadata.R
# - foundinpd/salmon.R
# - foundinpd/dgea.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - Participant metadata file (csv) available on GitHub repository 'XXX', at metadata/foundinpd/Participant_Status_17Sep2024.csv
# => Can be obtained from PPMI at https://ida.loni.usc.edu/web/ppmi-rnaseq-app (Download > Study Data > Subject Characteristics > Patient Status)
# - iPSC metadata file (csv) available on GitHub repository 'XXX', at metadata/foundinpd/iPSC_Catalog_Metadata_17Sep2024.csv
# => Can be obtained from PPMI at https://ida.loni.usc.edu/web/ppmi-rnaseq-app (Download > Study Data > Biospecimen > Bioseample Inventory > iPSC Catalog Metadata)
# - Salmon's quant.sf files available from PPMI upon request
# => Sequencing reads were mapped to custom human reference transcriptome (GRCh38.p12 - GENCODE Release 29 combined with LNCipedia version 5.2) using Salmon and shared by the PPMI via an Aspera transfer platform
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Iona Thomas-Wright (iona.thomaswright@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <-  "publications/Thomas-Wright_2026/data/foundinpd/make_sample_metadata" 

# Path to participant metadata file (csv)
# Obtained from PPMI at https://ida.loni.usc.edu/web/ppmi-rnaseq-app
# Download > Study Data > Subject Characteristics > Patient Status
participant_metadata_path <- "publications/Thomas-Wright_2026/metadata/foundinpd/Participant_Status_17Sep2024.csv"

# Path to iPSC metadata file (csv)
# Obtained from PPMI at https://ida.loni.usc.edu/web/ppmi-rnaseq-app
# Download > Study Data > Biospecimen > Bioseample Inventory > iPSC Catalog Metadata
ipsc_metadata_path <- "publications/Thomas-Wright_2026/metadata/foundinpd/iPSC_Catalog_Metadata_17Sep2024.csv"

# Path to directory containing Salmon quant files
# => Obtained from PPMI via Aspera transfer
salmon_quant_dir_path <- "/scratch/ben/rnaseq/seq_data/public/foundinpd/processed/RNAB/salmon_quant"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Construct participant metadata ------------------------------------------

participant_metadata <- read_csv(participant_metadata_path) %>% 
  select(patno = PATNO,
         disease_status = COHORT_DEFINITION,
         starts_with("ENRL")) %>% 
  mutate(patno = as.character(patno),
         disease_status = case_when(disease_status == "Healthy Control" ~ "CTRL",
                                    disease_status == "Parkinson's Disease" ~ "PD",
                                    disease_status %in% c("Prodromal", "SWEDD") ~ disease_status,
                                    TRUE ~ NA_character_) %>% 
           fct_relevel(c("CTRL", "Prodromal", "PD", "SWEDD"))) %>% 
  pivot_longer(cols = starts_with("ENRL"), names_to = "risk_name", values_to = "risk_value") %>% 
  mutate(risk_factor = map2_chr(.x = risk_value,
                                .y = risk_name,
                                .f = ~ ifelse(.x == 1, 
                                              str_remove(.y, "ENRL"), 
                                              NA_character_))) %>% 
  select(-c(risk_name, risk_value)) %>% 
  nest(data = risk_factor) %>% 
  mutate(risk_factor = map_chr(.x = data,
                               .f = ~ paste0(na.omit(.x$risk_factor), collapse = " / ")) %>% 
           {ifelse(. == "", NA_character_, .)}) %>% 
  select(-data)

# Construct iPSC metadata -------------------------------------------------

ipsc_metdata <- read_csv(ipsc_metadata_path) %>% 
  select(patno = PATNO,
         disease_status = COHORT,
         gender = GENDER,
         age_at_diagnosis = AGE_AT_DIAGNOSIS,
         age_at_sample_collection = AGE_AT_SAMPLE_COLLECTION,
         karyotype = KARYOTYPE) %>% 
  mutate(patno = as.character(patno),
         disease_status = case_when(disease_status == "Healthy Control" ~ "CTRL",
                                    disease_status == "Parkinson's Disease" ~ "PD",
                                    disease_status %in% c("Prodromal", "SWEDD") ~ disease_status,
                                    TRUE ~ NA_character_) %>% 
           fct_relevel(c("CTRL", "Prodromal", "PD", "SWEDD"))) %>% 
  unique() %>% 
  filter(karyotype != "Abnormal")

# Construct sample metadata -----------------------------------------------

sample_metadata <- tibble(sample_spe_dir = list.dirs(path = salmon_quant_dir_path, 
                                                     full.names = F, recursive = F)) %>% 
  mutate(files = setNames(
    object = map_chr(.x = sample_spe_dir,
                     .f = ~ list.files(path = file.path(salmon_quant_dir_path, .x),
                                       pattern = "quant.sf",
                                       full.names = TRUE)), 
    nm = sample_spe_dir),
    split_name = str_split(string = sample_spe_dir, pattern = "_"),
    patno = map_chr(.x = split_name, .f = ~ .x[[2]]) %>% 
      str_remove("^PPMI") %>% str_remove("B.*$"),
    cdi_number = map_chr(.x = split_name, .f = ~ .x[[3]]),
    patno_cdi_number = paste0(patno, "_", cdi_number),
    batch = map_chr(.x = split_name, .f = ~ .x[[2]]) %>% 
      str_extract("B.*$") %>% {ifelse(is.na(.), "B1", .)},
    div = map_chr(.x = split_name, .f = ~ .x[[4]]) %>% 
      str_replace("^[a-z]+", "DIV"),
    version = map_chr(.x = split_name, .f = ~ .x[[5]])) %>% 
  select(-split_name) %>% 
  left_join(ipsc_metdata, by = join_by(patno)) %>% 
  left_join(participant_metadata, by = join_by(patno, disease_status)) 

write_rds(sample_metadata, file.path(out_dir_path, "sample_metadata.rds"))

# Subset sample metadata 
# => For each combination of patno and DIV, keep only the sample-specific directories:
# - matching the "da" differentiation protocol
# - the sample with smallest version index among samples with smallest batch index
# => Keep only the samples with disease status CTRL, Prodromal or PD

filtered_files <- sample_metadata %>% 
  filter(str_detect(files, "_da\\d{1,2}_v\\d{1,}/quant.sf$")) %>%
  select(patno, batch, div, version, files) %>% 
  unique() %>% 
  arrange(patno, batch, div, version) %>% 
  nest(data = c(batch, version, files)) %>% 
  mutate(data = map(data, ~ .x[1,])) %>% 
  unnest(data) %>% 
  pull(files)

filtered_sample_metadata <- sample_metadata %>% 
  filter(files %in% filtered_files,
         disease_status %in% c("CTRL", "Prodromal", "PD"))

write_rds(filtered_sample_metadata, file.path(out_dir_path, "filtered_sample_metadata.rds"))


