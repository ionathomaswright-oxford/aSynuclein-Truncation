# File description --------------------------------------------------------

### Goal
# dgea.R performs differential gene expression analysis (DGEA) between biological conditions, focusing on SNCA.

### Execution steps
# x foundinpd/make_sample_metadata.R
# x foundinpd/make_feature_metadata.R
# x foundinpd/salmon.R
# ===>>> foundinpd/dgea.R

# => Any script listed before the current file must be executed beforehand 

### Requirements
# - filtered_sample_metadata.rds produced by make_sample_metadata.R
# - gene_metadata.csv produced by make_feature_metadata.R
# - transcript_metadata.csv produced by make_feature_metadata.R
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Iona Thomas-Wright (iona.thomaswright@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Thomas-Wright_2026/data/foundinpd/dgea"

# Path to sample metadata file (csv)
sample_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_sample_metadata/filtered_sample_metadata.rds"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_feature_metadata/v29/gene_metadata.csv"

# Path to transcript metadata file (csv)
transcript_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_feature_metadata/v29/transcript_metadata.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(rnaseqtools)
library(writexl)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Set n_workers for parallel processing
n_workers <- 12

# Define sample filtering parameters --------------------------------------

# Define DIV: select one from c("DIV0", "DIV25", "DIV65"), example for DIV65
div <- c("DIV65")

# Define DGE analysis parameters ------------------------------------------

# Define model_formula 
model_formula <- "~ disease_status"

# Define reduced_model_formula
reduced_model_formula <- ~ 1

# Define coef_var
coef_var <- "disease_status"

# Define coef_ref_lvl
coef_ref_lvl <- "CTRL"

# Define if protein-coding genes only
protein_coding_only <- TRUE

# Define adjusted p-value cutoff (FDR)
padj_dge <- 0.05

# Define min_cnt_excl threshold (exclusive; to be detected in n samples >= size of smaller group)
min_cnt_excl <- 0

# Load metadata -----------------------------------------------------------

# Load transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path)

# Load gene metadata
gene_metadata <- read_csv(gene_metadata_path)

# Import filtered sample metadata
sample_metadata <- read_rds(sample_metadata_path)

# Filter cell metadata ----------------------------------------------------

# Filter on DIV
sample_metadata <- sample_metadata[sample_metadata$div %in% div,]

# Set coef_var levels
sample_metadata[[coef_var]] = fct_relevel(sample_metadata[[coef_var]], coef_ref_lvl)

# Drop unused levels in categorical variables
sample_metadata <- sample_metadata %>%
  mutate(across(.cols = where(is.factor),
                .fns = droplevels))

# # Convert character variables present in model_formula to factor
model_vars <- colnames(sample_metadata)[map_lgl(.x = colnames(sample_metadata),
                                                .f = ~ str_detect(model_formula, .x))]

chr_model_vars <- model_vars[map_lgl(.x = model_vars,
                                     .f = ~ is.character(sample_metadata[[.x]]))]

sample_metadata <- sample_metadata %>% 
  mutate(across(.cols = all_of(chr_model_vars),
                .fns = as.factor))

# Summarise counts to gene level ------------------------------------------

# Construct txi object
txi_fpd_desq <- tximport(files = sample_metadata$files,
                         type = "salmon",
                         txOut = FALSE,
                         countsFromAbundance = "no",
                         tx2gene = transcript_metadata)

# DESeq2 workflow ---------------------------------------------------------

# Construct dds
dds <- DESeqDataSetFromTximport(txi = txi_fpd_desq,
                                colData = sample_metadata,
                                design = formula(model_formula))

# Subset to protein-coding genes only
if(protein_coding_only) {
  
  protein_coding_genes <- gene_metadata %>%
    filter(gene_type == "protein_coding") %>%
    pull(ensembl_gene_id_version)
  
  dds <- dds[rownames(dds) %in% protein_coding_genes,]
  
}

# Filter out low count genes
min_n_sample <- sample_metadata %>%
  group_by(.data[[coef_var]]) %>%
  add_count() %>%
  pull(n) %>%
  min()
keep <- rowSums(counts(dds) > min_cnt_excl) >= min_n_sample
dds <- dds[keep,]

# Test if DIV has an effect at any timepoint after the reference level time point (DIV0) 
# => Using likelihood ratio test
dds_lrt <- DESeq(object = dds, 
                 test = "LRT",
                 full = formula(model_formula),
                 reduced = formula(reduced_model_formula),
                 parallel = T,
                 BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))

res_lrt <- results(object = dds_lrt,
                   alpha = padj_dge)

res_lrt <- tidy_DESeqResults(input = res_lrt,
                             gene_id = "ensembl_gene_id_version")

res_lrt <- res_lrt %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) %>% 
  filter(gene_symbol == "SNCA") %>% 
  select(gene_symbol, everything())

# Test if specific levels of DIV have an effect compared to the reference level time point (DIV0)
# => Using Wald test
dds_wald <- DESeq(object = dds, 
                  test = "Wald",
                  parallel = T,
                  BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))

resultsNames(dds_wald)

# => Extract PD vs Control results
res_wald_PDvsControl <- results(object = dds_wald,
                                alpha = padj_dge,
                                name = "disease_status_PD_vs_CTRL")

res_wald_PDvsControl <- tidy_DESeqResults(input = res_wald_PDvsControl,
                                          gene_id = "ensembl_gene_id_version")

res_wald_PDvsControl <- res_wald_PDvsControl %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) %>% 
  filter(gene_symbol == "SNCA") %>% 
  select(gene_symbol, everything())

# => Extract DIV25 vs DIV0 results
res_wald_ProdromalvsControl <- results(object = dds_wald,
                                       alpha = padj_dge,
                                       name = "disease_status_Prodromal_vs_CTRL")

res_wald_ProdromalvsControl <- tidy_DESeqResults(input = res_wald_ProdromalvsControl,
                                                 gene_id = "ensembl_gene_id_version")

res_wald_ProdromalvsControl <- res_wald_ProdromalvsControl %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) %>% 
  filter(gene_symbol == "SNCA") %>% 
  select(gene_symbol, everything())

# Extract SNCA TPM values -------------------------------------------------

# TPMs are counts normalised for library size and transcript length but not for library composition
# => Can be used for visual representation but not for statistical inference

# Extract TPM
tpm <- txi_fpd_desq$abundance

# Retrieve SNCA Ensembl gene ID
ensembl_gene_id_snca <- gene_metadata %>% 
  filter(gene_symbol == "SNCA") %>% 
  pull(ensembl_gene_id_version)

# Subset SNCA TPMs
tpm_snca <- tpm[ensembl_gene_id_snca,]

tpm_snca <- tibble(sample_spe_dir = names(tpm_snca),
                   tpm_snca = unname(tpm_snca))

# Tidy SNCA TPM table
tpm_snca <- tpm_snca %>% 
  left_join(sample_metadata, 
            by = join_by(sample_spe_dir)) %>% 
  select(patno, div, disease_status, gender, 
         age_at_diagnosis, age_at_sample_collection, tpm_snca) %>% 
  mutate(log2_tpm1p_snca = log(x = tpm_snca+1L, base = 2))

write_xlsx(tpm_snca, file.path(out_dir_path, "disease_statusD65_snca_output_data.xlsx"))

# Extract SNCA DESeq2 normalized counts -----------------------------------

# DESeq2 normalised counts are counts normalised for library size and library composition but not for transcripts length
# => Can be used for visual representation but only for that of the same gene across multiple conditions

# Extract DESeq2 normlised counts
norm_cnt <- counts(dds_lrt, normalized = TRUE)

# Subset SNCA normalised counts
norm_cnt_snca <- norm_cnt[ensembl_gene_id_snca,]

norm_cnt_snca <- tibble(sample_spe_dir = names(norm_cnt_snca),
                        norm_cnt_snca = unname(norm_cnt_snca))

# Tidy SNCA normalised counts
norm_cnt_snca <- norm_cnt_snca %>% 
  left_join(sample_metadata, 
            by = join_by(sample_spe_dir)) %>% 
  select(patno, div, disease_status, gender, 
         age_at_diagnosis, age_at_sample_collection, norm_cnt_snca) 
