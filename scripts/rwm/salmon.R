# File description --------------------------------------------------------

### Goal
# salmon.R computes the relative abundance of SNCA transcripts across bulk RNA-seq samples.

### Execution steps
# x rwm/make_feature_metadata.R
# ===>>> rwm/salmon.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - sample_metadata.csv available on GitHub repository 'XXX', at metadata/rwm/sample_metadata.csv
# - gene_metadata.csv produced by make_feature_metadata.R
# - transcript_metadata.csv produced by make_feature_metadata.R
# - snca_metadata.csv produced by make_feature_metadata.R
# - Salmon's quant.sf files available from the Gene Expression Omnibus (GEO) under accession number GSE353879
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Iona Thomas-Wright (iona.thomaswright@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Thomas-Wright_2026/data/rwm/salmon"

# Path to sample metadata file (csv)
sample_metadata_path <- "publications/Thomas-Wright_2026/metadata/rwm/sample_metadata.csv"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Thomas-Wright_2026/data/rwm/make_feature_metadata/v46/gene_metadata.csv"

# Path to transcript metadata file (csv)
transcript_metadata_path <- "publications/Thomas-Wright_2026/data/rwm/make_feature_metadata/v46/transcript_metadata.csv"

# Path to snca metadata file (csv)
snca_metadata_path <- "publications/Thomas-Wright_2026/data/rwm/make_feature_metadata/v46/snca_metadata.csv"

# Path to directory containing Salmon quant files
salmon_quant_dir_path <- "publications/Thomas-Wright_2026/seq_data_Thomas-Wright_2026/"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Load feature metadata ---------------------------------------------------

# Load gene metadata
gene_metadata <- read_csv(gene_metadata_path)

# Load transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path) 

# Load SNCA metadata
snca_metadata <- read_csv(snca_metadata_path)

# Load sample metadata ----------------------------------------------------

sample_metadata <- read_csv(sample_metadata_path)

# Load Salmon's quantification results ------------------------------------

# Construct paths to Salmon's quant.sf files
salmon_quant_files <- setNames(object = map_chr(.x = sample_metadata$sample_id,
                                                .f = ~ file.path(salmon_quant_dir_path, paste0(.x, "_quant.sf"))),
                               nm = sample_metadata$line_name_generic)

# Construct txi object
txi <- tximport::tximport(files = salmon_quant_files,
                          type = "salmon",
                          txOut = TRUE,
                          tx2gene = transcript_metadata)

# Abundance of SNCA transcripts as log2(TPM+1) ----------------------------

# Construct SNCA TPM table for plotting
tpm_snca <- as_tibble(x = txi$abundance, rownames = "ensembl_transcript_id_version") %>%
  right_join(snca_metadata,
             by = join_by(ensembl_transcript_id_version)) %>%
  pivot_longer(cols = all_of(sample_metadata$line_name_generic),
               names_to = "line_name_generic",
               values_to = "tpm") %>%
  left_join(sample_metadata[, c("line_name", "line_name_generic", "genotype")],
            by = join_by(line_name_generic)) %>%
  mutate(log2_tpm1p = log(x = tpm+1L, base = 2),
         width_peptide = ifelse(is.na(width_peptide), 0, width_peptide)) %>%
  group_by(ensembl_transcript_id_version) %>%
  mutate(mean_log2_tpm1p = mean(log2_tpm1p, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(width_peptide, mean_log2_tpm1p) %>%
  mutate(transcript_id = fct_relevel(transcript_id, unique(.data$transcript_id))) %>%
  mutate(genotype = fct_recode(as_factor(genotype),
                               "Control" = "CTRL",
                               "SNCA-trip" = "TRIP"))

# Heatmaps of log2(TPM+1) for SNCA transcripts split by sample name or genotype
p_tpm_snca <- lapply(X = c("line_name_generic", "genotype", "gene_symbol"),
                     FUN = function(x) {
                       ggplot(data = tpm_snca,
                              mapping = aes(x = eval(sym(x)),
                                            y = transcript_id,
                                            fill = log2_tpm1p)) +
                         geom_tile(color = "white") +
                         scale_fill_gradientn(colours = c("white", "brown3")) +
                         ggpubr::theme_pubr(legend = "right") +
                         theme(axis.text.x = element_text(size = 14, face = "bold"),
                               axis.text.y = element_text(size = 12, face = "bold"),
                               axis.title.y = element_text(size = 10),
                               axis.line = element_line(linewidth = 0.8),
                               legend.title = element_text(size = 12,face = "bold")) +
                         labs(x = NULL,
                              y = NULL,
                              fill = quote(log[2] * "(TPM+1)"))
                     }) %>%
  setNames(nm = c("line_name_generic", "genotype", "gene_symbol"))

p_tpm_snca$gene_symbol +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_tpm_snca$genotype
p_tpm_snca$line_name_generic

# Write plots to file
walk(.x = c("line_name_generic", "genotype", "gene_symbol"),
     .f = ~ ggsave(filename = file.path(out_dir_path, paste0("snca_mrna_log2tpm1p_", .x, ".pdf")),
                   plot = p_tpm_snca[[.x]],
                   height = 7, width = 6))

# Relative abundance of SNCA transcripts as percentage of total -----------

# Extract Salmon's TPM values
tpm <- txi$abundance

# Construct list of SNCA transcripts
snca_transcripts_all <- transcript_metadata %>% 
  filter(gene_symbol == "SNCA")

snca_transcripts <- snca_transcripts_all$ensembl_transcript_id_version

# Rescale SNCA TPMs
resc_snca_tpm <- rescale_tpm(input = tpm, 
                             genes = snca_transcripts)

# Average rescaled SNCA TPMs by genotype
grp_df <- sample_metadata[, c("line_name_generic", "genotype")]

avg_resc_snca_tpm <- average_rescaled_tpm(input = resc_snca_tpm,
                                          grp_df = grp_df,
                                          sample_id = line_name_generic,
                                          gene_id = ensembl_transcript_id_version)

avg_resc_snca_tpm <- avg_resc_snca_tpm %>%
  left_join(unique(tpm_snca[, c("ensembl_transcript_id_version", "transcript_id")]),
            by = join_by(ensembl_transcript_id_version))

# Plot relative abundance of SNCA transcripts as percentage of total
tpm_rescaled_transcriptseq <- plot_rescaled_tpm(input = avg_resc_snca_tpm,
                                                feature_type = genes,
                                                feature_id = transcript_id,
                                                grp_id = genotype,
                                                y_label = "% of SNCA transcripts",
                                                colors = c(
                                                  "#C04343", "#C04344", "#C04345", "#AD486B", "#9A4D92", 
                                                  "#8753B9", "#7358DF", "#615DEA", "#4D62F3", "#3967F9", 
                                                  "#296CFC", "#296FFD", "#2940FF", "#2941FF", "#2942FF", 
                                                  "#2943FF", "#2944FF"))


# Define SNCA transcript sets
# => Each set includes SNCA transcripts with identical CDS
gs_df <- snca_metadata %>%
  select(ensembl_transcript_id_version, width_peptide) %>%
  mutate(width_peptide = ifelse(is.na(width_peptide),
                                "unknown",
                                width_peptide))

# Sum rescaled SNCA TPMs by transcript set and then average by genotype
avg_sum_resc_snca_tpm <- average_sum_rescaled_tpm(input = resc_snca_tpm,
                                                  grp_df = grp_df,
                                                  gs_df = gs_df,
                                                  sample_id = line_name_generic,
                                                  gene_id = ensembl_transcript_id_version)

avg_sum_resc_snca_tpm <- avg_sum_resc_snca_tpm %>%
  mutate(set_id = ifelse(width_peptide == "unknown",
                         "SNCA - unknown",
                         paste0("SNCA - ", width_peptide, "aa"))) %>%
  mutate(genotype = fct_recode(as_factor(genotype),
                               "Control" = "CTRL",
                               "SNCA-trip" = "TRIP"))

# Plot relative abundance of SNCA transcript sets as percentage of total
tpm_rescaled_protseq <- plot_rescaled_tpm(input = avg_sum_resc_snca_tpm,
                                          feature_type = genesets,
                                          feature_id = set_id,
                                          grp_id = genotype,
                                          y_label = "% of SNCA transcript universe",
                                          colors = c(
                                            "#C04343", "#9A4D92", "#8753B9", 
                                            "#615DEA", "#4D62F3", "#3967F9",
                                            "#296CFC", "#296FFD", "#2940FF")) +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face="bold"),
        axis.line = element_line(linewidth = 0.8),
        legend.text = element_text(size = 12,face = "bold"))

ggsave(filename = file.path(out_dir_path, "snca_mrna_prop_genotype.pdf"),
       plot = tpm_rescaled_protseq,
       width = 12,
       height = 10,
       device = cairo_pdf)
