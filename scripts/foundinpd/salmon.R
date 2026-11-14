# File description --------------------------------------------------------

### Goal
# salmon.R computes the relative abundance of SNCA transcripts across bulk RNA-seq samples.

### Execution steps
# x foundinpd/make_sample_metadata.R
# x foundinpd/make_feature_metadata.R
# ===>>> foundinpd/salmon.R
# - foundinpd/dgea.R

# => Any script listed before the current file must be executed beforehand 

### Requirements
# - filtered_sample_metadata.rds produced by make_sample_metadata.R
# - gene_metadata.csv produced by make_feature_metadata.R
# - transcript_metadata.csv produced by make_feature_metadata.R
# - snca_metadata.csv produced by make_feature_metadata.R
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Iona Thomas-Wright (iona.thomaswright@gmail.com)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Thomas-Wright_2026/data/foundinpd/salmon"

# Path to sample metadata file (csv)
sample_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_sample_metadata/filtered_sample_metadata.rds"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_feature_metadata/v29/gene_metadata.csv"

# Path to transcript metadata file (csv)
transcript_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_feature_metadata/v29/transcript_metadata.csv"

# Path to snca metadata file (csv)
snca_metadata_path <- "publications/Thomas-Wright_2026/data/foundinpd/make_feature_metadata/v29/snca_metadata.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)
library(writexl)
library(afex)
library(emmeans)

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

sample_metadata <- read_rds(sample_metadata_path)

# Load Salmon's quantification results ------------------------------------

# Construct txi object
txi <- tximport::tximport(files =  sample_metadata$files,
                          type = "salmon",
                          txOut = TRUE,
                          tx2gene = transcript_metadata)

# Exp1 Abundance of SNCA transcripts as log2(TPM+1) -----------------------

# Construct SNCA TPM tables for plotting (For all time points and DIV65 only)
tpm_snca <- as_tibble(x = txi$abundance, rownames = "ensembl_transcript_id_version") %>%
  right_join(snca_metadata,
             by = join_by(ensembl_transcript_id_version)) %>%
  pivot_longer(cols = all_of(sample_metadata$sample_spe_dir),
               names_to = "sample_spe_dir",
               values_to = "tpm") %>%
  left_join(sample_metadata[, c("sample_spe_dir", "cdi_number", "div", "disease_status")],
            by = join_by(sample_spe_dir)) %>%
  mutate(log2_tpm1p = log(x = tpm+1L, base = 2),
         width_peptide = ifelse(is.na(width_peptide), 0, width_peptide)) %>%
  group_by(ensembl_transcript_id_version) %>%
  mutate(mean_log2_tpm1p = mean(log2_tpm1p, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(width_peptide, mean_log2_tpm1p) %>%
  mutate(transcript_id = fct_relevel(transcript_id, unique(.data$transcript_id)))%>%
  mutate(disease_status = fct_recode(as_factor(disease_status),
                                     "Control" = "CTRL",
                                     "Prodromal" = "Prodromal",
                                     "PD" = "PD"))

tpm_snca_65 <- tpm_snca %>%
  filter(div == "DIV65")

# Exp1.1 Heatmaps of log2(TPM+1) for SNCA transcripts 
tpm_snca_tables <- list(all = tpm_snca,
                        div65 = tpm_snca_65)

vars <- c("div", "disease_status", "gene_symbol")

p_tpm_snca <- unlist(
  lapply(names(tpm_snca_tables), function(dataset_name) {
    
    df <- tpm_snca_tables[[dataset_name]]
    lapply(vars, function(var_name) {
      ggplot(data = df,
             mapping = aes(x = !!sym(var_name),
                           y = transcript_id,
                           fill = log2_tpm1p)) +
        geom_tile(color = "white") +
        scale_fill_gradientn(colours = c("white", "brown3"),
                             limits = c(0, 3),
                             oob = scales::squish) +
        ggpubr::theme_pubr(legend = "right") +
        theme(axis.text.x = element_text(size = 12, face = "bold"),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 10),
              axis.line = element_line(linewidth = 0.8),
              legend.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 12)) +
        ggpubr::rotate_x_text(angle = 0) +
        labs(x = NULL,
             y = NULL,
             fill = quote(log[2] * "(TPM+1)"))
    }) %>%
      setNames(paste0(vars, "_", dataset_name))
  }),
  recursive = FALSE
)

logTPM_pooledStatus_byDIV <- p_tpm_snca$div_all
logTPM_pooledStatus_byDIV
logTPM_pooledDIV_byStatus <- p_tpm_snca$disease_status_all
logTPM_pooledDIV_byStatus
logTPM_DIV65_byStatus <- p_tpm_snca$disease_status_div65
logTPM_DIV65_byStatus
logTPM_pooledDIV_pooledStatus <- p_tpm_snca$gene_symbol_all
logTPM_DIV65_pooledStatus <-p_tpm_snca$gene_symbol_div65

walk(.x = c("logTPM_pooledStatus_byDIV",
            "logTPM_pooledDIV_byStatus", 
            "logTPM_DIV65_byStatus", 
            "logTPM_pooledDIV_pooledStatus", 
            "logTPM_DIV65_pooledStatus"),
     .f = ~ ggsave(filename = file.path(out_dir_path, paste0("snca_mrna_", .x, ".pdf")),
                   plot = p_tpm_snca[[.x]],
                   height = 7, width = 6))

# Exp2 Relative abundance of SNCA transcripts as percentage of total --------


# Extract Salmon's TPM values
tpm <- txi$abundance

# Construct list of SNCA transcripts
snca_transcripts_all <- transcript_metadata %>% 
  filter(gene_symbol == "SNCA")

snca_transcripts <- snca_transcripts_all$ensembl_transcript_id_version

# Rescale SNCA TPMs
resc_snca_tpm <- rescale_tpm(input = tpm, 
                             genes = snca_transcripts) 

# Exp2.1 Relative abundance of SNCA transcripts including untranslated regions --------

# Average re-scaled SNCA TPMs 
grp_df <- sample_metadata[, c("sample_spe_dir", "disease_status","div")]

avg_resc_snca_tpm <- average_rescaled_tpm(input = resc_snca_tpm,
                                          grp_df = grp_df,
                                          sample_id = sample_spe_dir,
                                          gene_id = ensembl_transcript_id_version)

avg_resc_snca_tpm <- avg_resc_snca_tpm %>%
  left_join(unique(tpm_snca[, c("ensembl_transcript_id_version", "transcript_id")]),
            by = join_by(ensembl_transcript_id_version))

# Exp2.1.1 - Relative abundance DIV65 only, by Disease Status          

# Filter for DIV65 only
avg_resc_snca_tpm_65 <- avg_resc_snca_tpm  %>%
  filter(.data[["disease_status_div"]] %in% c("CTRL_DIV65","Prodromal_DIV65","PD_DIV65")) %>%
  rename_with(~ "disease_status", .cols = all_of("disease_status_div")) %>%
  mutate(PD_status = case_when(disease_status == "CTRL_DIV65" ~ "Control",
                               disease_status == "Prodromal_DIV65" ~ "Prodromal", 
                               disease_status == "PD_DIV65" ~ "PD"),
         PD_status = fct_relevel(PD_status, "Control", "Prodromal", "PD"))

# Plot as stacked bar graph out of 100%
tpm_rescaled_transcriptseq_65 <- plot_rescaled_tpm(input = avg_resc_snca_tpm_65,
                                                   feature_type = genes,
                                                   feature_id = transcript_id,
                                                   grp_id = PD_status,
                                                   y_label = "% of SNCA transcripts",
                                                   colors = c(
                                                     "#C04343", "#C04344", "#C04345", "#AD486B", "#9A4D92", 
                                                     "#8753B9", "#7358DF", "#615DEA", "#4D62F3", "#3967F9", 
                                                     "#296CFC", "#296FFD", "#2940FF", "#2941FF", "#2942FF", 
                                                     "#2943FF", "#2944FF")) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.line = element_line(linewidth = 0.8),
        legend.text = element_text(size = 10,face = "bold"))

ggsave(filename = file.path(out_dir_path, "snca_mrna_proportion_disease_status_D65.pdf"),
       plot = tpm_rescaled_transcriptseq_65,
       width = 12,
       height = 10,
       device = cairo_pdf)

# Exp2.1.2 - Relative abundance pooled across disease statuses, by DIV
# Average re-scaled SNCA TPM not pulling disease status
grp_df_div <- sample_metadata[, c("sample_spe_dir", "div")]

avg_resc_snca_tpm_div <- average_rescaled_tpm(input = resc_snca_tpm,
                                              grp_df = grp_df_div,
                                              sample_id = sample_spe_dir,
                                              gene_id = ensembl_transcript_id_version)

avg_resc_snca_tpm_div <- avg_resc_snca_tpm_div %>%
  left_join(unique(tpm_snca[, c("ensembl_transcript_id_version", "transcript_id")]),
            by = join_by(ensembl_transcript_id_version))

# Plot as stacked bar graph out of 100%
tpm_rescaled_transcriptseq_div <- plot_rescaled_tpm(input = avg_resc_snca_tpm_div,
                                                    feature_type = genes,
                                                    feature_id = transcript_id,
                                                    grp_id = div,
                                                    y_label = "% of SNCA transcripts",
                                                    colors = c(
                                                      "#C04343", "#C04344", "#C04345", "#AD486B", "#9A4D92", 
                                                      "#8753B9", "#7358DF", "#615DEA", "#4D62F3", "#3967F9", 
                                                      "#296CFC", "#296FFD", "#2940FF", "#2941FF", "#2942FF", 
                                                      "#2943FF", "#2944FF")) +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.line = element_line(linewidth = 0.8),
        legend.text = element_text(size = 10,face = "bold"))

ggsave(filename = file.path(out_dir_path, "snca_mrna_proportion_DIV.pdf"),
       plot = tpm_rescaled_transcriptseq_div,
       width = 12,
       height = 10,
       device = cairo_pdf)

# Exp2.2 Relative abundance of SNCA transcripts that code for the same protein sequence --------

# Define SNCA transcript sets and sum SNCA TPM data
# => Each set includes SNCA transcripts with identical CDS
gs_df <- snca_metadata %>%
  select(ensembl_transcript_id_version, width_peptide) %>%
  mutate(width_peptide = as.character(ifelse(is.na(width_peptide),
                                             "unknown",
                                             width_peptide)))

avg_sum_resc_snca_tpm <- average_sum_rescaled_tpm(input = resc_snca_tpm,
                                                  grp_df = grp_df,
                                                  gs_df = gs_df,
                                                  sample_id = sample_spe_dir,
                                                  gene_id = ensembl_transcript_id_version)

# Exp2.2.1 Relative abundance of coding sequences at D65, by disease status
# Filter for DIV65 only
avg_sum_resc_snca_tpm_65 <- avg_sum_resc_snca_tpm %>%
  mutate(set_id = ifelse(width_peptide == "unknown",
                         "SNCA - unknown",
                         paste0("SNCA - ", width_peptide, "aa"))) %>%
  filter(.data[["disease_status_div"]] %in% c("CTRL_DIV65","Prodromal_DIV65","PD_DIV65")) %>%
  rename_with(~ "disease_status", .cols = all_of("disease_status_div")) %>%
  mutate(PD_status = case_when(disease_status == "CTRL_DIV65" ~ "Control",
                               disease_status == "Prodromal_DIV65" ~ "Prodromal", 
                               disease_status == "PD_DIV65" ~ "PD"),
         PD_status = fct_relevel(PD_status, "Control", "Prodromal", "PD"))

# Plot relative abundance of SNCA coding sequences out of 100%  
tpm_rescaled_protseq_65 <- plot_rescaled_tpm(input = avg_sum_resc_snca_tpm_65,
                                             feature_type = genesets,
                                             feature_id = set_id,
                                             grp_id = PD_status,
                                             y_label = "% of SNCA transcript universe",
                                             colors = c(
                                               "#C04343", "#9A4D92", "#8753B9", 
                                               "#615DEA", "#4D62F3", "#3967F9",
                                               "#296CFC", "#296FFD", "#2940FF")) +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face="bold"),
        axis.line = element_line(linewidth = 0.8),
        legend.text = element_text(size = 10,face = "bold"))

ggsave(filename = file.path(out_dir_path, "snca_mrna_codingSeq_disease_status_D65.pdf"),
       plot = tpm_rescaled_protseq_65,
       width = 8,
       height = 8,
       device = cairo_pdf)

# Exp2.2.2 Relative abundance of coding sequences pooled across disease status, by DIV
# Average and sum re-scaled SNCA TMPs across disease statuses
avg_sum_resc_snca_tpm_div <- average_sum_rescaled_tpm(input = resc_snca_tpm,
                                                      grp_df = grp_df_div,
                                                      gs_df = gs_df,
                                                      sample_id = sample_spe_dir,
                                                      gene_id = ensembl_transcript_id_version)

avg_sum_resc_snca_tpm_div <- avg_sum_resc_snca_tpm_div %>%
  mutate(set_id = ifelse(width_peptide == "unknown",
                         "SNCA - unknown",
                         paste0("SNCA - ", width_peptide, "aa")))

# Plot relative abundance of SNCA coding sequences out of 100%  
tpm_rescaled_protseq_div <- plot_rescaled_tpm(input = avg_sum_resc_snca_tpm_div,
                                              feature_type = genesets,
                                              feature_id = set_id,
                                              grp_id = div,
                                              y_label = "% of SNCA transcript universe",
                                              colors = c(
                                                "#C04343", "#9A4D92", "#8753B9", 
                                                "#615DEA", "#4D62F3", "#3967F9",
                                                "#296CFC", "#296FFD", "#2940FF")) +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.line = element_line(linewidth = 0.8),
        legend.text = element_text(size = 10,face = "bold"))

ggsave(filename = file.path(out_dir_path, "snca_mrna_codingSeq_DIV.pdf"),
       plot = tpm_rescaled_protseq_div,
       width = 8,
       height = 8,
       device = cairo_pdf)

# Exp2.3 Plot relative abundance of SNCA coding sequences as a percentage of total for DIVs by separate disease statuses --------

# Set up the resc_snca_tpm tibble so I can subset from it
resc_snca_tpm_sampleData <- as_tibble(resc_snca_tpm, rownames = NA) %>% 
  rownames_to_column(var = "ensembl_transcript_id_version") %>%
  transpose_cnt_df(gene_id = "ensembl_transcript_id_version") %>%
  rename(sample_spe_dir = barcode) %>%
  left_join(sample_metadata[, c("sample_spe_dir", "div", "disease_status")]) 

# Cycle through a loop to: filter by disease status, average tpm values, plot them and save plots
Disease <- c("CTRL", "Prodromal", "PD")

for (i in Disease) {
  resc_snca_tpm_loop <- resc_snca_tpm_sampleData %>%
    filter(disease_status == i) %>%
    select(-div, -disease_status) %>%
    column_to_rownames("sample_spe_dir") 
  
  resc_snca_tpm_loop <- as.matrix(resc_snca_tpm_loop) %>%
    t()
  
  avg_sum_resc_snca_tpm_loop <- average_sum_rescaled_tpm(input = resc_snca_tpm_loop,
                                                         grp_df = grp_df_div,
                                                         gs_df = gs_df,
                                                         sample_id = sample_spe_dir,
                                                         gene_id = ensembl_transcript_id_version)
  
  avg_sum_resc_snca_tpm_loop <- avg_sum_resc_snca_tpm_loop  %>%
    mutate(set_id = ifelse(width_peptide == "unknown",
                           "SNCA - unknown",
                           paste0("SNCA - ", width_peptide, "aa")))
  
  tpm_rescaled_protseq_div <- plot_rescaled_tpm(input = avg_sum_resc_snca_tpm_loop,
                                                feature_type = genesets,
                                                feature_id = set_id,
                                                grp_id = div,
                                                y_label = "% of SNCA transcript universe",
                                                colors = c(
                                                  "#C04343", "#9A4D92", "#8753B9", 
                                                  "#615DEA", "#4D62F3", "#3967F9",
                                                  "#296CFC", "#296FFD", "#2940FF")) +
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face="bold"),
          axis.line = element_line(linewidth = 0.8),
          legend.text = element_text(size = 10,face = "bold"))
  
  ggsave(filename = file.path(out_dir_path, paste0("snca_mrna_codingSeq_", i, "_only.pdf")),
         plot = tpm_rescaled_protseq_div,
         width = 8,
         height = 8,
         device = cairo_pdf)
  
  assign(paste0("resc_sum_snca_tpm_", i), avg_sum_resc_snca_tpm_loop)
}

# Exp3 One and Two-way Anovas to look at the effect of DIV on 140aa transcript abundance --------

# Summing all the 140aa transcripts and making a DIV-independant sample ID column
resc_snca_tpm_notavg <- as_tibble(resc_snca_tpm, rownames = NA) %>% 
  rownames_to_column(var = "ensembl_transcript_id_version") %>% 
  transpose_cnt_df(gene_id = "ensembl_transcript_id_version") %>%
  rename(sample_spe_dir = barcode) %>%
  select(sample_spe_dir, 
         all_of(snca_metadata %>% 
                  filter(width_peptide == "140") %>% 
                  pull(ensembl_transcript_id_version)))

resc_snca_tpm_notavg_sum <- resc_snca_tpm_notavg %>%
  column_to_rownames(var = "sample_spe_dir") 

resc_snca_tpm_notavg_sum <- rowSums(resc_snca_tpm_notavg_sum) 

resc_snca_tpm_notavg_sum <- tibble(sample_spe_dir = names(resc_snca_tpm_notavg_sum),
                                   sum_tpm = unname(resc_snca_tpm_notavg_sum))

resc_snca_tpm_140sum <- left_join(resc_snca_tpm_notavg_sum, sample_metadata[, c("sample_spe_dir", "div", "disease_status")], 
                                  by = join_by(sample_spe_dir)) %>%
  mutate(four_digit_id = str_extract(sample_spe_dir, "(?<=_)(\\d{4})(?=_)"))

write_xlsx(resc_snca_tpm_140sum, file.path(out_dir_path, "resc_tpm_all_data.xlsx"))

# Exp3.1 One-way ANOVA for all disease statuses pooled
anova_result <- aov_ez(id = "four_digit_id",   
                       dv = "sum_tpm",           
                       within = "div",        
                       data = resc_snca_tpm_140sum)
print(anova_result)
emmeans(anova_result, pairwise ~ div, adjust = "bonferroni")

# Exp3.2 One-way ANOVA for disease status separately

resc_snca_tpm_140sum_disease <- resc_snca_tpm_140sum %>% 
  pivot_wider(names_from = disease_status, values_from = sum_tpm)

Disease <- c("CTRL", "Prodromal", "PD")

for (i in Disease) {
  anova_result_i <- aov_ez( id = "four_digit_id",   
                            dv = i,           
                            within = "div",        
                            data = resc_snca_tpm_140sum_disease)
  print(anova_result_i)
  emmeans_result_i <- emmeans(anova_result_i, pairwise ~ div, adjust = "bonferroni")
  print(emmeans_result_i)
}

# Exp3.2 Two-way ANOVA for an effect of disease status and DIV together (with interaction)
twoway_anova_result  <- aov_ez(id = "four_digit_id",   
                               dv = "sum_tpm",             
                               between = "disease_status",     
                               within = "div",           
                               data = resc_snca_tpm_140sum)

print(twoway_anova_result)
print(emmeans(twoway_anova_result, pairwise ~ div, adjust = "bonferroni"))
print(emmeans(twoway_anova_result, pairwise ~ disease_status, adjust = "bonferroni"))
print(emmeans(twoway_anova_result, pairwise ~ disease_status| div, adjust = "bonferroni"))

