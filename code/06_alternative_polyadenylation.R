# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(data.table)
library(biomaRt)

# Set file paths
conservation_file_path <- "data/mirna_variants_ocs.csv"
#MicroT database for microRNA-Target interactions within 3'UTRs and mre score of >0.7.
targets_file <- ""
# APA Database with filename atlas.clusters.2.0.GRCh38.96.tsv
apa_file_path <- ""

# Load and process data
targets <- fread(targets_file, sep = "\t")
conservation <- fread(conservation_file_path) %>%
  filter(type_mirna == "miRNA") %>%
  mutate(base_name = gsub("mir", "miR", gsub("-3p|-5p", "", mirna_name)),
         OCS = as.numeric(OCS))

# Define top/bottom 5% miRNAs by OCS
ocs_thresholds <- quantile(conservation$OCS, probs = c(0.05, 0.95), na.rm = TRUE)
top_miRNAs <- filter(conservation, OCS >= ocs_thresholds[2])$mirna_name
bottom_miRNAs <- filter(conservation, OCS <= ocs_thresholds[1])$mirna_name

# Filter targets
targets <- targets %>% filter(mirna %in% c(top_miRNAs, bottom_miRNAs))
top_targets <- filter(targets, mirna %in% top_miRNAs)
bottom_targets <- filter(targets, mirna %in% bottom_miRNAs)

# Load APA data and prepare one-APA and multi-APA gene groups
apa <- fread(apa_file_path) %>%
  filter(annotation == "TE") %>%
  select(1:13) %>%
  distinct() %>%
  drop_na()

genes_one_apa <- apa %>% group_by(gene_name) %>% filter(n() == 1) %>% ungroup()
genes_multi_apa <- apa %>%
  group_by(gene_name) %>%
  filter(n() > 1) %>%
  slice_max(frac_samples, n = 1) %>%
  ungroup()

# Gene symbol mapping
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(targets$ensembl_gene_id),
                  mart = mart)

targets <- left_join(targets, gene_map, by = "ensembl_gene_id") %>%
  filter(hgnc_symbol != "")

# Focus on targets in genes with one APA site
targets_one <- filter(targets, hgnc_symbol %in% genes_one_apa$gene_name)

# Unique genes bound by miRNAs, stratified by APA group
gene_counts <- data.frame()
for (thr in seq(0, 0.1, 0.005)) {
  tf <- filter(targets, mre_score >= thr) %>%
    left_join(apa, by = c("hgnc_symbol" = "gene_name"))
  
  one <- filter(tf, hgnc_symbol %in% genes_one_apa$gene_name)
  multi <- filter(tf, hgnc_symbol %in% genes_multi_apa$gene_name)
  
  row <- data.frame(
    mre_score_threshold = thr,
    top5_one_apa = n_distinct(filter(one, mirna %in% top_miRNAs)$hgnc_symbol),
    bottom5_one_apa = n_distinct(filter(one, mirna %in% bottom_miRNAs)$hgnc_symbol),
    top5_multi_apa = n_distinct(filter(multi, mirna %in% top_miRNAs)$hgnc_symbol),
    bottom5_multi_apa = n_distinct(filter(multi, mirna %in% bottom_miRNAs)$hgnc_symbol)
  )
  gene_counts <- rbind(gene_counts, row)
}
#Plot D
ggplot(gene_counts, aes(x = mre_score_threshold)) +
  geom_line(aes(y = top5_one_apa, color = "Top 5% - One APA")) +
  geom_line(aes(y = bottom5_one_apa, color = "Bottom 5% - One APA")) +
  geom_line(aes(y = top5_multi_apa, color = "Top 5% - Multiple APA")) +
  geom_line(aes(y = bottom5_multi_apa, color = "Bottom 5% - Multiple APA")) +
  labs(title = "MRE Distribution in APA Genes",
       x = "MRE Score", y = "Gene Count") +
  theme_minimal() +
  theme(legend.position = "top")

# Plot E
#Ratio of genes bound before APA site
binding_ratios <- data.frame()
for (thr in seq(0, 0.06, 0.005)) {
  tf <- filter(targets_one, mre_score >= thr) %>%
    left_join(genes_one_apa, by = c("hgnc_symbol" = "gene_name"))
  
  total_genes <- n_distinct(tf$hgnc_symbol)
  top_before <- n_distinct(filter(tf, mirna %in% top_miRNAs, start < rep)$hgnc_symbol)
  bottom_before <- n_distinct(filter(tf, mirna %in% bottom_miRNAs, start < rep)$hgnc_symbol)
  
  row <- data.frame(
    mre_score_threshold = thr,
    top_ratio = ifelse(total_genes > 0, top_before / total_genes, 0),
    bottom_ratio = ifelse(total_genes > 0, bottom_before / total_genes, 0),
    total = total_genes
  )
  binding_ratios <- rbind(binding_ratios, row)
}

ggplot(binding_ratios, aes(x = mre_score_threshold)) +
  geom_line(aes(y = top_ratio, color = "Top 5%")) +
  geom_line(aes(y = bottom_ratio, color = "Bottom 5%")) +
  geom_line(aes(y = total / max(total), color = "Scaled Total Genes")) +
  labs(title = "MRE relative to APA",
       x = "MRE Score", y = "Before Apa / Total") +
  theme_minimal() +
  theme(legend.position = "top")

