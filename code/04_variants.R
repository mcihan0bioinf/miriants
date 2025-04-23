#Load libraries
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(readr)
library(data.table)
library(patchwork)
library(tidyr)
library(pheatmap)

data <- fread("data/mirna_variants_ocs.csv")

#Filter mature miRNAs
miRNA_data <- data %>% filter(type_mirna == "miRNA") %>%
  mutate(position_relative = ifelse(Strand_mirna == "+", Start - Start_mirna + 1, End_mirna - Start + 1)) %>%
  mutate(across(starts_with("AF_"), ~ replace(., . == ".", 0) %>% as.numeric()))

#Plot A
#Set thresholds
snp_threshold <- 0.01
high_af_threshold <- 0.05
#Categorize SNP variants
miRNA_data$variant_type <- case_when(
  miRNA_data$Allele_freq >= high_af_threshold ~ "High AF (≥0.05)",
  miRNA_data$Allele_freq >= snp_threshold ~ "SNP (0.01 - 0.049)",
  TRUE ~ "Non-SNP (<0.01)"
)
plot_a <- ggplot(miRNA_data, aes(x = as.numeric(OCS), y = Allele_freq, color = variant_type)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("darkgreen", "red", "black"),
                     labels = c("Common Variants", "Rare Variants", "SNPs"),
                     name = "Variant Type") +
  geom_vline(xintercept = c(0.923, 0.702), linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = c(snp_threshold, high_af_threshold), linetype = "dashed", size = 0.3) +
  labs(title = "OCS vs. Allele Frequency for miRNA", x = "OCS", y = "Allele Frequency") +
  theme_minimal()

#Plot B
#Filter indels
indel_data <- miRNA_data %>% filter(Type == "indel")
#Categorize indel variants
indel_data$variant_type <- case_when(
  indel_data$Allele_freq >= high_af_threshold ~ "High AF (≥0.05)",
  indel_data$Allele_freq >= snp_threshold ~ "SNP (0.01 - 0.049)",
  TRUE ~ "Non-SNP (<0.01)"
)

plot_b <- ggplot(indel_data, aes(x = as.numeric(OCS), y = Allele_freq, color = variant_type)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("darkgreen", "red", "black"),
                     labels = c("Common Variants", "Rare Variants", "SNPs"),
                     name = "Variant Type") +
  geom_vline(xintercept = c(0.923, 0.702), linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = c(snp_threshold, high_af_threshold), linetype = "dashed", size = 0.3) +
  labs(title = "OCS vs. Allele Frequency for Indels", x = "OCS", y = "Allele Frequency") +
  theme_minimal()

#Plot C + D

#Group all miRNAs
overall_position_counts_all <- miRNA_data %>%
  mutate(group = "All miRNAs") %>%
  group_by(position_relative, group) %>%
  summarise(count = n(), .groups = "drop")

#Group top and bottom 5%
top_bottom_data <- miRNA_data %>%
  filter(Rank %in% c("Top", "Bottom")) %>%
  mutate(group = ifelse(Rank == "Top", "Top 5% OCS", "Bottom 5% OCS"))

overall_position_counts_top_bottom <- top_bottom_data %>%
  group_by(position_relative, group) %>%
  summarise(count = n(), .groups = "drop")

#Compute miRNA lengths
mirna_lengths <- miRNA_data %>%
  distinct(mirna_name, .keep_all = TRUE) %>%
  mutate(mirna_length = End_mirna - Start_mirna + 1) %>%
  select(mirna_name, mirna_length, Rank)

#Coverage function
compute_coverage <- function(df) {
  max_pos <- max(df$mirna_length, na.rm = TRUE)
  data.frame(position_relative = 1:max_pos,
             mirna_coverage = sapply(1:max_pos, function(p) sum(df$mirna_length >= p) / nrow(df)))
}

#Coverage by group
cov_all <- compute_coverage(mirna_lengths)
cov_top <- compute_coverage(filter(mirna_lengths, Rank == "Top"))
cov_bottom <- compute_coverage(filter(mirna_lengths, Rank == "Bottom"))

#Plot all grouped
plot_c <- ggplot() +
  geom_col(data = overall_position_counts_all, aes(x = position_relative, y = count), fill = "darkgrey", alpha = 0.8) +
  geom_line(data = cov_all, aes(x = position_relative, y = mirna_coverage * max(overall_position_counts_all$count)),
            color = "black", linetype = "dashed", size = 1) +
  scale_y_continuous(name = "Variant Count",
                     sec.axis = sec_axis(~ . / max(overall_position_counts_all$count), name = "miRNA coverage")) +
  labs(title = "All miRNAs", x = "Relative Position") +
  theme_minimal()

#Plot top vs bottom 5%
plot_d <- ggplot() +
  geom_col(data = overall_position_counts_top_bottom, aes(x = position_relative, y = count, fill = group), position = "dodge", alpha = 0.8) +
  geom_line(data = cov_top, aes(x = position_relative, y = mirna_coverage * max(overall_position_counts_top_bottom$count)),
            color = "darkblue", linetype = "dashed", size = 1) +
  geom_line(data = cov_bottom, aes(x = position_relative, y = mirna_coverage * max(overall_position_counts_top_bottom$count)),
            color = "darkred", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("darkblue", "darkred")) +
  scale_y_continuous(name = "Variant Count",
                     sec.axis = sec_axis(~ . / max(overall_position_counts_top_bottom$count), name = "miRNA coverage")) +
  labs(title = "Top vs Bottom 5% OCS", x = "Relative Position", fill = "Group") +
  theme_minimal()

#Plot E
#Define mutation sets
transitions <- c("A>G", "G>A", "C>T", "T>C")
transversions <- c("A>C", "A>T", "G>C", "G>T", "C>A", "C>G", "T>A", "T>G")

#Initialize results
results <- data.frame()

#Function to count mutation types by rank
count_by_rank <- function(data, mutation_set, rank_label) {
  sum(data$mutation[data$Rank == rank_label] %in% mutation_set, na.rm = TRUE)
}

#Iterate over AF thresholds
for (af_thresh in seq(0, 0.01, by = 0.0001)) {
  
  filtered_data <- data %>%
    filter(Allele_freq > af_thresh) %>%
    mutate(mutation = paste(Reference, Alternative, sep = ">"))
  
  regions <- list(
    "Overall" = filtered_data,
    "Seed" = filtered_data %>% filter(Seed == "Yes"),
    "Non-Seed" = filtered_data %>% filter(Seed == "No"),
    "Precursor" = filtered_data %>% filter(Seed == "No_primary")
  )
  
  for (region_name in names(regions)) {
    region_data <- regions[[region_name]]
    
    total_trans <- sum(region_data$mutation %in% transitions, na.rm = TRUE)
    total_transv <- sum(region_data$mutation %in% transversions, na.rm = TRUE)
    ti_tv_ratio <- ifelse(total_transv > 0, total_trans / total_transv, NA)
    
    rank_tab <- table(region_data$Rank)
    get_count <- function(rank) ifelse(rank %in% names(rank_tab), rank_tab[rank], 0)
    
    results <- rbind(results, data.frame(
      Allele_Freq = af_thresh,
      Region = region_name,
      Ti_Tv_Ratio = ti_tv_ratio,
      Transitions = total_trans,
      Transversions = total_transv,
      Bottom = get_count("Bottom"),
      Middle = get_count("Middle"),
      None = get_count("None"),
      Top = get_count("Top"),
      Bottom_Transitions = count_by_rank(region_data, transitions, "Bottom"),
      Middle_Transitions = count_by_rank(region_data, transitions, "Middle"),
      None_Transitions = count_by_rank(region_data, transitions, "None"),
      Top_Transitions = count_by_rank(region_data, transitions, "Top"),
      Bottom_Transversions = count_by_rank(region_data, transversions, "Bottom"),
      Middle_Transversions = count_by_rank(region_data, transversions, "Middle"),
      None_Transversions = count_by_rank(region_data, transversions, "None"),
      Top_Transversions = count_by_rank(region_data, transversions, "Top")
    ))
  }
}

#Label points: every 10th + AF = 0
label_points <- results %>%
  group_by(Region) %>%
  filter(Allele_Freq == 0 | row_number() %% 10 == 0)

plot_e <- ggplot(results, aes(x = Allele_Freq, y = Ti_Tv_Ratio, color = Region)) +
  geom_line(size = 1) +
  geom_point(data = label_points, size = 2) +
  geom_text(data = label_points, aes(label = Transitions), vjust = -0.5, size = 3.5, fontface = "bold") +
  labs(
    title = "Transition/Transversion Ratio Across miRNA Regions",
    x = "Allele Frequency",
    y = "Transition/Transversion Ratio"
  ) +
  scale_color_manual(
    name = "Region",
    values = c("Overall" = "black", "Seed" = "blue", "Non-Seed" = "red", "Precursor" = "darkgreen")
  ) +
  theme_minimal()

#Plot F
# Define weights for each component
SCS_weight <- 0.2  
NSCS_weight <- 0.1  
PCS_weight <- 0.6  
TVS_weight <- 0.1  

# Extract population columns
population_columns <- grep("^AF_joint", colnames(miRNA_data), value = TRUE)
population_columns=c("AF_joint_afr", "AF_joint_amr", "AF_joint_asj", 
                     "AF_joint_eas", "AF_joint_fin", "AF_joint_nfe", 
                     "AF_joint_mid", "AF_joint_sas", "AF_joint_remaining")
# Compute OCS for each population
ocs_results <- list()
for (pop in population_columns) {
  pop_miRNA_data <- miRNA_data %>%
    filter(get(pop) > 0) %>%
    mutate(Allele_freq = get(pop))
  population_ocs <- pop_miRNA_data %>%
    group_by(mirna_name) %>%
    summarise(
      SCS = ifelse(sum(Seed == "Yes") == 0, 1, 1 - sum(Allele_freq[Seed == "Yes"], na.rm = TRUE) / sum(Seed == "Yes")),
      NSCS = ifelse(sum(Seed == "No") == 0, 1, 1 - sum(Allele_freq[Seed == "No"], na.rm = TRUE) / sum(Seed == "No")),
      PCS = 1 - (n_distinct(Start) / (max(End_mirna) - min(Start_mirna) + 1)),
      TVS = 1 / (1 + n()),
      OCS = (SCS * SCS_weight) + (NSCS * NSCS_weight) + (PCS * PCS_weight) + (TVS * TVS_weight)
    ) %>%
    mutate(Population = pop)
  ocs_results[[pop]] <- population_ocs
}
final_ocs_results <- bind_rows(ocs_results)

# Prepare heatmap data
heatmap_miRNA_data <- final_ocs_results %>%
  select(mirna_name, Population, OCS) %>%
  complete(mirna_name, Population, fill = list(OCS = 1))

# Add total OCS values from original data
mirna_ocs <- miRNA_data %>%
  select(mirna_name, OCS) %>%
  distinct() %>%
  mutate(Population = "Allele_frequency")

heatmap_miRNA_data <- heatmap_miRNA_data %>%
  spread(key = Population, value = OCS) %>%
  replace(is.na(.), 1)

head(heatmap_miRNA_data)

#Filter most variable miRNAs
heatmap_miRNA_data$variance <- apply(heatmap_miRNA_data[, -1], 1, var)
top_20_miRNAs <- heatmap_miRNA_data %>%
  arrange(desc(variance)) %>%
  slice(1:20) %>%
  select(-variance)
heatmap_matrix <- as.matrix(top_20_miRNAs[, -1])
rownames(heatmap_matrix) <- top_20_miRNAs$mirna_name
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
heatmap_matrix <- t(heatmap_matrix)
colnames(heatmap_matrix) <- top_20_miRNAs$mirna_name

ocs_annotation <- miRNA_data %>%
  select(mirna_name, OCS) %>%
  distinct() %>%
  filter(mirna_name %in% colnames(heatmap_matrix))  # Keep only miRNAs in the heatmap
ocs_annotation$OCS <- as.numeric(ocs_annotation$OCS)
annotation_col <- data.frame(OCS = ocs_annotation$OCS)
rownames(annotation_col) <- ocs_annotation$mirna_name
heatmap_matrix=heatmap_matrix[population_columns,]
heatmap_matrix

pheatmap(
  heatmap_matrix,
  clustering_distance_rows = "euclidean",  # Cluster only miRNAs
  clustering_method = "complete",
  cluster_cols = T,
  cluster_rows=F,# Do not cluster populations (rows)
  scale = "none",  # No scaling
  show_rownames = TRUE,  # Show population names
  show_colnames = TRUE,  # Show miRNA names (column names)
  annotation_col = annotation_col,  # Add OCS annotation
  color = colorRampPalette(c("white", "brown3"))(10),
  main = "Top 100 Most Variable miRNAs (Clustered with OCS Annotation)",
  angle_col = 45,
  legend=T,
  annotation_legend = T
)

#Plot G
#Replace "." in AF columns with 0
data <- data %>%
  mutate(across(starts_with("AF_"), ~ replace(., . == ".", 0), .names = "{col}"))

#Select population and Seed columns
population_columns <- c("AF_joint_afr", "AF_joint_amr", "AF_joint_asj", 
                        "AF_joint_eas", "AF_joint_fin", "AF_joint_nfe", 
                        "AF_joint_mid", "AF_joint_sas", "AF_joint_remaining", "Allele_freq")

data_for_plotting <- data %>%
  select(Seed, all_of(population_columns))

#Reshape to long format
long_data <- data_for_plotting %>%
  gather(key = "Population", value = "Allele_freq", -Seed) %>%
  mutate(Allele_freq = as.numeric(Allele_freq)) %>%
  filter(Allele_freq > 0)

#Reorder Seed levels
long_data$Seed <- factor(long_data$Seed, levels = c("No_primary", "No", "Yes"))

#Plot violin + boxplot (log scale)
plot_g <- ggplot(long_data, aes(y = Allele_freq, x = Population, fill = Seed)) + 
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.9),
               outlier.shape = 16, outlier.size = 0.8, outlier.color = "black", alpha = 0.8) +
  scale_y_log10() +
  scale_fill_manual(values = c("darkred", "darkgreen", "darkblue")) +
  labs(title = "Allele Frequency Distribution by Seed Type across Populations",
       y = "Log10(Allele Frequency)", x = "Population") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

plot_a
plot_b
plot_c
plot_d
plot_e
#plot f
plot_g
