#Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)
library(tidyr)

#Define file paths
#TCGA Expression of microRNAs. 
expression_file <- ""
conservation_file <- "data/mirna_variants_ocs.csv"
#MicroT database for microRNA-Target interactions within 3'UTRs and mre score of >0.7.
filtered_targets_file <- ""

#Plot A
#Load and preprocess conservation data
conservation_data <- read.csv(conservation_file) %>%
  filter(type_mirna == "miRNA") %>%
  mutate(
    mirna_base = gsub("-3p|-5p", "", mirna_name),
    mirna_base = gsub("mir", "miR", mirna_base),
    OCS = as.numeric(OCS)
  )

#Load and preprocess TCGA expression data
tcga_expr <- read.csv(expression_file, row.names = 1) %>%
  rownames_to_column("mirna_name") %>%
  mutate(
    mean_expr = rowMeans(across(where(is.numeric)), na.rm = TRUE),
    mirna_base = gsub("-3p|-5p", "", mirna_name),
    mirna_base = gsub("mir", "miR", mirna_base)
  )

#Merge datasets by normalized miRNA names
tcga_combined <- inner_join(tcga_expr, conservation_data, by = "mirna_base") %>%
  select(mean_expr, mirna_name.x, OCS, mirna_base) %>%
  distinct()

#Compute top/bottom 5% OCS thresholds
top_5_threshold <- quantile(tcga_combined$OCS, 0.95, na.rm = TRUE)
bottom_5_threshold <- quantile(tcga_combined$OCS, 0.05, na.rm = TRUE)

#Classify OCS groups
tcga_combined <- tcga_combined %>%
  mutate(color = case_when(
    OCS >= top_5_threshold ~ "Top 5%",
    OCS <= bottom_5_threshold ~ "Bottom 5%",
    TRUE ~ "Other"
  ))

#Plot expression vs OCS
plot_a <- ggplot(tcga_combined, aes(x = mean_expr, y = OCS, color = color)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Top 5%" = "blue", "Bottom 5%" = "red", "Other" = "gray")) +
  labs(title = "TCGA Mean Expression vs. OCS", x = "Mean Expression", y = "OCS") +
  scale_x_continuous(labels = scales::label_scientific()) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", legend.title = element_blank())


#Plot B

# Prepare Conservation Comparison
conservation_comparison <- conservation_data %>%
  filter(grepl("3p|5p", mirna_name)) %>%
  mutate(base_name = gsub("-3p|-5p", "", mirna_name)) %>%
  group_by(base_name, variant_type = ifelse(grepl("3p", mirna_name), "3p", "5p")) %>%
  summarise(mean_OCS = mean(OCS, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = variant_type, values_from = mean_OCS) %>%
  filter(!is.na(`3p`) & !is.na(`5p`)) %>%
  mutate(difference = `3p` - `5p`)

conservation_summary <- conservation_comparison %>%
  group_by(base_name) %>%
  summarise(
    min_difference = min(difference, na.rm = TRUE),
    max_difference = max(difference, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(min_difference)

conservation_summary$y_index <- seq_len(nrow(conservation_summary))  

plot_b <- ggplot(conservation_summary, aes(x = min_difference, y = y_index)) +
  geom_line(aes(x = max_difference, y = y_index, 
                color = ifelse(max_difference > 0, "Positive", "Negative")), size = 1.2) +  
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +  
  labs(
    title = "3p vs 5p Conservation Differences",
    x = "Difference in Conservation Score",
    y = "Ordered microRNAs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_blank(),  # Hide y-axis labels for clarity
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    legend.position = "none"  # Remove legend if not needed
  ) +
  geom_vline(xintercept = 0, color = "grey50", size = 0.7)  # Reference line at zero

#Plot C

# Load and summarize microRNA target data
targets <- read.csv(filtered_targets_file, sep = "\t")

#Load conservation and target data
conservation_data <- conservation_data %>%
  filter(type_mirna == "miRNA", grepl("3p|5p", mirna_name)) %>%
  mutate(
    base_name = gsub("-3p|-5p", "", mirna_name),
    variant_type = ifelse(grepl("3p", mirna_name), "3p", "5p"),
    OCS = as.numeric(OCS)
  ) %>%
  group_by(base_name, variant_type) %>%
  summarise(mean_OCS = mean(OCS, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = variant_type, values_from = mean_OCS) %>%
  filter(!is.na(`3p`) & !is.na(`5p`)) %>%
  mutate(difference = `3p` - `5p`)

targets_summary <- read.csv(filtered_targets_file, sep = "\t") %>%
  filter(grepl("3p|5p", mirna)) %>%
  mutate(
    base_name = gsub("-3p|-5p", "", mirna),
    variant_type = ifelse(grepl("3p", mirna), "3p", "5p")
  ) %>%
  group_by(base_name, variant_type) %>%
  summarise(num_targets = n(), .groups = "drop") %>%
  pivot_wider(names_from = variant_type, values_from = num_targets, names_prefix = "targets_") %>%
  filter(!is.na(targets_3p) & !is.na(targets_5p))

#Merge and prepare data
conservation_comparison <- conservation_data %>%
  left_join(targets_summary, by = "base_name") %>%
  mutate(
    more_conserved = ifelse(`3p` > `5p`, "3p", "5p"),
    more_targets = ifelse(targets_3p > targets_5p, "3p", "5p"),
    consistent = (more_conserved == more_targets)
  ) %>%
  na.omit()

#Group into quantiles and count consistent/inconsistent cases
conservation_comparison <- conservation_comparison %>%
  mutate(abs_difference = abs(difference), group = ntile(abs_difference, 4))

group_counts <- conservation_comparison %>%
  group_by(group, consistent) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(consistent = factor(consistent, levels = c(TRUE, FALSE)))

#Create grouped bar plot
plot_c <- ggplot(group_counts, aes(x = factor(group), y = count, fill = consistent)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  labs(
    title = "True vs False Ratio Across 4 Groups",
    x = "Difference Groups (1 = Lowest, 4 = Highest)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "top")