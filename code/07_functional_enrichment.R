# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Load conservation data with precomputed OCS
conservation_file_path <- "data/mirna_variants_ocs.csv"
conservation_data <- fread(conservation_file_path) %>%
  filter(type_mirna == "miRNA") %>%
  mutate(OCS = as.numeric(OCS))

# Select top and bottom 124 miRNAs by OCS
top_mirnas <- conservation_data %>%
  arrange(desc(OCS)) %>%
  slice_head(n = 124) %>%
  pull(mirna_name)

bottom_mirnas <- conservation_data %>%
  arrange(OCS) %>%
  slice_head(n = 124) %>%
  pull(mirna_name)

# Run enrichment analysis
top_enrichment <- rba_mieaa_enrich(
  test_set = top_mirnas,
  sig_level = 1,
  mirna_type = "mature",
  test_type = "ORA",
  species = 9606
)

bottom_enrichment <- rba_mieaa_enrich(
  test_set = bottom_mirnas,
  sig_level = 1,
  mirna_type = "mature",
  test_type = "ORA",
  species = 9606
)

# Label groups
top_enrichment$Group <- "Top 5%"
bottom_enrichment$Group <- "Bottom 5%"

# Combine enrichment data
combined_enrichment <- bind_rows(top_enrichment, bottom_enrichment)

# Count over-/under-represented categories
overview_counts <- combined_enrichment %>%
  group_by(Group, Category, Enrichment) %>%
  summarise(Count = n(), .groups = "drop")

# Final summary plot
ggplot(overview_counts, aes(x = Category, y = Count, fill = Enrichment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Group, scales = "free_x") +
  coord_flip() +
  labs(
    title = "Counts of Over- and Under-represented Categories\nTop vs Bottom 5% Conserved miRNAs",
    x = "Category",
    y = "Count"
  ) +
  scale_fill_manual(values = c("over-represented" = "blue", "under-represented" = "red")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
