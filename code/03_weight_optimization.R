# Load Libraries
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)

#Load Data
#File storing microRNA variants
conservation_file <- "data/genomic_variants.csv"
#Expression data from GTEx and miRNATissueAtlas
gtex_file <- ""
atlas_file <- ""

conservation_data <- fread(conservation_file) %>%
  filter(type_mirna == "miRNA")

gtex_data <- fread(gtex_file)
atlas_data <- fread(atlas_file)

#Compute mean expression
compute_mean_expression <- function(data) {
  data %>%
    mutate(mean_expression = rowMeans(select(., -1)))
}

gtex_expression <- compute_mean_expression(gtex_data)
atlas_expression <- compute_mean_expression(atlas_data)

#Compute conservation scores
conservation_scores <- conservation_data %>%
  group_by(mirna_name) %>%
  summarize(
    SCS = ifelse(sum(Seed == "Yes") == 0, 1, 1 - mean(Allele_freq[Seed == "Yes"], na.rm = TRUE)),
    NSCS = ifelse(sum(Seed == "No") == 0, 1, 1 - mean(Allele_freq[Seed == "No"], na.rm = TRUE)),
    mirna_length = max(End_mirna) - min(Start_mirna) + 1,
    PCS = 1 - (n_distinct(Start) / mirna_length),
    TVS = 1 / (1 + n())
  ) %>%
  mutate(
    mirna_base = ifelse(
      str_count(mirna_name, "-") < 3,
      mirna_name,
      gsub("-3p|-5p|-\\d+$", "", mirna_name)
    )
  )

#Generate all possible weight combinations
generate_weights <- function() {
  x <- expand.grid(
    SCS_weight = seq(0.1, 0.7, by = 0.1),
    NSCS_weight = seq(0.1, 0.7, by = 0.1),
    PCS_weight = seq(0.1, 0.7, by = 0.1),
    TVS_weight = seq(0.1, 0.7, by = 0.1)
  ) * 10
  x$sum <- rowSums(x)
  x <- x[x$sum == 10 & x$SCS_weight > x$NSCS_weight, ]
  x <- x / 10
  x$sum <- NULL
  return(x)
}

#Compute overlap
calculate_overlap <- function(expression_data, conservation_scores, weights, threshold) {
  conservation_scores <- conservation_scores %>%
    mutate(
      OCS = SCS * weights$SCS_weight[1] +
        NSCS * weights$NSCS_weight[1] +
        PCS * weights$PCS_weight[1] +
        TVS * weights$TVS_weight[1]
    )
  
  key_column <- "mirna_name"
  
  top_conserved <- conservation_scores %>%
    arrange(desc(OCS)) %>%
    slice(1:ceiling(n() * threshold)) %>%
    pull(!!sym(key_column))
  
  top_expressed <- expression_data %>%
    arrange(desc(mean_expression)) %>%
    slice(1:ceiling(n() * threshold)) %>%
    distinct(!!sym(key_column), .keep_all = TRUE) %>%
    pull(!!sym(key_column))
  
  overlap <- length(intersect(top_conserved, top_expressed))
  precision <- overlap / length(top_conserved)
  recall <- overlap / length(top_expressed)
  f1_score <- 2 * (precision * recall) / (precision + recall + 1e-10)
  
  return(data.frame(
    SCS_weight = weights$SCS_weight[1],
    NSCS_weight = weights$NSCS_weight[1],
    PCS_weight = weights$PCS_weight[1],
    TVS_weight = weights$TVS_weight[1],
    threshold = threshold,
    overlap = overlap,
    precision = precision,
    recall = recall,
    f1_score = f1_score
  ))
}

#Evaluate the expression datasets
evaluate_dataset <- function(expression_data, conservation_scores, dataset_name) {
  weights_grid <- generate_weights()
  thresholds <- seq(0.05, 0.25, by = 0.05)
  
  map_dfr(thresholds, function(threshold) {
    map_dfr(seq_len(nrow(weights_grid)), function(i) {
      weights <- weights_grid[i, ]
      calculate_overlap(expression_data, conservation_scores, weights, threshold) %>%
        mutate(dataset = dataset_name)
    })
  })
}

evaluate <- function(expression_data, conservation_scores, dataset_name) {
  results <- evaluate_dataset(expression_data, conservation_scores, dataset_name)
  return(results)
}

results_gtex <- evaluate(gtex_expression, conservation_scores, "GTEx")
results_atlas <- evaluate(atlas_expression, conservation_scores, "Atlas")

#Compute deviations
compute_deviation_and_best <- function(results_gtex, results_atlas) {
  combined_results <- bind_rows(results_gtex, results_atlas) %>%
    group_by(SCS_weight, NSCS_weight, PCS_weight, TVS_weight, threshold) %>%
    summarize(avg_f1 = mean(f1_score), .groups = "drop") %>%
    group_by(threshold) %>%
    mutate(
      mean_f1 = mean(avg_f1),
      sd_f1 = sd(avg_f1),
      deviation = (avg_f1 - mean_f1) / (sd_f1 + 1e-10)
    ) %>%
    ungroup()
  
  deviation_results <- combined_results %>%
    group_by(SCS_weight, NSCS_weight, PCS_weight, TVS_weight) %>%
    summarize(deviation = abs(deviation), .groups = "drop")
  
  best_result <- deviation_results %>%
    top_n(1, abs(deviation))
  
  print("Best Weight Combination Based on Greatest Average Deviation from the Mean:")
  print(best_result)
  
  combined_results %>%
    mutate(weight_combination = paste(SCS_weight, NSCS_weight, PCS_weight, TVS_weight, sep = "-")) %>%
    ggplot(aes(x = weight_combination, y = deviation, color = as.factor(threshold))) +
    geom_point() +
    labs(
      title = "Deviation from Mean F1-Score for Weight Combinations",
      x = "Weight Combination (SCS-NSCS-PCS-TVS)",
      y = "Deviation (z-score)",
      color = "Threshold"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

x=compute_deviation_and_best(results_gtex, results_atlas)
