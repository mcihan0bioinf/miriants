# Load required libraries
library(Biostrings)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(readr)
library(biomaRt)

# Load validated miRNA-gene interactions
interactions <- read_csv("validated_interactions.csv")

# Load and convert miRNA seed sequences from RNA to DNA
mirna_seeds <- readRNAStringSet("mirna_seed.fasta")
mirna_seeds_dna <- as(DNAStringSet(RNAStringSet(mirna_seeds)), "DNAStringSet")
mirna_seeds_dna <- mirna_seeds_dna[names(mirna_seeds_dna) %in% interactions$miRNA]

# Function to fetch 3' UTR genomic locations using biomaRt
get_utr_genomic_locations <- function(gene_list, batch_size = 50) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  utr_data_list <- list()
  
  for (i in seq(1, length(gene_list), by = batch_size)) {
    batch_genes <- gene_list[i:min(i + batch_size - 1, length(gene_list))]
    tryCatch({
      utr_batch <- getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "transcript_start", "transcript_end", "strand", "3utr"),
        filters = "hgnc_symbol",
        values = batch_genes,
        mart = mart
      )
      utr_data_list[[length(utr_data_list) + 1]] <- utr_batch
    }, error = function(e) message("Error: ", e$message))
  }
  
  do.call(rbind, utr_data_list)
}

# Get and process UTR sequences
target_genes <- unique(interactions$`Target Gene`)
utr_data <- get_utr_genomic_locations(target_genes) %>%
  filter(`3utr` != "Sequence unavailable") %>%
  mutate(utr_length = nchar(`3utr`)) %>%
  arrange(hgnc_symbol, desc(utr_length)) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

# Function to find genomic binding sites for a miRNA-gene pair
find_genomic_matches <- function(mirna_name, seed_seq, target_gene, utr_data) {
  gene_data <- utr_data[utr_data$hgnc_symbol == target_gene, ]
  if (nrow(gene_data) == 0) return(NULL)
  
  chrom <- gene_data$chromosome_name
  strand <- gene_data$strand
  utr_seq <- gene_data$`3utr`
  utr_len <- nchar(utr_seq)
  
  utr_genomic_start <- if (strand == 1) gene_data$transcript_end - utr_len + 1 else gene_data$transcript_start
  utr_genomic_end <- if (strand == 1) gene_data$transcript_end else gene_data$transcript_start + utr_len - 1
  
  matches <- matchPattern(seed_seq, DNAString(utr_seq))
  if (length(matches) == 0) return(NULL)
  
  match_df <- as.data.frame(ranges(matches))
  match_df$genomic_start <- if (strand == 1) utr_genomic_start + match_df$start - 1 else utr_genomic_end - match_df$end + 1
  match_df$genomic_end <- if (strand == 1) utr_genomic_start + match_df$end - 1 else utr_genomic_end - match_df$start + 1
  match_df$miRNA <- mirna_name
  match_df$Target_Gene <- target_gene
  match_df$Chromosome <- chrom
  match_df$Strand <- strand
  
  match_df
}

# Run match finding for all validated miRNA-gene pairs
all_matches <- list()
for (i in seq_len(nrow(interactions))) {
  miRNA <- interactions$miRNA[i]
  gene <- interactions$`Target Gene`[i]
  if (!(miRNA %in% names(mirna_seeds_dna))) next
  
  result <- find_genomic_matches(miRNA, mirna_seeds_dna[[miRNA]], gene, utr_data)
  if (!is.null(result)) {
    all_matches[[paste(miRNA, gene, sep = "_")]] <- result
  }
}

# Combine all results and prepare BED format
final_results <- do.call(rbind, all_matches)
final_results$chr <- paste0("chr", final_results$Chromosome)
final_results$strand <- ifelse(final_results$Strand == -1, "-", "+")
final_results$info <- paste(final_results$miRNA, final_results$Target_Gene, final_results$chr,
                            final_results$genomic_start, final_results$genomic_end,
                            final_results$strand, sep = "_")

bed_df <- final_results[, c("chr", "genomic_start", "genomic_end", "strand", "info")]

# Save as BED file
write.table(bed_df, "mirna_genomic_binding_sites.bed", col.names = FALSE,
            row.names = FALSE, quote = FALSE, sep = "\t")
