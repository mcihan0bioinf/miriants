library(data.table)
library(dplyr)
library(tidyr)

# Load target site variants
targets <- fread("data/all_target_variants.csv")
# Load miRNA binding site BED data
mirna_bs_bed <- fread("data/mirna_genomic_binding_sites.bed")

colnames(targets) <- c(
  "Chr", "Start", "ID", "Reference", "Alternative", "Quality", "Check", "Allele_count", "Allele_total", "Allele_freq",
  "AC_genomes", "AN_genomes", "AF_genomes", "AC_exomes", "AN_exomes", "AF_exomes", "grpmax_joint", "AF_joint_afr",
  "AF_joint_amr", "AF_joint_asj", "AF_joint_eas", "AF_joint_fin", "AF_joint_nfe", "AF_joint_mid", "AF_joint_sas",
  "AF_joint_remaining", "AF_exomes_afr", "AF_exomes_amr", "AF_exomes_asj", "AF_exomes_eas", "AF_exomes_fin", "AF_exomes_nfe",
  "AF_exomes_mid", "AF_exomes_sas", "AF_exomes_remaining", "AF_genomes_afr", "AF_genomes_amr", "AF_genomes_asj", "AF_genomes_eas",
  "AF_genomes_fin", "AF_genomes_nfe", "AF_genomes_mid", "AF_genomes_sas", "AF_genomes_remaining"
)
# Keep only high-quality variants
targets <- targets[Check == "PASS"]
# Rename columns for joining
setnames(mirna_bs_bed, old = c("V1", "V2", "V3"), new = c("Chr", "StartRange", "EndRange"))
targets[, position := Start]

# Join to associate variants with miRNA binding sites
targets_with_mirna <- mirna_bs_bed[targets, on = .(Chr, StartRange <= Start, EndRange >= Start), nomatch = NA]
targets_with_mirna <- targets_with_mirna[, !c("StartRange", "EndRange", "V4"), with = FALSE]

# Extract metadata from BED annotation
targets_split <- separate(
  targets_with_mirna, col = V5,
  into = c("mirna_name", "gene", "chromosome", "target_start", "target_end", "target_strand"),
  sep = "_", remove = TRUE
)

# Convert to numeric
targets_split$target_start <- as.integer(targets_split$target_start)
targets_split$target_end <- as.integer(targets_split$target_end)

# Compute relative binding site position
targets_split$bs_position <- ifelse(
  targets_split$target_strand == "+",
  targets_split$position - targets_split$target_start + 1,
  targets_split$target_end - targets_split$position + 1
)

# Load miRNA variants affecting the seed region
mirna_variants <- fread("data/mirna_variants_ocs.csv")
mirna_variants <- mirna_variants[Seed == "Yes"]
mirna_variants[, bs_position := position - 1]

# Merge miRNA and target variants by miRNA name and relative position
pairs <- merge(
  targets_split, mirna_variants,
  by = c("mirna_name", "bs_position"),
  suffixes = c("_TAR", "_MIR")
)

# Filter for reference alleles of the same purine/pyrimidine class
filtered_pairs <- pairs %>%
  filter(
    (Reference_MIR %in% c("A", "T") & Reference_TAR %in% c("A", "T")) |
      (Reference_MIR %in% c("C", "G") & Reference_TAR %in% c("C", "G"))
  ) %>%
  filter(
    nchar(Reference_TAR) == 1,
    nchar(Reference_MIR) == 1,
    nchar(Alternative_TAR) == 1,
    nchar(Alternative_MIR) == 1
  )

# Keep strand-aware matched reference allele pairs
filtered_pairs_strand <- filtered_pairs %>%
  filter(
    (Strand_mirna == target_strand & Reference_TAR == Reference_MIR) |
      (Strand_mirna != target_strand & Reference_TAR != Reference_MIR)
  )

# Filter for compensatory or identical alternative alleles
filtered_pairs_final <- filtered_pairs_strand %>%
  filter(
    (Strand_mirna == target_strand & Alternative_MIR == Alternative_TAR) |
      (Strand_mirna != target_strand & (
        (Alternative_MIR == "A" & Alternative_TAR == "T") |
          (Alternative_MIR == "T" & Alternative_TAR == "A") |
          (Alternative_MIR == "C" & Alternative_TAR == "G") |
          (Alternative_MIR == "G" & Alternative_TAR == "C")
      ))
  )
