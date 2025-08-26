# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load variants
targets <- fread("combined_target_variants.tsv")
colnames(targets) <- c("Chr","Start","ID","Reference","Alternative","Quality","Check","Allele_count","Allele_total","Allele_freq",
                       "AC_genomes","AN_genomes","AF_genomes","AC_exomes","AN_exomes","AF_exomes","grpmax_joint","AF_joint_afr",
                       "AF_joint_amr","AF_joint_asj","AF_joint_eas","AF_joint_fin","AF_joint_nfe","AF_joint_mid","AF_joint_sas",
                       "AF_joint_remaining","AF_exomes_afr","AF_exomes_amr","AF_exomes_asj","AF_exomes_eas","AF_exomes_fin","AF_exomes_nfe",
                       "AF_exomes_mid","AF_exomes_sas","AF_exomes_remaining","AF_genomes_afr","AF_genomes_amr","AF_genomes_asj","AF_genomes_eas",
                       "AF_genomes_fin","AF_genomes_nfe","AF_genomes_mid","AF_genomes_sas","AF_genomes_remaining")
targets <- targets[targets$Check == "PASS",]

# Load miRNA binding sites
mirna_bs_bed <- fread("mirna_pairs_genomic.bed")

# Harmonize chromosomes
setDT(targets); setDT(mirna_bs_bed)
mirna_bs_bed[, Chr := sub("^chr","",chromosome)]
targets[, Chr := sub("^chr","",Chr)]

# Cast types and add position
mirna_bs_bed[, `:=`(target_start=as.integer(target_start),
                    target_end=as.integer(target_end))]
targets[, position := Start]

# Interval join
targets_with_mirna <- mirna_bs_bed[
  targets, on=.(Chr, target_start <= Start, target_end >= Start), nomatch=0L
]

# Parse name field
targets_split <- separate(
  targets_with_mirna,
  col=name,
  into=c("mirna_name","gene","chromosome_name","target_start_name","target_end_name","target_strand_name"),
  sep="_", remove=FALSE
)
setDT(targets_split)
targets_split[, `:=`(
  target_start_name = as.integer(target_start_name),
  target_end_name   = as.integer(target_end_name)
)]
targets_split[, bs_position := ifelse(target_strand=="+",
                                      position - target_start + 1L,
                                      target_end - position + 1L)]

# Load conservation data
mirna_variants <- fread("genome_variants_ocs.csv")
mirna_variants <- mirna_variants[Seed=="Yes",]
mirna_variants$bs_position <- mirna_variants$position-1

# Merge
pairs <- merge(targets_split, mirna_variants,
               by=c("mirna_name","bs_position"),
               suffixes=c("_TAR","_MIR"))

# Filter
filtered_pairs <- pairs %>%
  filter(
    (Reference_MIR %in% c("A","T") & Reference_TAR %in% c("A","T")) |
      (Reference_MIR %in% c("C","G") & Reference_TAR %in% c("C","G"))
  )
filtered_pairs <- filtered_pairs[nchar(Reference_TAR)==1 &
                                   nchar(Reference_MIR)==1 &
                                   nchar(Alternative_TAR)==1 &
                                   nchar(Alternative_MIR)==1,]

# Strand-aware filtering
filtered_pairs_strand <- filtered_pairs %>%
  filter(
    (Strand_mirna==target_strand & Reference_TAR==Reference_MIR) |
      (Strand_mirna!=target_strand & Reference_TAR!=Reference_MIR)
  )
filtered_pairs_final <- filtered_pairs_strand %>%
  filter(
    (Strand_mirna==target_strand & Alternative_MIR==Alternative_TAR) |
      (Strand_mirna!=target_strand & (
        (Alternative_MIR=="A" & Alternative_TAR=="T") |
          (Alternative_MIR=="T" & Alternative_TAR=="A") |
          (Alternative_MIR=="C" & Alternative_TAR=="G") |
          (Alternative_MIR=="G" & Alternative_TAR=="C")
      ))
  )

# Synonymous/non-synonymous
variant_data <- filtered_pairs_strand %>%
  mutate(synonymous = ifelse(do.call(paste0,.) %in% do.call(paste0,filtered_pairs_final),"yes","no"))
mirna_variant_counts <- variant_data %>%
  group_by(mirna_name, OCS, synonymous) %>%
  summarise(count=n(), .groups="drop") %>%
  pivot_wider(names_from=synonymous, values_from=count, values_fill=list(count=0)) %>%
  rename(synonymous=yes, non_synonymous=no) %>%
  mutate(ratio_syn_nonSyn = ifelse(non_synonymous>0, synonymous/non_synonymous, NA))


# Allele frequency fix
filtered_pairs_final <- filtered_pairs_final %>%
  mutate(across(matches("^(AF_|Allele)"), ~ as.numeric(replace(., .==".","0"))))

# Build plotting data
plotting_data <- filtered_pairs_final[,c("Chr_MIR","Start_mirna","End_mirna",
                                         "Chr_TAR","target_start","target_end","target_strand",
                                         "Allele_freq_MIR","Allele_freq_TAR",
                                         "grpmax_joint_MIR","grpmax_joint_TAR",
                                         "OCS","mirna_name","Strand_mirna","gene",
                                         "bs_position","position_TAR","Start",
                                         "Reference_TAR","Alternative_TAR","Reference_MIR","Alternative_MIR","seed_class")]
plotting_data$Allele_freq_MIR <- as.numeric(plotting_data$Allele_freq_MIR)
plotting_data$Allele_freq_TAR <- as.numeric(plotting_data$Allele_freq_TAR)
plotting_data <- plotting_data %>%
  mutate(snv_type = case_when(
    Allele_freq_TAR > 0.01 & Allele_freq_MIR > 0.01 ~ "snp_both",
    Allele_freq_TAR < 0.01 & Allele_freq_MIR < 0.01 ~ "rare_both",
    TRUE ~ "mixed"
  ))

# Collapse seeds
plotting_data[, seed_class := if ("7mer-A1" %in% seed_class & "7mer-m8" %in% seed_class) "8mer" else seed_class,
              by=setdiff(names(plotting_data),"seed_class")]
seed_rank <- c("6mer"=1,"7mer-A1"=2,"7mer-m8"=2,"8mer"=3)
plotting_data[, rank := seed_rank[seed_class]]
plotting_data_collapsed <- plotting_data[, .SD[which.max(rank)],
                                         by=setdiff(names(plotting_data),c("seed_class","rank"))]

# Export
write.csv(plotting_data_collapsed,"compensatory_pairs_all.csv",row.names=FALSE)

