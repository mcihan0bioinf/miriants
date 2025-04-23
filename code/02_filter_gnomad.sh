#!/bin/bash
# This script extracts genomic variant attributes for the microRNA coordinates from miRBase from the gnomAD for downstream analysis.
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}
# Paths
input_file="gnomad.joint.v4.1.sites.${chr}.vcf.bgz"
#input bed file with genomic coordinates for microRNA targets, mature microRNA sequence or precursor microRNAs
bed_file=""
output_file="joint_variants_${chr}.tsv"

# Extract and filter VCF
bcftools view -R $bed_file $input_file | \
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/BOTH_FILTERED\t%INFO/EXOMES_FILTERED\t%INFO/GENOMES_FILTERED\t%INFO/AC_joint\t%INFO/AN_joint\t%INFO/AF_joint\t%INFO/AC_genomes\t%INFO/AN_genomes\t%INFO/AF_genomes\t%INFO/AC_exomes\t%INFO/AN_exomes\t%INFO/AF_exomes\t%INFO/AF_joint_sas\t%INFO/AC_joint_sas\t%INFO/AN_joint_sas\t%INFO/AF_genomes_sas\t%INFO/AC_genomes_sas\t%INFO/AN_genomes_sas\t%INFO/AF_exomes_sas\t%INFO/AC_exomes_sas\t%INFO/AN_exomes_sas\t%INFO/AF_joint_eas\t%INFO/AC_joint_eas\t%INFO/AN_joint_eas\t%INFO/AF_genomes_eas\t%INFO/AC_genomes_eas\t%INFO/AN_genomes_eas\t%INFO/AF_exomes_eas\t%INFO/AC_exomes_eas\t%INFO/AN_exomes_eas\t%INFO/AF_joint_nfe\t%INFO/AC_joint_nfe\t%INFO/AN_joint_nfe\t%INFO/AF_genomes_nfe\t%INFO/AC_genomes_nfe\t%INFO/AN_genomes_nfe\t%INFO/AF_exomes_nfe\t%INFO/AC_exomes_nfe\t%INFO/AN_exomes_nfe\t%INFO/AF_joint_fin\t%INFO/AC_joint_fin\t%INFO/AN_joint_fin\t%INFO/AF_genomes_fin\t%INFO/AC_genomes_fin\t%INFO/AN_genomes_fin\t%INFO/AF_exomes_fin\t%INFO/AC_exomes_fin\t%INFO/AN_exomes_fin\t%INFO/AF_joint_afr\t%INFO/AC_joint_afr\t%INFO/AN_joint_afr\t%INFO/AF_genomes_afr\t%INFO/AC_genomes_afr\t%INFO/AN_genomes_afr\t%INFO/AF_exomes_afr\t%INFO/AC_exomes_afr\t%INFO/AN_exomes_afr\t%INFO/AF_joint_amr\t%INFO/AC_joint_amr\t%INFO/AN_joint_amr\t%INFO/AF_genomes_amr\t%INFO/AC_genomes_amr\t%INFO/AN_genomes_amr\t%INFO/AF_exomes_amr\t%INFO/AC_exomes_amr\t%INFO/AN_exomes_amr\t%INFO/AF_joint_asj\t%INFO/AC_joint_asj\t%INFO/AN_joint_asj\t%INFO/AF_genomes_asj\t%INFO/AC_genomes_asj\t%INFO/AN_genomes_asj\t%INFO/AF_exomes_asj\t%INFO/AC_exomes_asj\t%INFO/AN_exomes_asj\t%INFO/AF_joint_mid\t%INFO/AC_joint_mid\t%INFO/AN_joint_mid\t%INFO/AF_genomes_mid\t%INFO/AC_genomes_mid\t%INFO/AN_genomes_mid\t%INFO/AF_exomes_mid\t%INFO/AC_exomes_mid\t%INFO/AN_exomes_mid\t%INFO/grpmax_joint\n' > $output_file

# Print message to track progress
echo "Filtered and processed ${chr} -> $output_file"

#After all jobs finish, combine the files and add a header
if [ $SLURM_ARRAY_TASK_ID -eq 23 ]; then
  combined_file="combined_joint_variants.tsv"
  #Extract the header from one of the output files
  first_file=$(ls joint_variants_*.tsv | head -n 1)
  head -n 1 "$first_file" > "$combined_file"
  # Concatenate all data, skipping headers
  tail -n +2 -q joint_variants_*.tsv >> "$combined_file"
  echo "Combined all chromosome files into $combined_file"
fi

