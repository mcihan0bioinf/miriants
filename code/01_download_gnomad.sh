#!/bin/bash
# Download gnomAD v4.1.0 data
# Base URL for downloading files
BASE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/"
# Array of chromosomes to download (1 to 22, X, Y)
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
# Loop through each chromosome and download the corresponding file
for chr in "${chromosomes[@]}"; do
    file="gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
    wget "${BASE_URL}${file}" 
    wget "${BASE_URL}${file}.tbi"
done
