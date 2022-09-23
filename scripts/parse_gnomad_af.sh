#!/usr/bin/bash

# Download (or process on the cloud using `hail`) gnomAD variant dataset (v2 GRCh38 liftover)
# https://gnomad.broadinstitute.org/downloads#exac-variants
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi

# Create a file with positions of interest in `bed` format
# Append REF/ALT for downstream processing
awk -F'\t' '{if (NR!=1) print "chr"$1"\t"$2"\t"$4\t"$5"}' ../data/final_mutant_list.tsv > ../data/final_mutant_list_grch38_pos.txt

# Query gnomAD variant dataset for the target sites
echo "CHROM\tPOS\tID\tREF\tALT\tFILTER\tAC\tAN\tAF\tAF_popmax\tpopmax\tallele_type\tvariant_type\tn_alt_alelles\tvep" > ../data/final_mutant_list_gnomad_allele_freq.txt
bcftools view -R ../data/final_mutant_list_grch38_pos.txt ~/tmp/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AF%INFO/AF_popmax\t%INFO/popmax\t%INFO/allele_type\t%INFO/variant_type\t%INFO/n_alt_alleles\t%INFO/vep\n' >> ../data/final_mutant_list_gnomad_allele_freq.txt
gzip ../data/final_mutant_list_gnomad_allele_freq.txt

# Match by observed allele

