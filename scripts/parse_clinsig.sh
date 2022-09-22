#!/usr/bin/bash

# ClinVar variants (reference assembly: GRCh38)
# Source directory: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

# Download ClinVar VCF file (56 MB) and index
cd data/
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20220917.vcf.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20220917.vcf.gz.tbi

# Extract relevant INFO fields
echo 'CHROM\tPOS\tID\tREF\tALT\tALLELEID\tCLNREVSTAT\tCLNSIG\tCLNSIGCONF\tCLNVC\tCLNVI\tGENEINFO\tMC\n' > clinvar_20220917.filtered.tsv
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ALLELEID\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\t%INFO/CLNVC\t%INFO/CLNVI\t%INFO/GENEINFO\t%INFO/MC\n' clinvar_20220917.vcf.gz >> clinvar_20220917.filtered.tsv

# Extract records only for missense variants found in filtered subset of ACMG genes
head -1 clinvar_20220917.filtered.tsv > clinvar.acmg.tsv
awk '{print $0}' acmg.filtered.genes.list \| 
rg -f - clinvar_20220917.filtered.tsv \|
rg 'missense_variant' \|
rg 'single_nucleotide_variant' >> clinvar.acmg.tsv

# ClinVar weekly updates: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
# Download assembly-specific variant annotation (Release date: 2022-09-19) (136 MB)
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# Extract HGVS annotations
zcat < variant_summary.txt.gz | head -n 1 | cut -f 1,3,5,7 > variant_summary.hgvs.tsv
awk '{split($0,v,/:/); print v[1]}' acmg.filtered.genes.list  \| 
rg -zf - variant_summary.txt.gz \| 
rg 'GRCh38' \| 
rg 'single nucleotide variant' \| 
rg '\(p.' \| 
cut -f 1,3,5,7 >> variant_summary.hgvs.tsv