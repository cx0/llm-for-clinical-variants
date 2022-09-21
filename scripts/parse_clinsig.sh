#!/usr/bin/bash

# ClinVar weekly updates: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
# Download assembly-specific variant annotation (Release date: 2022-09-19)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

awk '{print "\t"$0"\t"}' ../data/acmg.filtered.genes.list \| 
rg -zf - variant_summary.txt.gz \| 
rg 'GRCh38' \| 
rg 'single nucleotide variant' \| 
rg '\(p.' \| 
cut -f 1,3,5,7 > ../data/clinvar_parsed.txt

# append header with the relevant columns
sed -i '.bak' '1s/^/\#AlleleID\tName\tGeneSymbol\tClinicalSignificance\'$'\n/g' ../data/clinvar_parsed.txt
