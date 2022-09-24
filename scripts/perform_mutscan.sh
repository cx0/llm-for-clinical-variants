#!/bin/bash
python scripts/esm_mut_scan.py  --model-location esm1v_t33_650M_UR90S_1 --seq-input data/acmg.filtered.genes.mut_scan_input.fa  --bg-mut-input data/final_mutant_list_mut_scan_input.tsv --mutation-col MUTATION --gene-col GENE --output-path results/ --scoring-strategy wt-marginals --mut-chunks 2
