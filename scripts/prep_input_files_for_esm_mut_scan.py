import sys
import pandas as pd
from utils import AminoAcidsDF, FROM_UNIPROT, read_fasta

TO_IUPAC = AminoAcidsDF.set_index("three_letter")["iupac"].to_dict()

def to_iupac_mut(mut):
    res_1, pos, res_2 = mut[:3], mut[3:-3], mut[-3:]
    if res_1 not in TO_IUPAC or res_2 not in TO_IUPAC:
        return "NA"
    return f"{TO_IUPAC[res_1]}{pos}{TO_IUPAC[res_2]}"


def main():
    var_fname= 'data/final_mutant_list.slim.tsv'
    var_df = pd.read_csv(var_fname, sep='\t')
    res_df = var_df.assign(
        MUTATION=lambda df: df["AA_CHANGE"].apply(to_iupac_mut),
        GENE=lambda df:df["GENEINFO"].apply(lambda x: x.split(":")[0])
    )
    res_df = res_df[res_df["MUTATION"] != "NA"]
    res_df.to_csv("data/final_mutant_list_mut_scan_input.tsv", sep="\t")
    seq_dict = read_fasta("data/acmg.filtered.genes.fa")
    with open("data/acmg.filtered.genes.mut_scan_input.fa", "w") as fl:
        for seqd, seq in seq_dict.items():
            uniprot_id = seqd.split()[0].split("|")[-1]
            fl.write(f">{FROM_UNIPROT[uniprot_id]}\n{seq}\n")


if __name__ == "__main__":
    main()