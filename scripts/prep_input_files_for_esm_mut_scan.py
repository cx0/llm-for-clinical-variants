import sys
import pandas as pd
from utils import AminoAcidsDF, FROM_UNIPROT, read_fasta

TO_IUPAC = AminoAcidsDF.set_index("three_letter")["iupac"].to_dict()

def to_iupac_mut(mut):
    r = mut.split()[-1]
    r = r.replace("(", "").replace(")", "").replace("p.", "")
    res_1, pos, res_2 = r[:3], r[3:-3], r[-3:]
    if res_1 not in TO_IUPAC or res_2 not in TO_IUPAC:
        return "NA"
    return f"{TO_IUPAC[res_1]}{pos}{TO_IUPAC[res_2]}"


def main():
    clinvar_fname = 'data/clinvar_parsed.txt'
    clinvar_df = pd.read_csv(clinvar_fname, sep='\t')
    res_df = clinvar_df.assign(
        Mutation=lambda df: df["Name"].apply(to_iupac_mut)
    )
    res_df = res_df[res_df["Mutation"] != "NA"]
    res_df.to_csv("data/clinvar_mut_scan_input.tsv", sep="\t")
    seq_dict = read_fasta("data/acmg.filtered.genes.fa")
    with open("data/acmg.filtered.genes.mut_scan_input.fa", "w") as fl:
        for seqd, seq in seq_dict.items():
            uniprot_id = seqd.split()[0].split("|")[-1]
            fl.write(f">{FROM_UNIPROT[uniprot_id]}\n{seq}\n")


if __name__ == "__main__":
    main()