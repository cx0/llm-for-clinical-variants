import requests
import utils

if __name__ == '__main__':
    seq_dict = utils.read_fasta('data/acmg.filtered.genes.fa')
    for seqd in seq_dict.keys():
        uniprot_id = seqd.split()[0].split("|")[-1]
        af_id = seqd.split()[0].split("|")[1]
        with open(f'structures/{utils.FROM_UNIPROT[uniprot_id]}.pdb', 'wb') as fl:
            fl.write(
                requests.get(f'https://alphafold.ebi.ac.uk/files/AF-{af_id}-F1-model_v3.pdb').content
            )