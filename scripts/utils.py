from io import StringIO
import pandas as pd
from Bio import SeqIO

AminoAcidsDF = pd.read_csv(
StringIO("""iupac	three_letter	name
A	Ala	Alanine
C	Cys	Cysteine
D	Asp	Aspartic Acid
E	Glu	Glutamic Acid
F	Phe	Phenylalanine
G	Gly	Glycine
H	His	Histidine
I	Ile	Isoleucine
K	Lys	Lysine
L	Leu	Leucine
M	Met	Methionine
N	Asn	Asparagine
P	Pro	Proline
Q	Gln	Glutamine
R	Arg	Arginine
S	Ser	Serine
T	Thr	Threonine
V	Val	Valine
W	Trp	Tryptophan
Y	Tyr	Tyrosine
X	Ter	Termination"""),
sep='\t'
)

def read_fasta(fname):
    seq_dict = {}
    for record in SeqIO.parse(fname, "fasta"):
        seq_dict[record.description] = str(record.seq)
    return seq_dict

def uniprot_ids_to_gene_names():
    return {
        k.split()[0].split("|")[-1]: [x.split("=")[-1] for x in k.split() if x.startswith("GN=")][0]
        for k in read_fasta("data/acmg.filtered.genes.fa")
    }

FROM_UNIPROT = uniprot_ids_to_gene_names()
TO_UNIPROT = {v:k for k, v in FROM_UNIPROT.items()}
