import re

gene_list = '''ACTA2	ENST00000224784
ACTC1	ENST00000290378
ACVRL1	ENST00000388922
BAG3	ENST00000369085
BMPR1A	ENST00000372037
BTD	ENST00000643237
CASQ2	ENST00000261448
DES	ENST00000373960
DSC2	ENST00000280904
ENG	ENST00000373203
GAA	ENST00000302262
GLA	ENST00000218516
HFE	ENST00000357618
HNF1A	ENST00000257555
KCNQ1	ENST00000155840
LDLR	ENST00000558518
LMNA	ENST00000677389
MAX	ENST00000264444
MEN1	ENST00000450708
MLH1	ENST00000231790
MSH2	ENST00000233146
MUTYH	ENST00000456914
MYL2	ENST00000228841
MYL3	ENST00000292327
NF2	ENST00000338641
OTC	ENST00000039007
PCSK9	ENST00000302118
PKP2	ENST00000340811
PMS2	ENST00000265849
PRKAG2	ENST00000287878
PTEN	ENST00000371953
RB1	ENST00000267163
RPE65	ENST00000262340
SDHAF2	ENST00000301761
SDHB	ENST00000375499
SDHC	ENST00000367975
SDHD	ENST00000375549
SMAD3	ENST00000327367
SMAD4	ENST00000342988
STK11	ENST00000326873
TGFBR1	ENST00000374994
TGFBR2	ENST00000295754
TMEM127	ENST00000258439
TMEM43	ENST00000306077
TNNC1	ENST00000232975
TNNI3	ENST00000344887
TNNT2	ENST00000656932
TP53	ENST00000546954
TPM1	ENST00000403994
TRDN	ENST00000334268
TTR	ENST00000237014
VHL	ENST00000286428
WT1	ENST00000621533'''

enst2gene = dict([a.split('\t')[::-1] for a in gene_list.split("\n")])

# Download FASTA alignments of 99 vertebrate genomes with human for CDS regions
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/alignments/knownGene.multiz100way.exonAA.fa.gz
# FASTA alignment is provided at the exon level - we need to concatenate. 

fasta_file = open('/Users/onuralp/Desktop/knownGene.multiz100way.exonAA.fa', 'r')
lines = fasta_file.readlines()

protein2sequence = {}

for line in lines:
    header = re.search('^>(\S+)', line)

    if header:
        enst_id, species = re.search('^>(ENST\d+)\.\d+(\_\S+)\_\d+\_\d+', line).groups()
        protein_id = enst_id + species

        if protein2sequence.get(protein_id) is None:
            protein2sequence[protein_id] = ""
    else:
        protein2sequence[protein_id] += line.strip()

# write MSA for each protein to file
for k,v in protein2sequence.items():
    enst_id = k.split('_')[0]
    
    if enst_id in list(enst2gene.keys()):
        with open('/Users/onuralp/Desktop/' + enst2gene[enst_id]+'.txt', 'a') as f:
            f.write(">" + k + "\n" + v + "\n")

