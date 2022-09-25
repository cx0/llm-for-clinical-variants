import argparse
import os
import pathlib
import string

import torch

from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import itertools
from typing import List, Tuple
import numpy as np
from utils import read_fasta


def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    # This is an efficient way to delete lowercase characters and insertion characters from a string
    deletekeys = dict.fromkeys(string.ascii_lowercase)
    deletekeys["."] = None
    deletekeys["*"] = None

    translation = str.maketrans(deletekeys)
    return sequence.translate(translation)


def read_msa(filename: str, nseq: int) -> List[Tuple[str, str]]:
    """ Reads the first nseq sequences from an MSA file, automatically removes insertions.
    The input file must be in a3m format (although we use the SeqIO fasta parser)
    for remove_insertions to work properly."""

    msa = [
        (record.description, remove_insertions(str(record.seq)))
        for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)
    ]
    return msa


def create_parser():
    parser = argparse.ArgumentParser(
        description="Perform a mutational scan with an ensemble of ESM-1v models."  # noqa
    )

    # fmt: off
    parser.add_argument(
        "--model-location",
        type=str,
        help="PyTorch model file OR name of pretrained model to download (see README for models)",
        nargs="+",
    )
    parser.add_argument(
        "--bg-mut-input",
        type=pathlib.Path,
        help="CSV file containing the list of background mutations",
    )
    parser.add_argument(
        "--seq-input",
        type=pathlib.Path,
        help="FASTA file with sequences, named by gene"
    )
    parser.add_argument(
        "--mutation-col",
        type=str,
        default="mutant",
        help="column in the background mutation file with mutations given as 'AiB'"
    )
    parser.add_argument(
        "--gene-col",
        type=str,
        default="gene",
        help="Name of gene column in mutation file"
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        help="Output path containing the predictions",
    )
    parser.add_argument(
        "--offset-idx",
        type=int,
        default=1,
        help="Offset of the mutation positions in `--mutation-col`"
    )
    parser.add_argument(
        "--max-seq-len",
        type=int,
        default=1000,
        help="max seq length, skip anything above"
    )
    parser.add_argument(
        "--mut-chunks",
        type=int,
        default=5,
        help="Number of background mutations to process at a time, tweak for avail gpu mem"
    )
    parser.add_argument(
        "--scoring-strategy",
        type=str,
        default="wt-marginals",
        choices=["wt-marginals", "pseudo-ppl", "masked-marginals"],
        help=""
    )
    parser.add_argument(
        "--msa-path",
        type=pathlib.Path,
        help="path to MSA in a3m format (required for MSA Transformer)"
    )
    parser.add_argument(
        "--msa-samples",
        type=int,
        default=400,
        help="number of sequences to select from the start of the MSA"
    )
    # fmt: on
    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
    return parser


def parse_mutation(mut, offset_idx):
    wt, idx, mt = mut[0], int(mut[1:-1]) - offset_idx, mut[-1]
    return wt, idx, mt


def mut_scan(args, model_name, bgmut, sequence, token_probs, alphabet, offset_idx):
    scores = [
            [
                model_name,
                bgmut,
                token,
                seq_idx + offset_idx,
                (token_probs[1 + seq_idx, token_idx] 
                 - token_probs[1 + seq_idx, alphabet.get_idx(residue)]
                ).cpu().numpy()
            ]
            for token, token_idx in alphabet.to_dict().items()
            for seq_idx, residue in enumerate(sequence)
    ]
    return pd.DataFrame(scores, columns=["model", args.mutation_col, "residue", "sequence_idx", "score"])



def get_all_background_sequences(sequence, bg_mutations, offset_idx):
    data = []
    for bgmut in bg_mutations:
        wt, idx, mt = parse_mutation(bgmut, offset_idx)
        data.append( (bgmut, sequence[:idx] + mt + sequence[idx+1:]) )
    return data


def compute_pppl(row, sequence, model, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]

    # modify the sequence
    sequence = sequence[:idx] + mt + sequence[(idx + 1) :]

    # encode the sequence
    data = [
        ("protein1", sequence),
    ]

    batch_converter = alphabet.get_batch_converter()

    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)

    # compute probabilities at each position
    log_probs = []
    for i in range(1, len(sequence) - 1):
        batch_tokens_masked = batch_tokens.clone()
        batch_tokens_masked[0, i] = alphabet.mask_idx
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens_masked.cuda())["logits"], dim=-1)
        log_probs.append(token_probs[0, i, alphabet.get_idx(sequence[i])].item())  # vocab size
    return sum(log_probs)


def outfname(gene):
    return args.output_path / f'mutscan_{gene}.tsv'


def main(args):
    # Load the background mutation data
    all_mut_df = pd.read_csv(args.bg_mut_input, sep='\t')
    seq_dict = read_fasta(args.seq_input)
    seq_dict = {
        gene: seq for gene, seq in seq_dict.items() 
        if len(seq) <= args.max_seq_len and not os.path.exists(outfname(gene))
    }

    # inference for each model
    for model_location in args.model_location:
        model, alphabet = pretrained.load_model_and_alphabet(model_location)
        model.eval()
        if torch.cuda.is_available() and not args.nogpu:
            model = model.cuda()
            print("Transferred model to GPU")
        # inference across all sequences
        for gene, sequence in tqdm(seq_dict.items()):
            results = pd.DataFrame()
            mut_df = all_mut_df[all_mut_df[args.gene_col] == gene]
            
            batch_converter = alphabet.get_batch_converter()

            data = get_all_background_sequences(sequence, mut_df[args.mutation_col], args.offset_idx)
            for i in tqdm(range(0, len(data), args.mut_chunks), leave=False):
                batch_labels, batch_strs, batch_tokens = batch_converter(data[i:i+args.mut_chunks])
                try:
                    with torch.no_grad():
                        if args.nogpu:
                            token_probs = torch.log_softmax(model(batch_tokens)["logits"], dim=-1)
                        else:
                            token_probs = torch.log_softmax(model(batch_tokens.cuda())["logits"], dim=-1)
                except:
                    print(f"Could not fit chunk {i} of {gene} into memory")

                if args.scoring_strategy == "wt-marginals": 
                    for j in range(token_probs.shape[0]):
                        res_df = mut_scan(
                                args,
                                model_location, 
                                data[i][0], 
                                sequence,
                                token_probs[j, :, :], 
                                alphabet, 
                                args.offset_idx
                            )
                        results =  pd.concat([
                            results,
                            mut_df.merge(
                                res_df,
                                left_on=args.mutation_col,
                                right_on=args.mutation_col,
                                how="right"
                            )
                        ]
                        )
                elif args.scoring_strategy == "masked-marginals":
                    all_token_probs = []
                    for i in tqdm(range(batch_tokens.size(1))):
                        batch_tokens_masked = batch_tokens.clone()
                        batch_tokens_masked[0, i] = alphabet.mask_idx
                        with torch.no_grad():
                            token_probs = torch.log_softmax(
                                model(batch_tokens_masked.cuda())["logits"], dim=-1
                            )
                        all_token_probs.append(token_probs[:, i])  # vocab size
                    token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
                    mut_df[model_location] = mut_df.apply(
                        lambda row: mut_scan(
                            row[args.mutation_col],
                            sequence,
                            token_probs,
                            alphabet,
                            args.offset_idx,
                        ),
                        axis=1,
                    )
                elif args.scoring_strategy == "pseudo-ppl":
                    tqdm.pandas()
                    mut_df[model_location] = mut_df.progress_apply(
                        lambda row: compute_pppl(
                            row[args.mutation_col], sequence, model, alphabet, args.offset_idx
                        ),
                        axis=1,
                    )

            results.to_csv(outfname(gene), sep='\t')


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
