import argparse
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
        "--sequence",
        type=str,
        help="Base sequence to which mutations will be applied",
    )
    parser.add_argument(
        "--bg-mut-input",
        type=pathlib.Path,
        help="CSV file containing the list of background mutations",
    )
    parser.add_argument(
        "--mutation-col",
        type=str,
        default="mutant",
        help="column in the background mutation file with mutations given as 'AiB'"
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        help="Output file containing the predictions",
    )
    parser.add_argument(
        "--offset-idx",
        type=int,
        default=0,
        help="Offset of the mutation positions in `--mutation-col`"
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


def mut_scan(model_name, bgmut, sequence, token_probs, alphabet, offset_idx):
    wt, idx, mt = parse_mutation(bgmut, offset_idx)

    mt_encoded = alphabet.get_idx(mt)
    scores = [
            [
                model_name,
                bgmut,
                token,
                seq_idx + offset_idx,
                token_probs[1 + i, ] - token_probs[1 + i, alphabet.get_idx(sequence[i])]
            ]
            for token, seq_idx in alphabet.to_dict().items()
    ]
    return pd.DataFrame(scores, columns=["model", "background_mutation", "residue", "sequence_idx", "score"])



def get_all_background_sequences(sequence, bg_mutations, offset_idx):
    data = []
    for bgmut in bg_mutations:
        wt, idx, mt = parse_mutation(bgmut, offset_idx)
        data.append( (bgmut, args.sequence[:idx] + mt + args.sequence[idx+1:]) )
    return data


def compute_pppl(row, sequence, model, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"

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


def main(args):
    # Load the background mutation data
    df = pd.read_csv(args.bg_mut_input)
    results = pd.DataFrame()

    # inference for each model
    for model_location in args.model_location:
        model, alphabet = pretrained.load_model_and_alphabet(model_location)
        model.eval()
        if torch.cuda.is_available() and not args.nogpu:
            model = model.cuda()
            print("Transferred model to GPU")

        batch_converter = alphabet.get_batch_converter()

        if isinstance(model, MSATransformer):
            data = [read_msa(args.msa_path, args.msa_samples)]
            assert (
                args.scoring_strategy == "masked-marginals"
            ), "MSA Transformer only supports masked marginal strategy"

            batch_labels, batch_strs, batch_tokens = batch_converter(data)

            all_token_probs = []
            for i in tqdm(range(batch_tokens.size(2))):
                batch_tokens_masked = batch_tokens.clone()
                batch_tokens_masked[0, 0, i] = alphabet.mask_idx  # mask out first sequence
                with torch.no_grad():
                    token_probs = torch.log_softmax(
                        model(batch_tokens_masked.cuda())["logits"], dim=-1
                    )
                all_token_probs.append(token_probs[:, 0, i])  # vocab size
            token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
            df[model_location] = df.apply(
                lambda row: mut_scan(
                    row[args.mutation_col], args.sequence, token_probs, alphabet, args.offset_idx
                ),
                axis=1,
            )

        else:
            data = get_all_background_sequences(args.sequence, df[args.mutation_col])
            batch_labels, batch_strs, batch_tokens = batch_converter(data)

            if args.scoring_strategy == "wt-marginals":
                for i in range(batch_tokens.shape[0]):
                    with torch.no_grad():
                        token_probs = torch.log_softmax(model(batch_tokens.cuda()[i, :, :])["logits"], dim=-1)
                    results =  results.append(mut_scan(model_location, data[i][0], args.sequence, token_probs, alphabet, offset_idx))
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
                df[model_location] = df.apply(
                    lambda row: mut_scan(
                        row[args.mutation_col],
                        args.sequence,
                        token_probs,
                        alphabet,
                        args.offset_idx,
                    ),
                    axis=1,
                )
            elif args.scoring_strategy == "pseudo-ppl":
                tqdm.pandas()
                df[model_location] = df.progress_apply(
                    lambda row: compute_pppl(
                        row[args.mutation_col], args.sequence, model, alphabet, args.offset_idx
                    ),
                    axis=1,
                )

    df.to_csv(args.dms_output)


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
