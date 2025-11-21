import enum
import os
from multiprocessing import Process, Pool
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import tqdm

from utils import release_mutations_from_two_seqs


PATH_TO_MULAL = "./data/omicron/mulal2wuhan.fasta.sample"
PATH_TO_REF = "./data/raw/ideal_table.csv"
PATH_TO_OMICRON_MUTATIONS = "./data/omicron/omicron_signatures.txt"
NSEQS = 133493
THREADS = 24

reader = SeqIO.parse(PATH_TO_MULAL, 'fasta')
ref = next(reader)


def get_nums(rec):
    seq1, seq2 = str(ref.seq), str(rec.seq)
    mutations = release_mutations_from_two_seqs(seq1, seq2)
    return mutations, rec.name


def release_mutations():
    with Pool(THREADS) as p:
        muts = p.map(get_nums, reader)
    return muts


def load_specific_mutations(path: str):
    omutations = dict()
    with open(path) as fin:
        for line in fin:
            protein, mutations = line.strip().split(": ")
            mutations = mutations.split(", ")
            omutations[protein] = mutations
    return omutations


def add_nsp_to_ref(ref: pd.DataFrame, inplace=False):
    if not inplace:
        ref = ref.copy()
    nsp_positions = [
        (266, 805), (806, 2719), (2720, 8554), (8555, 10054), 
        (10055, 10972), (10973, 11842), (11843, 12091), (12092, 12685), 
        (12686, 13024), (13025, 13441),
        (900_000, 900_000), # nsp11, that we don't use
        (13442, 16236), # nsp12
        (16237, 18039), (18040, 19620), (19621, 20658), (20659, 21552), 
    ]
    ref["NspName"] = None
    for i, (start, end) in enumerate(nsp_positions, 1):
        ref["NspName"].at[(ref.Pos >= start) & (ref.Pos <= end)] = f"nsp{i}"
    
    return ref



def load_ref_annot(path_to_ref):
    cols = [
        'Pos', 'GenName', 'GenType', 'RefNuc', 'AltNuc', 
        'RefCodon', 'AltCodon', 'NucInCodon', 'RefAa', 'AltAa',
    ]
    ref = pd.read_csv(path_to_ref, usecols=cols)[cols]

    ref["AaInProtein"] = None
    for gene, c in ref[ref.GenType == "translated"].GenName.value_counts().iteritems():
        if gene == "ORF1ab":
            continue
        genelen = c / 9
        assert genelen == int(genelen)
        genelen = int(genelen)
        aa_numbers = np.repeat(np.arange(1, genelen + 1), 9)
        assert aa_numbers[-1] == genelen and aa_numbers.shape[0] == c
        ref["AaInProtein"].at[(ref.GenType == "translated") & (ref.GenName == gene)] = aa_numbers

    add_nsp_to_ref(ref, inplace=True)



    # TODO drop columns and add nsp labels and aa positions
    # TODO from table to set of distinct signatures

    return ref


def process_one_mut(mut: List[tuple]) -> float:
    """ check one record of mutations and return prob of "omicronity" """

    for pos, refnuc, altnuc in mut:
        Pos = pos + 1

    return 0.0 
    

def check_mutations(muts):
    ref_annot = load_ref_annot(PATH_TO_REF)

    for mut, name in muts:
        pass  # TODO


if __name__ == "__main__":
    load_specific_mutations(PATH_TO_OMICRON_MUTATIONS)
    load_ref_annot(PATH_TO_REF)
