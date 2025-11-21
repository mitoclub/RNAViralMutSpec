"""
authors of one paper said that all sars-cov-2 genomes are different from Wuhan reference
only in 30-40 substitutions. We must check it
"""

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import tqdm

from utils import release_mutations_from_two_seqs

PATH_TO_DATA = '../data/omicron/'
PATH_TO_MULAL = PATH_TO_DATA + 'mulal.filtered.fasta'
PATH_TO_DISTRIBUTION = PATH_TO_DATA + 'mutation_number_distribution.json'
PATH_TO_DISTRIBUTION_DF = '../data/mutation_number_distribution.csv'


def count_distribution(max_seqs=-1):
    reader = SeqIO.parse(PATH_TO_MULAL, 'fasta')

    ref = next(reader)
    assert ref.id == "external/NC_045512.2"

    # collection_of_mutations = []
    nsubst_distribution = []
    for i, record in tqdm.tqdm(enumerate(reader), total=164688):
        cur_seq = record.seq
        mutations = release_mutations_from_two_seqs(str(ref.seq), str(cur_seq))
        # collection_of_mutations.append(mutations)
        nsubst_distribution.append(len(mutations))
        if i == max_seqs:
            break

    return nsubst_distribution


def write_to_file(max_seqs=-1):
    nsubst_distribution = count_distribution(max_seqs)

    with open(PATH_TO_DISTRIBUTION, 'w') as fout:
        json.dump(nsubst_distribution, fout)


def read_from_file():
    with open(PATH_TO_DISTRIBUTION, 'r') as fin:
        nsubst_distribution = json.load(fin)
    return nsubst_distribution


def plot_bar(nsubst_distribution):
    cnt = pd.Series(nsubst_distribution).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(list(map(str, cnt.index.values)), cnt.values, )
    plt.xticks(fontsize=6, rotation=90)
    plt.yscale('log')
    plt.show()


def build_df(nsubst_distribution):
    reader = SeqIO.parse(PATH_TO_MULAL, 'fasta')
    ref = next(reader)
    assert ref.id == "external/NC_045512.2"

    nsubst_distribution_df = []
    for i, record in tqdm.tqdm(enumerate(reader), total=164688):
        nsubst_distribution_df.append((
            record.id, nsubst_distribution[i]
        ))
    df = pd.DataFrame(nsubst_distribution_df, columns=['name', 'difference'])
    return df


# write_to_file()
nsubst_distribution = read_from_file()

# df = build_df(nsubst_distribution)
# df.to_csv(PATH_TO_DISTRIBUTION_DF, index=None)

plot_bar(nsubst_distribution)
