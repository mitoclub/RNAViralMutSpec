# usage: python3 {sys.argv[0]} path_to_input_fasta [path_to_filtered_output_fasta]

"""
drop seqs from initial fasta file before multiple alignnment
"""

import os
import sys

import tqdm

from utils import read_fasta_generator

MIN_SEQ_LEN = 29001
PATH_TO_NAMES = "data/gisaid/allowable_names.txt"

known_seq_hashes = set()
known_names = set()
BASE_CHARSET = {'\n', 'A', 'C', 'G', 'T'}


def read_allowable_names(path):
    allowable_names = set()    
    with open(path, 'r') as fin:
        for line in fin:
            name = line.rstrip('\n')
            allowable_names.add(name)
    return allowable_names


allowable_names = read_allowable_names(PATH_TO_NAMES)


def filter_NNN_and_dublicates_and_write(input_filepath, fout, drop_N=True):
    for header, seq in read_fasta_generator(input_filepath):
        # print(header)
        seq_name = header.strip(">\n")
        seq_hash = hash(seq)

        if seq_name not in allowable_names:
            continue
        if drop_N:
            for char in seq:
                if char not in BASE_CHARSET:
                    continue
        if len(seq.replace('\n', '')) < MIN_SEQ_LEN:
            continue
        if seq_name in known_names or seq_hash in known_seq_hashes:
            continue

        fout.write(header.replace(" ", "_"))
        fout.write(seq)

        known_names.add(seq_name)
        known_seq_hashes.add(seq_hash)


def parse_args():
    if len(sys.argv) == 1:
        print(
            "usage: python3 {} path_to_input_fasta|path_to_dir "
                "[path_to_filtered_output_fasta]".format(sys.argv[0])
        )
        exit(1)
        
    inp_fasta_path = sys.argv[1]
    if len(sys.argv) > 2:
        out_fasta_path = sys.argv[2]
    else:
        out_fasta_path = inp_fasta_path.rstrip('/').replace('.fasta', '')
        out_fasta_path = out_fasta_path + ".filtered.fasta"
    
    assert os.path.exists(inp_fasta_path), "No such file {}".format(inp_fasta_path)
    return inp_fasta_path, out_fasta_path


def main():
    inp_fasta_path, out_fasta_path = parse_args()
    print(inp_fasta_path, out_fasta_path, '\n')

    with open(out_fasta_path, 'w') as fout:
        if os.path.isdir(inp_fasta_path):
            basenames = os.listdir(inp_fasta_path)
            for i, file in enumerate(basenames):
                filepath = os.path.join(inp_fasta_path, file)
                print("{:04} Processing {}".format(i, filepath))
                filter_NNN_and_dublicates_and_write(filepath, fout)
                # if i == 5:
                #     break
        else:
            filter_NNN_and_dublicates_and_write(inp_fasta_path, fout)


if __name__ == "__main__":
    main()
