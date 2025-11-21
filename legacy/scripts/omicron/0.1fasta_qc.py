import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

PATH_TO_FASTA = "./data/omicron/sequences.fasta"
PATH_TO_OUT = "./data/omicron/sequences.filtered.fasta"

MIN_SEQ_LEN = 29001
N_SHARE_CUTOFF = 0.01


def main():
    sequences = SeqIO.parse(PATH_TO_FASTA, "fasta")
    ndrops = 0
    rec: SeqRecord = None
    with open(PATH_TO_OUT, "w") as fout:
        for rec in tqdm(sequences, total=262000):
            seq = str(rec.seq)
            n_count = seq.count("N")
            lenght = len(seq)
            if lenght < MIN_SEQ_LEN or n_count / lenght > N_SHARE_CUTOFF:
                ndrops += 1
                continue

            rec.description = rec.description.replace(" ", "_")
            rec.id = rec.name = rec.description
            SeqIO.write(rec, fout, "fasta")

    print("Dropped {} seqs".format(ndrops), file=sys.stderr)


if __name__ == "__main__":
    main()
