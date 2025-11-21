import random

import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from utils import count_two_seqs_diff

PATH_TO_MULAL_IN = "./data/mulal.fasta"
MAX_MUT_NUM = 80
DROP_PROB = 0.82
MAX_DEL_FRACTION = 1 / 30

NSEQS = 134000  # number of records in input fasta


def get_mut_num(rec: SeqRecord, ref: SeqRecord) -> int:
    """ return number of mutations relatively to refseq """
    mut_num = count_two_seqs_diff(str(ref.seq), str(rec.seq))
    return mut_num


def mulal_filtrator(inpath: str, max_mut: int, drop_prob: float, del_frac: float) -> SeqRecord:
    """
    generator of filtered records

    filtration:
        - drop seqs highly different from refseq
        - drop $DROP_PROB seqs by random 

    params:
        inpath - path to input multiple alingnment fasta
    """
    reader = SeqIO.parse(inpath, "fasta")
    ref = next(reader)
    ndropped = 0
    for rec in tqdm(reader, total=NSEQS):
        if random.random() < drop_prob:
            ndropped += 1
            continue

        if str(rec.seq).count("-") / len(rec.seq) > del_frac:
            ndropped += 1
            continue
        mut_num = get_mut_num(rec, ref)
        if mut_num > max_mut:
            ndropped += 1
            continue

        yield rec
    print("Dropped: {} seqs".format(ndropped))


def fasta_writer(seqs, handle):
    """ write records to file

    params:
        seqs - iterable of seqs
        handle - filepath or file handle
    """
    SeqIO.write(seqs, handle, "fasta-2line")


@click.command("qc")
@click.option("--inpath", help="path to input mulal fasta")
@click.option("--outpath", default=None, help="path to output mulal fasta, default: base(inpath).filtered.fasta")
@click.option("--max-mut", default=MAX_MUT_NUM, help="maximum allovable num of mutations")
@click.option("--drop-prob", default=DROP_PROB, help="probability to drop random sequence")
@click.option("--del-frac", default=MAX_DEL_FRACTION, help="maximum fraction of deletions (-) in sequence")
def main(inpath: str, outpath: str, max_mut: int, drop_prob: float, del_frac: float):
    outpath = outpath or inpath.replace("fasta", "filtered.fasta")
    filtered_seqs = mulal_filtrator(inpath, max_mut, drop_prob, del_frac)
    fasta_writer(filtered_seqs, outpath)


if __name__ == "__main__":
    main()
