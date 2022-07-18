import os
from typing import List

import click
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def filter_data(data: pd.DataFrame, idx: set) -> pd.DataFrame:
    data = data.copy()
    data = data[data.date.str.contains(
        "\d{4}(-\d\d){1,2}")].sort_values("date")
    data = data[data.strain.isin(idx)]
    data["month"] = data.date.str.slice(0, 7)
    return data


def split_by_date(df: pd.DataFrame, months, size=10000) -> List[pd.DataFrame]:
    date_intervals = [
        months[:4],     # head
        months[10:11],  # intermediate
        months[-2:],    # tail
    ]
    print(df.month.value_counts().sort_index())
    splits = [
        df[df.month.isin(x)].sample(size) for x in date_intervals
    ]
    return splits


# @click.command("releaser")
# @click.argument("mulal", required=True, type=click.Path(exists=True))
# @click.argument("metadata", required=True, type=click.Path(exists=True))
# @click.argument("out", required=True, type=click.Path(exists=False))
def main(mulal, metadata, out):
    # assert not os.path.exists(out), "directory alerady exist"
    os.mkdir(out)

    reader = SeqIO.parse(mulal, "fasta")
    idx = set()
    for rec in reader:
        idx.add(rec.name)

    reader = SeqIO.parse(mulal, "fasta")  # update iterator

    cols = ["strain", "gisaid_epi_isl", "date"]
    metadata = pd.read_csv(metadata, usecols=cols)
    filtered = filter_data(metadata, idx)
    months = filtered.month.unique()

    n_splits = 3
    df_splits = split_by_date(filtered, months)
    name_splits = [set(d["strain"].values) for d in df_splits]
    record_splits = [[] for _ in range(n_splits)]
    for rec in reader:
        for i in range(n_splits):
            if rec.name in name_splits[i]:
                record_splits[i].append(rec)
                continue

    for i in range(n_splits):
        fasta_out = os.path.join(out, f"mulal_split_{i+1}.fasta")
        metadata_out = os.path.join(out, f"metadata_{i+1}.csv")

        df_splits[i].to_csv(metadata_out, index=None)
        SeqIO.write(record_splits[i], fasta_out, "fasta-2line")


if __name__ == "__main__":
    # main()

    mulal_path = "./data/mulal.fasta"
    metadata_path = "./data/gisaid/metadata_clean.csv"
    out = "./data/splits"
    main(mulal_path, metadata_path, out)
