"""
read tree and multiple alignment from prank output and
write json file that contains all horisontal substitutions from tree

json format:
[
    (parent node,
    child node,
    substitution [(position_0_based, parent_nucl, child_nucl), ...]), ...
]
"""

import json
import os
from queue import Queue
import sys

import click
from ete3 import PhyloTree
import pandas as pd

from utils import node_parent, FastaStorage, release_mutations_from_two_seqs


# TODO use temporary mysql database instead of storage dict 


def get_mutations(node, storage: FastaStorage) -> tuple:
    parent = node_parent(node)
    if parent is None:
        print("WoW", file=sys.stderr)
        return
    seq_of_parent = storage.get_sequence(parent.name)
    seq_of_child = storage.get_sequence(node.name)
    assert len(seq_of_child) == len(seq_of_parent), (
        "parent and child seq lenghts aren't equal"
    )
    mutations = release_mutations_from_two_seqs(
        seq_of_parent,
        seq_of_child,
    )
    return parent.name, node.name, mutations


def bfs_to_extract_mutspec(path_to_tree: str, stor: FastaStorage) -> list:
    """ BFS for extraction """
    tree = PhyloTree(path_to_tree, format=1)
    print("tree readed", file=sys.stderr)
    discovered_nodes = set()
    Q = Queue()
    Q.put(tree)
    discovered_nodes.add(tree.name)

    overall_mutations = []
    while not Q.empty():
        cur_node = Q.get()
        for child in cur_node.children:
            Q.put(child)

        if cur_node.name not in discovered_nodes:
            discovered_nodes.add(cur_node.name)

            # main process starts here
            cur_mutations_of_one = get_mutations(cur_node, stor)
            overall_mutations.append(cur_mutations_of_one)
    print(
        f"n_discovered_nodes: {len(discovered_nodes)}, "
        f"n_overall_mutations: {len(overall_mutations)}",
        file=sys.stderr
    )
    return overall_mutations


def save_mutations_to_json(mutations: list, path: str):
    with open(path, 'w') as fout:
        json.dump(mutations, fout)


def save_mutations_to_df(mutations: list, path: str):
    mut_collection = []
    for parent, child, pair_substs in mutations:
        for pos, src_nucl, mut_nucl, src_context, mut_context in pair_substs:
            mut_collection.append((
                pos, src_nucl, mut_nucl, src_context, mut_context, parent, child
            ))
    cols = ["pos", "parent_nucl", "child_nucl", "parent_nucl_context", 
            "child_nucl_context", "parent_node", "child_node"]
    df = pd.DataFrame(mut_collection, columns=cols)
    df.to_csv(path, index=None)


@click.command(
    "extractor", 
    help="extract mutations from phylognenetic tree and dump it to json;  "
         "run script after prank completion, i.e. after ancestors reconstruction"
)
@click.option("--tree", required=True, help="path to tree in newick format")
@click.option("--fasta", required=True, help="path to fasta, containing every tree node sequence")
@click.option("--out-csv", required=True, help="path to output csv")
@click.option("--out-json", default=None, help="path to output json")
def extract_mutations(tree: str, fasta: str, out_csv: str, out_json: str):
    stor = FastaStorage(fasta)
    mutations = bfs_to_extract_mutspec(tree, stor)
    save_mutations_to_df(mutations, out_csv)
    if out_json is not None:
         save_mutations_to_json(mutations, out_json)


if __name__ == '__main__':
    extract_mutations()
