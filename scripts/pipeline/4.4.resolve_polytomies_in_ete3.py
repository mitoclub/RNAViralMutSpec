#!/usr/bin/python3
# USAGE: ./resolve_polytomies_in_ete3.py simple_tree [resolved_tree]

# http://etetoolkit.org/docs/latest/reference/reference_tree.html?highlight=resolve#ete3.TreeNode.resolve_polytomy

from ete3 import PhyloTree
import sys


def main():
    assert len(sys.argv) > 1, "need arguments for processing\nUSAGE: ./resolve_polytomies_in_ete3.py simple_tree [resolved_tree]"
    path_to_tree = sys.argv[1]
    tree = PhyloTree(path_to_tree, format=9)
    tree.resolve_polytomy()

    outfile = sys.argv[2] if len(sys.argv) > 2 else path_to_tree + ".resolved"
    tree.write(format=9, outfile=outfile)


if __name__ == "__main__":
    main()

# path = "./multiple_alignment_gisaid.fasta"
# names, seqs = [], []

# with open(path) as fin:
#     for line in fin:
#         if line.startswith('>'):
#             names.append(line.strip())
#         else:
#             seqs.append(line.strip())
