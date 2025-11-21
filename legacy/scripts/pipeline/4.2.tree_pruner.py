"""
Чтобы удалить ноду в бифуркации (или зонтике) и сохранить действительное расстояние
до оставшегося единственного (или нет) узла без промежуточной ноды, надо обратиться
по ссылке к этой ноде и вызвать метод `.delete(preserve_branch_length=True)`

"""

from random import randint
import sys

import click
from ete3 import PhyloTree

from utils import node_parent

DEFAULT_TREE_FORMAT = 1
MAX_TREE_LEN = 55_000
APPLICATION_NAME = "pipeline"


def umbreallas2singles(tree: PhyloTree):
    """delete all umbreallas; automatically resolve polytomies
    (maybe not in other cases)
    """
    leaves = tree.get_leaves()
    leaves_parents = set(map(node_parent, leaves))

    for preleaf in leaves_parents:
        cur_leaves = []
        for child in preleaf.get_children():
            if len(child) == 1:
                # this child is leaf
                cur_leaves.append(child)
        if len(cur_leaves) == 1:
            continue

        accounted_dists = set()
        for leaf in cur_leaves:
            d = leaf.dist
            if d in accounted_dists:
                leaf.delete(preserve_branch_length=True)
            else:
                accounted_dists.add(d)


def is_approx_equal(n1, n2):
    if n1 == n2:
        return True
    max_, min_ = max(n1, n2), min(n1, n2)
    cutoff = .25  # max_ >4 times more than min_ => is not equal
    if min_ / max_ > cutoff:
        return True
    return False


def is_binary_tree(tree: PhyloTree):
    nchildrens_set = set()
    for node in tree.iter_descendants():
        nchildrens_set.add(len(node.get_children()))
    if nchildrens_set == {0, 2}:
        return True
    return False


def filter_parents(nodes: set) -> list:
    """ 
    return only preterminal nodes that have only terminal childs 
    """
    pre_terminals = []
    for node in nodes:
        subtree_len = len(node)
        n_childrens = len(node.get_children())
        if subtree_len == n_childrens:
            pre_terminals.append(node)
    return pre_terminals


def random_deletion_in_bifurcations(tree: PhyloTree):
    """ delete random nodes in biburcations """
    if not is_binary_tree(tree):
        raise ValueError("tree must be binary")

    leaves = tree.get_leaves()
    preleaves = set(map(node_parent, leaves))
    preterminals = filter_parents(preleaves)

    for node in preterminals:
        childrens = node.get_children()
        assert len(childrens) == 2, "tree is not binary"

        if is_approx_equal(childrens[0].dist, childrens[1].dist):
            idx = randint(0, 1)
            childrens[idx].delete(preserve_branch_length=True)


def pruning(path_to_newick: str, max_tree_len=MAX_TREE_LEN, frmt=DEFAULT_TREE_FORMAT):
    tree = PhyloTree(path_to_newick, format=frmt)
    umbreallas2singles(tree)
    print("Umbrellas resolved")
    i = 0
    while len(tree) > max_tree_len:
        i += 1
        random_deletion_in_bifurcations(tree)  # but it's single object
        if i > 10:
            print("10 iterations of random deletions passed, increase cutoff")
            break
    return tree, i



@click.command("pruner", help="resolve 'umbrellas' and prune tree to contain less than `max-tree-len` nodes")
@click.option("--inpath", required=True, help="path to tree in newick format")
@click.option("--outpath", default=None, help="path to output newick, default: `inpath`.pruned")
@click.option("--max-tree-len", default=MAX_TREE_LEN, show_default=True, help="max allowed tree lenght after prunning")
def main(inpath: str, outpath: str, max_tree_len: int):
    outpath = outpath or f"{inpath}.pruned"

    print(f"Run command: ./tree_pruner.py {inpath} {outpath}\n")
    pruned_tree, n_iterations = pruning(inpath, max_tree_len=max_tree_len)
    print(
        f"{len(pruned_tree)} leaves in pruned tree after umbreallas deletion "
        f"and {n_iterations} iterations of random bifurcations deletion\n"
    )
    pruned_tree.write(outfile=outpath, format=DEFAULT_TREE_FORMAT)
    print("Prunned tree wrote done.\n")


if __name__ == "__main__":
    main()
