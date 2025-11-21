#!/usr/bin/python3

import sys
import pandas as pd


def get_mapper(path_to_meta):
    meta_df = pd.read_csv(path_to_meta, sep='\t')
    mapper = {sname: uid for sname, uid in
              meta_df[['Strain', 'gisaid_epi_isl']].to_dict('split')['data']}
    return mapper


def tree_formatter(path_to_tree, mapper: dict = None):
    mapper = mapper
    with open(path_to_tree) as fin:
        raw_tree = fin.read()
    new_tree = raw_tree
    for sname, uid in mapper.items():
        new_tree = new_tree.replace(sname, uid)
    return new_tree


def main():
    try:
        assert len(sys.argv) == 3, 'need 2 arguments'
        path_to_tree = sys.argv[1]
        path_to_meta = sys.argv[2]
        default_mapper = get_mapper(path_to_meta)
        tree = tree_formatter(path_to_tree, default_mapper)
        sys.stdout.write(tree)
    except Exception as e:
        sys.stderr.write(str(type(e)) + "\n")
        sys.stderr.write(str(e) + "\n")


if __name__ == '__main__':
    main()
