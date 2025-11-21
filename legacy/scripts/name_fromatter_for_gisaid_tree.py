#!/usr/bin/python3

from functools import partial
import re
import sys

DEFAULT_TREE_PATH = "../data/multiple_alignment_gisaid.fasta.tre-simple"
PATTERN_TO_EXTRACT_SHORT_NAME = re.compile(
    r"hCoV-19/([^\(\),]+?)\|(EPI_ISL_\d+?)\|[\d-]+?([,\)])"
)


def name_formatter(node_name, what_part_to_use=0):
    """
    :param what_part_to_use: - 0 if short name is needed, 1 if ID is needed
    """
    parts = node_name.groups()
    return "".join([parts[what_part_to_use], parts[-1]])


def names_replacer(string, what_part_to_use):
    new_string = re.sub(PATTERN_TO_EXTRACT_SHORT_NAME,
                        partial(name_formatter,
                                what_part_to_use=what_part_to_use),
                        string)
    return new_string


def tree_formatter(path_to_tree, what_part_to_use):
    with open(path_to_tree) as fin:
        raw_tree = fin.read()
    new_tree = names_replacer(raw_tree, what_part_to_use)
    return new_tree


def main():
    try:
        path_to_tree = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TREE_PATH
        what_part_to_use = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        tree = tree_formatter(path_to_tree, what_part_to_use)
        sys.stdout.write(tree)
    except Exception as e:
        print(e)


if __name__ == '__main__':
    main()


NODES_NAMES = """\
external/NC_045512.2,
hCoV-19/USA/MI-MDHHS-SC21541/2020|EPI_ISL_516268|2020-03-14,
hCoV-19/Algeria/G0638_2264/2020|EPI_ISL_418241|2020-03-02,
hCoV-19/Algeria/G0640_2265/2020|EPI_ISL_418242|2020-03-08,
hCoV-19/Algeria/G0860_2262/2020|EPI_ISL_420037|2020-03-02,
hCoV-19/Andorra/202552/2020|EPI_ISL_539496|2020-03-13,
hCoV-19/Argentina/C121/2020|EPI_ISL_420600|2020-03-07,
hCoV-19/Argentina/C1374/2020|EPI_ISL_420599|2020-03-18,
hCoV-19/Argentina/C3013/2020|EPI_ISL_420598|2020-03-22,
hCoV-19/Argentina/Heritas-HG023/2020|EPI_ISL_615121|2020-05-26,
hCoV-19/Argentina/Heritas_HG001/2020|EPI_ISL_476496|2020-04-22,\
"""


def _test_pattern_recognize_correctly():
    found_instances = re.findall(PATTERN_TO_EXTRACT_SHORT_NAME, NODES_NAMES)
    print(found_instances)


def _test_tree_formatter_works_correctly():
    extracted_units = names_replacer(NODES_NAMES, 1)
    print(extracted_units)


def _test_compare_chars_counts():
    _chars = "(),"
    with open(DEFAULT_TREE_PATH) as fin:
        base_tree = fin.read()
    sname_tree = tree_formatter(DEFAULT_TREE_PATH, 0)
    uid_tree = tree_formatter(DEFAULT_TREE_PATH, 1)

    for char in _chars:
        correct_cnt = base_tree.count(char)
        sname_cnt = sname_tree.count(char)
        uid_cnt = uid_tree.count(char)

        print(repr(char), "\t", correct_cnt, sname_cnt, uid_cnt)


def test_check_chars_in_sname():
    _chars = "(),"
    with open(DEFAULT_TREE_PATH) as fin:
        base_tree = fin.read()
    
    found_names = re.findall(PATTERN_TO_EXTRACT_SHORT_NAME, base_tree)

    for sname, uid, end_char in found_names:
        assert all([c not in sname for c in _chars])
        assert all([x in '0123456789' for x in uid.lstrip("EPI_ISL_")])
        assert end_char in _chars
