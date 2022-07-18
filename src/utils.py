import json
import sys
import time
from queue import Queue
from typing import Tuple
from collections import defaultdict

import numpy as np
import networkx as nx
import pandas as pd
from rna_tools.SecondaryStructure import parse_vienna_to_pairs
from ete3 import PhyloTree
import tqdm

PRANK_OUTPUT_PATH = "../data/mulal_gisaid_2021-01-22.filtered.twice.fasta.prank.anc"
TREE_PATH = PRANK_OUTPUT_PATH + ".dnd"
FASTA_PATH = PRANK_OUTPUT_PATH + ".fas"

REFSEQ_PATH = "../data/covid_ref.fasta"
SUBSTITUTIONS_PATH = "../data/overall_mutations_with_context.json"
SEC_STRUCT_PATH_v1 = "../data/structure_data/SARSCoV2-RNA.ss"
SEC_STRUCT_PATH_from_site = "../data/structure_data/SARSCoV2-RNA_from_site_0based.ss"


def read_fasta_generator(filepath):
    """read fasta without deleting '\n' from line ends to write that
    in the next step
    """
    with open(filepath) as fin:
        seq = header = ''
        for line in fin:
            if line.startswith(">"):
                if seq != '':
                    yield header, seq
                header = line
                seq = ''
            else:
                seq += line
        yield header, seq


# def get_sequence(node_name: str) -> str:
#     """read fasta file online and write seqs to dict if another one are needed"""
#     if node_name in fasta_storage:
#         seq = fasta_storage[node_name]
#         return seq

#     for new_node_name, new_seq in fasta_generator:
#         new_node_name = new_node_name.strip(">\n")
#         new_seq = new_seq.replace("\n", "")
#         fasta_storage[new_node_name] = new_seq
#         if new_node_name == node_name:
#             return new_seq


class FastaStorage:
    def __init__(self, path_to_fasta: str):
        self.path_to_fasta = path_to_fasta
        self.fasta_generator = read_fasta_generator(path_to_fasta)
        self.fasta_storage = dict()
    
    def get_sequence(self, node_name: str) -> str:
        """read fasta file online and write seqs to dict if another one are needed"""
        if node_name in self.fasta_storage:
            seq = self.fasta_storage[node_name]
            return seq

        for new_node_name, new_seq in self.fasta_generator:
            new_node_name = new_node_name.strip(">\n")
            new_seq = new_seq.replace("\n", "")
            self.fasta_storage[new_node_name] = new_seq
            if new_node_name == node_name:
                return new_seq


def node_parent(node):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def trim_two_seqs(seq1: str, seq2: str) -> Tuple[int]:
    """there are '------' in the start and end of seqs. Drop it!"""
    n = len(seq1)
    start_pos, stop_pos = 0, n
    for start_pos in range(n):
        if not (seq1[start_pos] == '-' or seq2[start_pos] == "-"):
            break

    for stop_pos in range(n - 1, -1, -1):
        if not (seq1[stop_pos] == '-' or seq2[stop_pos] == "-"):
            break
    stop_pos += 1
    return start_pos, stop_pos


def extract_context(seq: str, pos: int) -> str:
    """ extract context from given seq around pos (2 nucl)

    if seq = "ATCGACT" and pos = 3, return "tcGac"
    """
    prefix = seq[pos - 2: pos].lower()
    center = seq[pos].upper()
    suffix = seq[pos + 1: pos + 3].lower()
    context = prefix + center + suffix
    return context


def release_mutations_from_two_seqs(parent_seq: str, child_seq: str) -> list:
    assert isinstance(parent_seq, str) and isinstance(child_seq, str), (
        "input strings must be str-type"
    )
    start_pos, stop_pos = trim_two_seqs(parent_seq, child_seq)
    assert stop_pos - start_pos > 25000  # TODO check доп проверка того, чтоб секи были нормальные

    mutations = []
    for pos in range(start_pos, stop_pos):
        sourse_nucl = parent_seq[pos]
        mutant_nucl = child_seq[pos]
        if sourse_nucl != mutant_nucl:
            sourse_context = extract_context(parent_seq, pos)
            mutant_context = extract_context(child_seq, pos)
            mutations.append((
                pos, sourse_nucl, mutant_nucl, sourse_context, mutant_context
            ))
    return mutations


def count_two_seqs_diff(parent_seq: str, child_seq: str) -> int:
    start_pos, stop_pos = trim_two_seqs(parent_seq, child_seq)
    diff = 0
    for pos in range(start_pos, stop_pos):
        sourse_nucl = parent_seq[pos]
        mutant_nucl = child_seq[pos]
        if sourse_nucl != mutant_nucl:
            diff += 1
    return diff


def dict_from_pairs(pairs: list):
    """there is a list of pairs, that will be converted to
    dict, where each element of pair is key and value at the same time:
    {k: v, v: k} for any (k, v) in pairs
    """
    d1 = {k: v for k, v in pairs}
    d2 = {v: k for k, v in pairs}
    d1.update(d2)
    return d1


def assign_ss_types(start, stop, paired_pos_dict):
    """ return the type of nucleotides in ss;

    type variants:
    0: free
    1: stem
    2: hairpin loop
    3: bulge loop or internal loop

    @param ss_idx, int - index of one secondary cluster from output of `read_ss_file` func
    @param nucl_idx, int - index of nucleotide in the ss
    """
    n = stop - start + 1
    inside_loop = False
    ss_types = [0] * n
    cur_loop_stop = -1

    for pos in range(n):
        if pos in paired_pos_dict:
            if not inside_loop:
                cur_loop_stop = paired_pos_dict[pos]
                inside_loop = True

            cur_ss_type = 1
            last_stem_pos = pos
        else:
            if inside_loop:
                next_stem_pos = pos + 1
                # here we find next stem pos
                while next_stem_pos not in paired_pos_dict:
                    next_stem_pos += 1

                if paired_pos_dict[last_stem_pos] == next_stem_pos:
                    cur_ss_type = 2
                else:
                    cur_ss_type = 3
            else:
                cur_ss_type = 0

        ss_types[pos] = cur_ss_type

        if pos == cur_loop_stop:
            inside_loop = False

    return ss_types


def read_ss_file(filepath, from_site=False):
    """read ss file in multiple fasta format and
    return list of clusters, where there are:
    - cluster header,
    - start_pos (0-based),
    - stop_pos (0-based),
    - ss_pairs of paired nucleotides (0-based)
    - dict of pairs in both direction of pair (very usefull)
    - ss sequence"""

    if "from_site" in filepath:
        from_site = True

    clusters = []
    with open(filepath) as fin:
        for line in fin:
            if line == "\n":
                break
            if line.startswith(">"):
                header = line.strip()
            else:
                ss_seq = line.strip()
                ss_pairs = parse_vienna_to_pairs(ss_seq)
                ss_pairs_0_based = [(x - 1, y - 1) for x, y in ss_pairs[0]]
                paired_pos_dict = dict_from_pairs(ss_pairs_0_based)
                start_pos, stop_pos = map(
                    int, header.lstrip(">consensus_").lstrip(">motif_").split("-"))
                if not from_site:
                    start_pos -= 1
                    stop_pos -= 1
                ss_types = assign_ss_types(
                    start_pos, stop_pos, paired_pos_dict)
                clusters.append((
                    header,
                    start_pos,
                    stop_pos,
                    ss_pairs_0_based,
                    paired_pos_dict,
                    ss_types,
                    ss_seq,
                ))
    return clusters


def ss_graph_constructor(ss_cluster):
    n = ss_cluster[2] - ss_cluster[1] + 1
    G = nx.path_graph(n)
    pairs = ss_cluster[3]
    for v1, v2 in pairs:
        G.add_edge(v1, v2)
    return G


def full_genome_pairs_constructor(ss_clusters):
    full_genome_pairs = defaultdict(lambda: None)  # 0-based

    for ss in ss_clusters:
        start = ss[1]
        pairs = ss[4]
        for k, v in pairs.items():
            full_genome_pairs[k + start] = v + start

    full_genome_pairs = dict(full_genome_pairs)
    return full_genome_pairs


def full_ss_graph_constructor(ss_clusters, refseq_len=29903):
    full_genome_pairs = full_genome_pairs_constructor(ss_clusters)
    n = refseq_len
    G = nx.path_graph(n)
    for v1, v2 in full_genome_pairs.items():
        G.add_edge(v1, v2)
    return G


def primary_distance(pos1, pos2):
    """@param pos1, pos2: int[0..29902], positions of genome using to calculate"""
    return abs(pos1 - pos2)


def secondary_distance(pos1, pos2, full_genome_graph=None):
    """calculate distance using secondary structure

    @param pos1, pos2: int[0..29902], positions of genome using to calculate
    distance between them

    distance calculated by dijkstra_path_length in genome graph
    """
    if full_genome_graph is None:
        ss_clusters = read_ss_file(SEC_STRUCT_PATH_from_site)
        full_genome_graph = list(map(ss_graph_constructor, ss_clusters))

    dist = nx.dijkstra_path_length(full_genome_graph, pos1, pos2)
    return dist


def is_paired(nucl1, nucl2):
    """check if given nucleotides can form Watson-Crick pairing"""
    nset = {nucl1, nucl2}
    if nset == {'A', 'T'} or nset == {'C', 'G'}:
        return True
    return False


def create_full_genome_pairs(ss_clusters):
    full_genome_pairs = defaultdict(lambda: None)  # 0-based
    for ss in ss_clusters:
        start = ss[1]
        pairs = ss[4]
        for k, v in pairs.items():
            full_genome_pairs[k + start] = v + start
    full_genome_pairs = dict(full_genome_pairs)
    return full_genome_pairs


def print_ss_fasta_to_draw(ss_clusters, idx1, idx2=None):
    refseq = next(read_fasta_generator(REFSEQ_PATH))[1].replace('\n', '')
    idx2 = idx2 or idx1 + 1
    for i in range(idx1, idx2):
        cl = ss_clusters[i]
        print(cl[0])
        print(refseq[cl[1]: cl[2] + 1])
        print(cl[-1])


def determine_substitution_type(
    subst: tuple,
    parent_node: str,
    child_node: str,
    full_genome_pairs: dict
):
    """ looking at secondary structure and return type of subst:
    {complementary, not_complementary, not_in_sh}

    @param subst: tuple(pos, parent_nucl, child_nucl)
    to adress `ss_clusters`

    return stype (int): substitution type, one of:
    {
    0: комплементарные замены, приводящие к формированию комплементарных взаимодействий
    1: некомплементарные замены, разрушающие уже имеющиеся комплементарные взаимодействия
    2: замены, которые не восстанавливают комп. взаимодействия, хотя в референсе взаимодействие было
    3: вне вторичных взаимодействий
    4: еще есть комплементарные, которые заново сформированы за 1 итерацию
    }
    """
    pos, parent_nucl, child_nucl = subst[0], subst[1], subst[2]

    if pos in full_genome_pairs:
        # данный код выполняется, если
        # замутировала та, что находится в спаренном состоянии в рефсеке
        # нужно проверить ее тип, разрушила она парную связь иль нет
        the_pair_pos = full_genome_pairs[pos]

        parent_seq = get_sequence(parent_node)
        child_seq = get_sequence(child_node)

        parent_nucl_of_paired = parent_seq[the_pair_pos]
        child_nucl_of_paired = child_seq[the_pair_pos]

        # check if parent pair was existed
        if is_paired(parent_nucl, parent_nucl_of_paired):
            if not is_paired(child_nucl, child_nucl_of_paired):
                stype = 1  # parent yes, child no
            else:
                # здесь присходит парная замена с сохранением вторичной связи
                #                 print('DIN-DON!', file=sys.stderr)
                stype = 4  # самый редкий, должно быть
        else:
            # родительские нуклы не были спарены, спарены ли детские?
            if is_paired(child_nucl, child_nucl_of_paired):
                stype = 0  # parent no, child yes
            else:
                stype = 2  # parent no, child no
    else:
        stype = 3

    return stype
