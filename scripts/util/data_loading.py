from typing import Dict
import json

from Bio import Phylo


def load_phylogenetic_tree(tree_file: str, tree_format: str = "newick"):
    """
    Load phylogenetic trees using Bio.Phylo python package.
    :param tree_file: the file containing the phylogenetic tree
    :param tree_format: a format in: [newick, nexus, nexml, phyloxml, cdao]
    :return: the loaded phylogenetic tree
    """
    with open(tree_file) as tree_handle:
        tree = Phylo.read(tree_handle, tree_format)
    return tree


def load_json(json_file: str) -> Dict:
    """
    Load a json file
    :param json_file: the json file to load
    :return: the json file as a dictionary
    """
    with open(json_file) as json_handle:
        json_dict = json.load(json_handle)
    return json_dict
