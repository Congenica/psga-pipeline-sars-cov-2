import re
from typing import List


def get_samples_from_tree(tree) -> List[str]:
    """
    Return the samples from the phylogenetic tree structure
    :param tree: the phylogenetic tree (e.g. nwk tree)
    :return: a list of samples
    """
    return [terminal.name for terminal in tree.get_terminals()]


def sortkey_mutation_by_position(mutation: str) -> int:
    """
    Return the position of a mutation: C123T => 123
    """
    return int(re.search(r"\d+", mutation)[0])
