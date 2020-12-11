from typing import Dict, List

import click

from scripts.db.database import session_handler
from scripts.db.models import Sample
from scripts.util.load_nextstrain_data import (
    get_samples_from_tree,
    sortkey_mutation_by_position,
)
from scripts.util.data_loading import (
    load_json,
    load_phylogenetic_tree,
)

SampleNTMutation = Dict[str, List[str]]


def restructure_loaded_mutations(loaded_mutations: Dict[str, Dict[str, List[str]]]) -> SampleNTMutation:
    """
    Restructure the loaded nucletide mutations to simplify iterations. This structure contains all nodes
    :param loaded_mutations: the loaded dictionary from the Nextstrain nt_muts.json file
    :return: the restructured dictionary
    """
    return {node: muts_dict["muts"] for node, muts_dict in loaded_mutations.items()}


def get_mutations_per_sample(samples: List[str], loaded_mutations: SampleNTMutation, tree) -> SampleNTMutation:
    """
    Return a dictionary with structure { SAMPLE : [MUTATIONS] }, where the mutations per sample include
    the ancestor nodes
    :param samples: the list of samples to process
    :param loaded_mutations: a dictionary with structure { SAMPLE : [MUTATIONS] }.
    :param tree: the phylogenetic tree
    :return: the structure
    """
    sample_mutations: SampleNTMutation = {sample: [] for sample in samples}
    for sample in samples:
        for node in tree.get_path(sample):
            mutations = loaded_mutations[node.name]
            if mutations:
                sample_mutations[sample] = list(set(sample_mutations[sample]) | set(mutations))
    return sample_mutations


def load_sample_nt_mutations_to_db(session, sample_name: str, mutations: str) -> None:
    sample = session.query(Sample).filter_by(lab_id=sample_name).one_or_none()

    if not sample:
        sample = Sample(lab_id=sample_name)
        session.add(sample)

    sample.nucleotide_muts = mutations


def format_sample_mutations_dict(sample_mutations: SampleNTMutation) -> Dict[str, str]:
    """
    Format the sample_mutations. For any sample the aggregated changes should be reported as
    a colon separated list eg C3168A,T3654A,C12885T. The list is sorted by position
    :param sample_mutations: the structure containing the mutations per sample
    :return: the formatted structure
    """
    sample_mutations_formatted = {}
    for sample, mutations in sample_mutations.items():
        sample_mutations_formatted[sample] = ",".join(sorted(mutations, key=sortkey_mutation_by_position))
    return sample_mutations_formatted


def print_sample_mutations(sample_mutations: Dict):
    for sample in sample_mutations:
        print(f"\t{sample} => {sample_mutations[sample]}")


@click.command()
@click.option(
    "--nt-muts-json",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Nextstrain JSON file containing the nucleotide mutations",
)
@click.option(
    "--tree-nwk",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Nextstrain NWK file containing the phylogenetic tree of viral gene mutations (newick format)",
)
def load_nextstrain_nt_muts_data(nt_muts_json: str, tree_nwk: str) -> None:
    """
    Load Nextstrain nucleotide mutations to the database
    """

    tree = load_phylogenetic_tree(tree_nwk)
    loaded_mutations = load_json(nt_muts_json)["nodes"]
    loaded_mutations = restructure_loaded_mutations(loaded_mutations)

    samples = get_samples_from_tree(tree)
    print(f"Samples: {samples}")

    sample_mutations = get_mutations_per_sample(samples, loaded_mutations, tree)
    print("Raw sample-mutation structure:")
    print_sample_mutations(sample_mutations)

    sample_mutations_formatted = format_sample_mutations_dict(sample_mutations)
    print("Formatted sample-mutation structure:")
    print_sample_mutations(sample_mutations_formatted)

    with session_handler() as session:
        for sample_name, nt_mutations in sample_mutations_formatted.items():
            load_sample_nt_mutations_to_db(session, sample_name, nt_mutations)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_nextstrain_nt_muts_data()
