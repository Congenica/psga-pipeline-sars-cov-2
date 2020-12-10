from typing import Dict, List
import json

from natsort import natsorted
import click
from Bio import Phylo

from scripts.db.database import session_handler
from scripts.db.models import Sample

SampleGeneMutation = Dict[str, Dict[str, List[str]]]


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


def get_samples_from_tree(tree) -> List[str]:
    """
    Return the samples from the phylogenetic tree structure
    :param tree: the phylogenetic tree (e.g. nwk tree)
    :return: a list of samples
    """
    return [terminal.name for terminal in tree.get_terminals()]


def get_mutations_per_gene_per_sample(samples: List[str], loaded_mutations: Dict, tree) -> SampleGeneMutation:
    """
    Build a dictionary with structure { SAMPLE : { KEY : [MUTATIONS] } }, where KEY can be (gene|amino acid)
    :param samples: the list of samples to process
    :param loaded_mutations: a dictionary with structure { KEY : [MUTATIONS] }.
    :param tree: the phylogenetic tree
    :return: the structure
    """
    sample_gene_mutations: SampleGeneMutation = {sample: {} for sample in samples}
    for sample in samples:
        for node in tree.get_path(sample):
            for gene, mutations in loaded_mutations[node.name]["aa_muts"].items():
                if mutations:
                    if gene in sample_gene_mutations[sample]:
                        sample_gene_mutations[sample][gene] = list(
                            set(sample_gene_mutations[sample][gene]) | set(mutations)
                        )
                    else:
                        sample_gene_mutations[sample][gene] = mutations
    return sample_gene_mutations


def load_sample_aa_mutations_to_db(session, sample_name: str, mutations: str) -> None:
    sample = session.query(Sample).filter_by(lab_id=sample_name).one_or_none()

    if not sample:
        sample = Sample(lab_id=sample_name)
        session.add(sample)

    sample.amino_acid_muts = mutations


def format_sample_gene_mutations_dict(sample_gene_mutations: SampleGeneMutation) -> Dict[str, str]:
    """
    Format the sample_gene_mutations. For any sample the aggregated changes should be reported as
    a semi-colon separated list eg ORF1a:E37D,S86T;ORF6:L98P;S:G222R
    :param sample_gene_mutations: the structure containing the mutations per gene per sample
    :return: the formatted structure
    """
    sample_gene_mutations_formatted = {}
    for sample, gene_mutations in sample_gene_mutations.items():
        sample_gene_mutations_formatted[sample] = ";".join(
            natsorted([f"{gene}:{','.join(natsorted(mutations))}" for gene, mutations in gene_mutations.items()])
        )
    return sample_gene_mutations_formatted


def print_sample_gene_mutations(sample_gene_mutations: Dict):
    for sample in sample_gene_mutations:
        print(f"\t{sample} => {sample_gene_mutations[sample]}")


@click.command()
@click.option(
    "--aa-muts-json",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Nextstrain JSON file containing the amino acid mutations",
)
@click.option(
    "--tree-nwk",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Nextstrain NWK file containing the phylogenetic tree of viral gene mutations (newick format)",
)
def load_nextstrain_aa_muts_data(aa_muts_json: str, tree_nwk: str) -> None:
    """
    Load Nextstrain amino acid mutations to the database
    """

    tree = load_phylogenetic_tree(tree_nwk)
    loaded_mutations = load_json(aa_muts_json)["nodes"]

    samples = get_samples_from_tree(tree)
    print(f"Samples: {samples}")

    sample_gene_mutations = get_mutations_per_gene_per_sample(samples, loaded_mutations, tree)
    print("Raw sample-gene-mutation structure:")
    print_sample_gene_mutations(sample_gene_mutations)

    sample_gene_mutations_formatted = format_sample_gene_mutations_dict(sample_gene_mutations)
    print("Formatted sample-gene-mutation structure:")
    print_sample_gene_mutations(sample_gene_mutations_formatted)

    with session_handler() as session:
        for sample_name, aa_mutations in sample_gene_mutations_formatted.items():
            load_sample_aa_mutations_to_db(session, sample_name, aa_mutations)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_nextstrain_aa_muts_data()
