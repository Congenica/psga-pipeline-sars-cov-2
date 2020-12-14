from typing import Dict, List

from natsort import natsorted
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

SampleGeneAAMutation = Dict[str, Dict[str, List[str]]]


def get_mutations_per_gene_per_sample(samples: List[str], loaded_mutations: Dict, tree) -> SampleGeneAAMutation:
    """
    Build a dictionary with structure { SAMPLE : { KEY : [MUTATIONS] } }, where KEY can be (gene|amino acid)
    :param samples: the list of samples to process
    :param loaded_mutations: a dictionary with structure { KEY : [MUTATIONS] }.
    :param tree: the phylogenetic tree
    :return: the structure
    """
    sample_gene_mutations: SampleGeneAAMutation = {sample: {} for sample in samples}
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


def format_sample_gene_mutations_dict(sample_gene_mutations: SampleGeneAAMutation) -> Dict[str, str]:
    """
    Format the sample_gene_mutations. For any sample the aggregated changes should be reported as
    a semi-colon separated list eg ORF1a:E37D,S86T;ORF6:L98P;S:G222R
    :param sample_gene_mutations: the structure containing the mutations per gene per sample
    :return: the formatted structure
    """
    sample_gene_mutations_formatted = {}
    for sample, gene_mutations in sample_gene_mutations.items():
        sample_gene_mutations_formatted[sample] = ";".join(
            natsorted(
                [
                    f"{gene}:{','.join(sorted(mutations, key=sortkey_mutation_by_position))}"
                    for gene, mutations in gene_mutations.items()
                ]
            )
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
@click.option(
    "--discard-sample-ids",
    type=str,
    default="NC_045512.2",  # Wuhan sample
    help="List of sample ids to be discarded. IDs are separated by colon (e.g. s-1,s-2,s-3)",
)
def load_nextstrain_aa_muts_data(aa_muts_json: str, tree_nwk: str, discard_sample_ids: str) -> None:
    """
    Load Nextstrain amino acid mutations to the database
    """

    tree = load_phylogenetic_tree(tree_nwk)
    loaded_mutations = load_json(aa_muts_json)["nodes"]

    # list of samples without those to discard
    samples = list(set(get_samples_from_tree(tree)) - set(discard_sample_ids.split(",")))

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
