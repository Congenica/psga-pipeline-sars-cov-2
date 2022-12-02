from pathlib import Path
from typing import Dict
import click

from scripts.util.data_loading import load_yaml, write_yaml

# script for generating sars-cov-2 typing data YAML file from DATA_SOURCE (see variable below)
# Assumptions: the repo is cloned locally

PHE_LABEL = "phe-label"
WHO_LABEL = "who-label"
ALTERNATE_NAMES = "alternate-names"
BELONGS_TO_LINEAGE = "belongs-to-lineage"
DESCRIPTION = "description"
INFORMATION_SOURCES = "information-sources"
GENE = "gene"
AMINO_ACID_CHANGE = "amino-acid-change"

COVERAGE = "coverage"
COVERAGE_VALUE = 0.7  # default for all
VARIANTS = "variants"

DATA_SOURCE = "Generated from: https://github.com/phe-genomics/variant_definitions"


def processing(yaml_input_dir: Path, output_yaml_path: Path, output_txt_path: Path) -> None:
    """
    Load an extract of variant definitions and store it into one YAML file.
    """
    variant_definitions = {}
    input_paths = [f for f in yaml_input_dir.glob("*") if f.is_file() and f.suffix in [".yml", ".yaml"]]
    for input_path in input_paths:
        input_dict = load_yaml(input_path)
        phe_label = input_dict[PHE_LABEL]
        # extract fields of interest
        variants: Dict = {}
        for var in input_dict[VARIANTS]:
            if GENE in var and AMINO_ACID_CHANGE in var:
                gene = var[GENE]
                amino_acid_change = var[AMINO_ACID_CHANGE]
                if gene not in variants:
                    variants[gene] = []
                variants[gene].append(amino_acid_change)

        # use .get() so that if the key does not exist, the field is set to None
        variant_definitions[phe_label] = {
            WHO_LABEL: input_dict.get(WHO_LABEL),
            ALTERNATE_NAMES: input_dict.get(ALTERNATE_NAMES),
            BELONGS_TO_LINEAGE: input_dict.get(BELONGS_TO_LINEAGE),
            DESCRIPTION: input_dict[DESCRIPTION],
            INFORMATION_SOURCES: input_dict[INFORMATION_SOURCES],
            COVERAGE: COVERAGE_VALUE,
            VARIANTS: variants,
        }

    write_yaml(variant_definitions, output_yaml_path)
    with open(output_txt_path, "w") as fpath:
        fpath.write(DATA_SOURCE)


@click.command()
@click.option(
    "--yaml-input-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    required=True,
    help="The location to the variant definition yaml files",
)
@click.option(
    "--output-yaml-file",
    type=click.Path(dir_okay=False, writable=True),
    default="SARS-CoV-2.types.yaml",
    help="YAML file containing the variant definitions",
)
def extract_typing_data(
    yaml_input_dir: str,
    output_yaml_file: str,
) -> None:
    output_yaml_path = Path(output_yaml_file)
    output_txt_path = output_yaml_path.with_suffix(".txt")
    processing(Path(yaml_input_dir), output_yaml_path, output_txt_path)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    extract_typing_data()
