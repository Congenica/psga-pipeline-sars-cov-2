import shutil
from pathlib import Path
from typing import Any, Callable, Dict
import re
import csv
import pickle as Pickle
import tempfile
import click
from git import Repo
from Bio import SeqIO
import ahocorasick

from scripts.primer_cols import PRIMER_INDEX_COLS, PRIMER_NAME, FASTA_PATH, PICKLE_PATH, TOTAL_NUM_PRIMER

SARS_COV_2 = "SARS-CoV-2"
EPI2ME_LABS = "epi2me-labs"
QUICK_LAB = "quick-lab"
SUPPORTED_PATHOGENS = [SARS_COV_2]

SCAFFOLD, START, END, NAME, SCORE, STRAND = "scaffold", "start", "end", "name", "score", "strand"
BED_COLS = [
    SCAFFOLD,
    START,
    END,
    NAME,
    SCORE,
    STRAND,
]

REFERENCE = "reference"
SCHEME = "scheme"
BED = "bed"
FASTA = "fasta"
PICKLE = "pickle"

SCHEME_TYPE = Dict[str, Dict[str, Path]]


def create_automaton(input_path: Path, filetype: str = FASTA) -> ahocorasick.Automaton:
    """
    Store the reads in an Aho-Corasick automaton trie structure for primer look up

    :param input_path: the path of the input file to create the automaton for. File to be accessible with SeqIO.
    :param filetype: the type of file (e.g. FASTA)
    :return: the Aho-Corasick automaton object
    """
    automaton = ahocorasick.Automaton()
    with open(input_path, "r") as fd:
        for record in SeqIO.parse(fd, filetype):
            seq = str(record.seq)
            automaton.add_word(seq, seq)
    return automaton


def dump_pickle(output_path: Path, object_to_pickle: Any) -> None:
    """
    Dump a Python object to file using pickle

    :param output_path: the path of the pickle output file
    :param object_to_pickle: the Python object to dump using pickle
    """
    with open(output_path, "wb") as pickle_out:
        Pickle.dump(object_to_pickle, pickle_out)


def _decorate_with_new_line(method: Callable) -> Callable:
    """
    A decorator to automatically add a new line
    :param method: a method to write something
    :return: the decorated method
    """

    def decorated(text):
        method(f"{text}\n")

    return decorated


def extract_primer_sequences(
    ref_fasta: Path,
    scheme_bed: Path,
    scheme_fasta: Path,
    primer: str,
    strands_in_name: list[str] = None,
    left_strand: str = None,
) -> None:
    """
    Extract primer sequences from ref_fasta using coordinates in scheme_bed
    and generate a scheme_fasta containing these sequences.
    Both left and right primer are needed in PCR sequecing, but here we only need the left primer,
    which aligns with the forward read, in order to detect the primers in samples.
    * left primer: forward read
    * right primer: reverse/complement read

    Some scheme bed files do not have `strand` (+, -).
    In these cases, the strand is obtained from the sequence name (e.g. _LEFT, _RIGHT).
    An error is raised if the strand cannot be inferred

    :param ref_fasta: the reference FASTA file path
    :param scheme_bed: the BED file containing coordinates of scheme primers
    :param scheme_fasta: the FASTA file containing sequences of scheme primers
    :param primer: the name of the primer (e.g. ARTIC_V4)
    :param strands_in_name: a list containing the strand names. These can be present in the primer name in the BED file
    :param left_strand: the name of the left strand, if specified in the primer name in the BED file
    """

    if not strands_in_name:
        strands_in_name = ["LEFT", "RIGHT"]
        left_strand = "LEFT"

    if not (ref_fasta and ref_fasta.is_file()):
        raise FileNotFoundError(f"File {ref_fasta} was not found")

    if not (scheme_bed and scheme_bed.is_file()):
        raise FileNotFoundError(f"File {scheme_bed} was not found")

    # load the reference sequence
    # the reference fasta contains 1 single record
    for seq_record in SeqIO.parse(ref_fasta, FASTA):
        ref_sequence = str(seq_record.seq)

    with open(scheme_bed, newline="") as bedfile, open(scheme_fasta, "w") as outputfile:
        write = _decorate_with_new_line(outputfile.write)

        def _write_primer(primer: str, bed_row: dict[str, str], seq_idx: int) -> None:
            start = int(bed_row[START])
            end = int(bed_row[END])
            primer_sequence = ref_sequence[start:end]
            write(f">{primer}_primer_seq_{seq_idx}")
            write(primer_sequence)

        bed_reader = csv.DictReader(bedfile, fieldnames=BED_COLS, delimiter="\t")
        for seq_idx, bed_row in enumerate(bed_reader):
            # we only extract the forward primer, which matches the forward read
            # and, if present, is stored at the beginning of the read
            if bed_row[STRAND]:
                if bed_row[STRAND] == "+":
                    # forward strand, extract the sequence
                    _write_primer(primer, bed_row, seq_idx)
            else:
                # bed_row[STRAND] is None, infer the strand from name if possible or raise error
                try:
                    inferred_strand = next(s for s in strands_in_name if s in bed_row[NAME])
                    if inferred_strand == left_strand:
                        _write_primer(primer, bed_row, seq_idx)
                except StopIteration as unknown_strand:
                    raise ValueError(
                        f"Unknown strand in {scheme_bed}. Cannot extract primer sequences"
                    ) from unknown_strand


def prepare_dest_files(dest_schemes_path: Path, scheme_name: str, pathogen: str, version_path: Path) -> SCHEME_TYPE:
    """
    Prepare the primers destination files (e.g. reference fasta, scheme BED/FASTA/PICKLE)

    :param dest_schemes_path: the path storing the destination primer scheme files
    :param scheme_name: the name of the primer scheme
    :param pathogen: the name of the pathogen
    :param version_path: the path containing the primer version
    :return: object storing fasta and pickle paths for each primer scheme
    """
    schemes = {}
    version_name = version_path.name
    # e.g. V4.1 => V4-1
    version_name = version_name.replace(".", "-").replace("_", "-")
    dest_scheme_path = dest_schemes_path / scheme_name / SARS_COV_2 / version_name
    dest_scheme_path.mkdir(parents=True, exist_ok=True)
    primer_files = [
        f"{pathogen}.{REFERENCE}.{FASTA}",
        f"{pathogen}.{SCHEME}.{BED}",
        f"{pathogen}.{SCHEME}.{FASTA}",
        f"{pathogen}.{SCHEME}.{PICKLE}",
    ]

    for primer_file in version_path.iterdir():
        if any(primer_file.name.endswith(x) for x in primer_files):
            shutil.copy2(primer_file, dest_scheme_path)

    ref_fasta = dest_scheme_path / f"{pathogen}.{REFERENCE}.{FASTA}"
    scheme_bed = dest_scheme_path / f"{pathogen}.{SCHEME}.{BED}"
    # file to be generated
    scheme_name_version = f"{scheme_name}_{version_name}"
    scheme_fasta_path = dest_scheme_path / f"{pathogen}.{SCHEME}.{FASTA}"
    scheme_pickle_path = dest_scheme_path / f"{pathogen}.{SCHEME}.{PICKLE}"
    extract_primer_sequences(ref_fasta, scheme_bed, scheme_fasta_path, scheme_name_version)
    # create an ahocorasick automaton object so that this is can
    # be loaded quickly for all samples
    automaton = create_automaton(scheme_fasta_path, FASTA)
    # Serialise (pickle) the ahocorasick automaton object
    dump_pickle(scheme_pickle_path, automaton)
    schemes[scheme_name_version] = {
        FASTA: scheme_fasta_path,
        PICKLE: scheme_pickle_path,
    }
    return schemes


def raise_if_not_dir(source_schemes_path: Path) -> None:
    """
    Raise FileNotFoundError if source_schems_path is not a directory

    :param source_scheme_path: the path to the source primer schemes
    """
    if not source_schemes_path.is_dir():
        raise FileNotFoundError("Path to the source scheme was not found. Check input arguments")


def process_sars_cov_2_epi2me_labs_primers(
    base_path: Path, source_schemes: str, dest_schemes_path: Path
) -> SCHEME_TYPE:
    """
    Process SARS-CoV-2 primer schemes in the epi2me-labs repo

    :param base_path: the path to the epi2me-labs cloned repository
    :param source_schemes: the subpath to the primer schemes within the cloned repository
    :param dest_schemes_path: the path storing the destination primer scheme files
    :return: object storing fasta and pickle paths for each primer scheme
    """
    schemes = {}
    source_schemes_path = base_path / source_schemes
    raise_if_not_dir(source_schemes_path)

    for source_scheme_path in source_schemes_path.iterdir():
        scheme_name = source_scheme_path.name
        for version_path in source_scheme_path.iterdir():
            schemes.update(prepare_dest_files(dest_schemes_path, scheme_name, SARS_COV_2, version_path))
    return schemes


def process_sars_cov_2_quick_lab_primers(base_path: Path, source_schemes: str, dest_schemes_path: Path) -> SCHEME_TYPE:
    """
    Process SARS-CoV-2 ARTIC primer schemes in the quick-lab repo

    :param base_path: the path to the quick-lab cloned repository
    :param source_schemes: the subpath to the primer schemes within the cloned repository
    :param dest_schemes_path: the path storing the destination primer scheme files
    :return: object storing fasta and pickle paths for each primer scheme
    """
    schemes = {}
    source_schemes_path = base_path / source_schemes
    raise_if_not_dir(source_schemes_path)

    for scheme_length_dir in source_schemes_path.iterdir():
        # these schemes are separated by length
        if scheme_length_dir.name not in ["400", "1200"]:
            continue
        for source_scheme_path in scheme_length_dir.iterdir():
            # these schemes are all ARTIC
            scheme_name = "ARTIC"
            version = source_scheme_path.name
            # fix name inconsistency (e.g. SARs-CoV-2 which should be SARS-CoV-2)
            for file_path in source_scheme_path.iterdir():
                fixed_filename = re.sub(SARS_COV_2, SARS_COV_2, file_path.name, flags=re.IGNORECASE)
                file_path.rename(Path(file_path.parent, fixed_filename))

            # rename fasta file for consistency with epi2me-labs
            source_fasta = Path(source_scheme_path / f"MN908947.3.{FASTA}")
            source_fasta.rename(Path(source_scheme_path, f"{SARS_COV_2}.{REFERENCE}.{FASTA}"))

            # rename bed file for consistency with epi2me-labs. Also reformat it so that artic minion can process it:
            # - remove primer sequences from bed file
            # - reformat name as artic minion:align_trim() expects: <pathogen>_<pair>_<side>
            #   whereas quick-lab stores the name as: <pathogen>_<length>_<pair>_<side>_<??>
            source_bed = Path(source_scheme_path / f"{SARS_COV_2}_{version}.primer.{BED}")
            dest_bed = Path(source_scheme_path, f"{SARS_COV_2}.{SCHEME}.{BED}")
            with open(source_bed, newline="") as inputfile, open(dest_bed, "w", newline="") as outputfile:
                bed_reader = csv.DictReader(inputfile, fieldnames=BED_COLS, delimiter="\t")
                bed_writer = csv.DictWriter(outputfile, fieldnames=BED_COLS, delimiter="\t")
                for record in bed_reader:
                    # drop the sequence
                    record.pop(None)
                    pathogen, _, pair, side, _ = record[NAME].split("_")
                    record[NAME] = f"{pathogen}_{pair}_{side}"
                    bed_writer.writerow(record)
            source_bed.unlink()

            schemes.update(prepare_dest_files(dest_schemes_path, scheme_name, SARS_COV_2, source_scheme_path))
    return schemes


def generate_primer_index_file(pathogen: str, dest_schemes_path: Path, schemes: SCHEME_TYPE) -> None:
    """
    Generate an index summarising the primer fasta, pickle and primer number

    :param pathogen: the pathogen name
    :param dest_schemes_path: the path storing the destination primer scheme files
    :param schemes: object storing fasta and pickle paths for each primer scheme
    """
    scheme_index_path = dest_schemes_path / f"{pathogen}_primer_index.csv"
    with open(scheme_index_path, "w", newline="") as outputfile:
        writer = csv.DictWriter(outputfile, fieldnames=PRIMER_INDEX_COLS)
        writer.writeheader()
        for scheme in sorted(schemes):
            fasta_file = schemes[scheme][FASTA]
            pickle_file = schemes[scheme][PICKLE]
            fasta_path = Path("/") / Path(fasta_file).relative_to(scheme_index_path.parent.parent)
            pickle_path = Path("/") / Path(pickle_file).relative_to(scheme_index_path.parent.parent)
            num_primers = len(list(SeqIO.parse(fasta_file, FASTA)))
            writer.writerow(
                {
                    PRIMER_NAME: scheme,
                    FASTA_PATH: fasta_path,
                    PICKLE_PATH: pickle_path,
                    TOTAL_NUM_PRIMER: num_primers,
                }
            )
            print(f"Generated primer fasta and pickle files:\n - {fasta_path}\n - {pickle_path}")
    print(f"Generated primer index containing number of primers: {scheme_index_path}")


ORGANISE_PRIMERS = {
    SARS_COV_2: {
        EPI2ME_LABS: process_sars_cov_2_epi2me_labs_primers,
        QUICK_LAB: process_sars_cov_2_quick_lab_primers,
    },
}


def prepare_primers(
    source_name: str,
    primer_repo: str,
    repo_commit: str,
    pathogen: str,
    source_schemes: str,
    dest_schemes_path: Path,
) -> SCHEME_TYPE:
    """
    An entry point function for cloning the repository and invoking the repository-specific function
    responsible for organising the primer data.

    :param source_name: the name of the source (e.g. 'epi2me_labs')
    :param primer_repo: the repository URL storing the primers
    :param repo_commit: the repository commit to checkout
    :param pathogen: the pathogen name
    :param source_schemes: the subpath to the primer schemes within the cloned repository
    :param dest_schemes_path: the path storing the destination primer scheme files
    :return: object storing fasta and pickle paths for each primer scheme
    """
    # use a temporary directory for cloning the repository
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Cloning repo {primer_repo} on {repo_commit}")
        repo_path = Path(tmpdir)
        repo = Repo.clone_from(primer_repo, repo_path)
        repo.git.checkout(repo_commit)
        return ORGANISE_PRIMERS[pathogen][source_name](repo_path, source_schemes, dest_schemes_path)


@click.command()
@click.option(
    "--dependencies-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="The input file containing the primers repositories / commits",
)
@click.option(
    "--dest-schemes",
    type=click.Path(exists=False, file_okay=False, dir_okay=True, writable=True),
    required=True,
    help="output directory to store the primer schemes (e.g. ./primer_schemes)",
)
@click.option(
    "--pathogen",
    type=click.Choice(SUPPORTED_PATHOGENS, case_sensitive=False),
    default=SARS_COV_2,
    help="The name of the pathogen",
)
def fetch_primers(
    dependencies_file: str,
    dest_schemes: str,
    pathogen: str,
) -> None:
    """
    Fetch the SARS-CoV-2 primers and pre-process the primer sequence fasta files.

    e.g.
    python fetch_primers.py \
        --dependencies-file ../../docker/sars_cov_2/primers.sources \
        --dest-schemes ../../data/sars_cov_2/primer_schemes \
        --pathogen sars-cov-2

    :param dependencies_file: the file storing the primer sources to use
    :param dest_schemes: the path storing the destination primer scheme files
    :param pathogen: the pathogen name
    """
    dest_schemes_path = Path(dest_schemes)
    shutil.rmtree(dest_schemes_path, ignore_errors=True)

    schemes = {}

    with open(dependencies_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            source_name = row["source_name"]
            primer_repo = row["repository"]
            repo_commit = row["commit"]
            source_schemes = row["source_schemes"]
            schemes.update(
                prepare_primers(source_name, primer_repo, repo_commit, pathogen, source_schemes, dest_schemes_path)
            )

    generate_primer_index_file(pathogen, dest_schemes_path, schemes)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    fetch_primers()
