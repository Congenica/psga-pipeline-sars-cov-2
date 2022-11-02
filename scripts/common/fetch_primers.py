import shutil
from pathlib import Path
from typing import List
from collections import Counter
import csv
import tempfile
import click
from git import Repo
from Bio import SeqIO
from Bio.Seq import Seq

SARS_COV_2 = "SARS-CoV-2"
SUPPORTED_PATHOGENS = [SARS_COV_2]
START, END, NAME, STRAND = "start", "end", "name", "strand"
BED_COLS = [
    "scaffold",
    START,
    END,
    NAME,
    "score",
    STRAND,
]
PRIMER_KEY, STRAND_KEY = "primers", "strand"


def _decorate_with_new_line(method):
    """
    Automatically add new line
    """

    def decorated(text):
        method(f"{text}\n")

    return decorated


def get_version(deps_path: Path, key: str) -> str:
    """
    Return the version for the key in deps_path
    """
    with open(deps_path, newline="") as inputfile:
        reader = csv.DictReader(inputfile, fieldnames=["key", "value"], delimiter="=")
        for row in reader:
            if row["key"] == key:
                return row["value"]

    raise ValueError(f"Key: {key} was not found in {deps_path}")


def extract_primer_sequences(
    ref_fasta: Path,
    scheme_bed: Path,
    scheme_fasta: Path,
    primer: str,
    strands_in_name: List[str] = None,
    right_strand: str = None,
) -> Counter:
    """
    Extract primer sequences from ref_fasta using coordinates in scheme_bed
    and generate a scheme_fasta containing these sequences.

    Some scheme bed files do not have `strand`.
    In these cases, the strand is obtained from the sequence name (e.g. _LEFT, _RIGHT).
    An error is raised if the strand cannot be inferred
    """

    if not strands_in_name:
        strands_in_name = ["LEFT", "RIGHT"]
        right_strand = "RIGHT"

    if not (ref_fasta and ref_fasta.is_file()):
        raise FileNotFoundError(f"File {ref_fasta} was not found")

    if not (scheme_bed and scheme_bed.is_file()):
        raise FileNotFoundError(f"File {scheme_bed} was not found")

    # load the reference sequence
    # the reference fasta contains 1 single record
    for seq_record in SeqIO.parse(ref_fasta, "fasta"):
        ref_sequence = str(seq_record.seq)

    # count the number of primers and strands as a basic QC
    primer_counter = Counter({PRIMER_KEY: 0, STRAND_KEY: 0})

    with open(scheme_bed, newline="") as bedfile, open(scheme_fasta, "w") as outputfile:
        write = _decorate_with_new_line(outputfile.write)
        bed_reader = csv.DictReader(bedfile, fieldnames=BED_COLS, delimiter="\t")
        write(f">{primer}")
        for bed_row in bed_reader:
            primer_counter.update([PRIMER_KEY])
            start = int(bed_row[START])
            end = int(bed_row[END])
            primer_sequence = ref_sequence[start:end]
            if bed_row[STRAND]:
                if bed_row[STRAND] == "-":
                    # RIGHT strand
                    primer_sequence = Seq(primer_sequence).reverse_complement()
                    primer_counter.update([STRAND_KEY])
            else:
                # bed_row[STRAND] is None, infer the strand from name if possible or raise error
                try:
                    inferred_strand = next(s for s in strands_in_name if s in bed_row[NAME])
                    if inferred_strand == right_strand:
                        primer_sequence = Seq(primer_sequence).reverse_complement()
                        primer_counter.update([STRAND_KEY])
                except StopIteration as unknown_strand:
                    raise ValueError(
                        f"Unknown strand in {scheme_bed}. Cannot extract primer sequences"
                    ) from unknown_strand

            write(primer_sequence)

    return primer_counter


def organise_sars_cov_2_primers(source_schemes_path: Path, dest_schemes_path: Path) -> None:
    """
    Organise the primer schemes in order to execute `artic minion` correctly to:
    dest_schemes_path/<SCHEME_NAME>/<PATHOGEN>/<SCHEME_VERSION>

    This function also adds the fasta file of the primer sequences and their index file
    """
    print("Copy primers:")
    scheme_fastas = []
    for source_scheme_path in source_schemes_path.iterdir():
        scheme_name = source_scheme_path.name
        for version_path in source_scheme_path.iterdir():
            version_name = version_path.name
            # e.g. V4.1 => V4
            version_name = version_name.replace(".", "-")
            dest_scheme_path = dest_schemes_path / scheme_name / SARS_COV_2 / version_name
            dest_scheme_path.mkdir(parents=True, exist_ok=True)
            for primer_file in version_path.iterdir():
                shutil.copy2(primer_file, dest_scheme_path)

            ref_fasta = dest_scheme_path / f"{SARS_COV_2}.reference.fasta"
            scheme_bed = dest_scheme_path / f"{SARS_COV_2}.scheme.bed"
            # file to be generated
            scheme_fasta = dest_scheme_path / f"{SARS_COV_2}.scheme.fasta"
            primer_counter = extract_primer_sequences(
                ref_fasta, scheme_bed, scheme_fasta, f"{scheme_name}_{version_name}"
            )
            print(
                f"\t{version_path} to {dest_scheme_path}; "
                f"primers: {primer_counter[PRIMER_KEY]}; "
                f"reverse complemented primers: {primer_counter[STRAND_KEY]}"
            )
            scheme_fastas.append(scheme_fasta)

    scheme_fasta_index_path = dest_schemes_path / f"{SARS_COV_2}_primer_fasta_index.txt"
    with open(scheme_fasta_index_path, "w") as outputfile:
        write = _decorate_with_new_line(outputfile.write)
        for fasta in sorted(scheme_fastas):
            fasta_path = Path("/") / Path(fasta).relative_to(scheme_fasta_index_path.parent.parent)
            write(fasta_path)


ORGANISE_PRIMERS = {
    SARS_COV_2: organise_sars_cov_2_primers,
}


@click.command()
@click.option(
    "--dependencies-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="The input file containing the primer_scheme_commit",
)
@click.option(
    "--dest-schemes",
    type=click.Path(exists=False, file_okay=False, dir_okay=True, writable=True),
    required=True,
    help="output directory to store the primer schemes (e.g. ./primer_schemes)",
)
@click.option(
    "--primer-repo",
    type=str,
    default="https://github.com/epi2me-labs/wf-artic",
    help="The repository containing the primer schemes",
)
@click.option(
    "--primer-repo-commit-key",
    type=str,
    default="primer_schemes_commit",
    help="The key name in --ncov-deps file containing the commit to checkout",
)
@click.option(
    "--source-schemes",
    type=str,
    default="data/primer_schemes",
    help="The path to the primer schemes in the github repository",
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
    source_schemes: str,
    primer_repo: str,
    primer_repo_commit_key: str,
    pathogen: str,
) -> None:
    """
    Fetch the SARS-CoV-2 primers and pre-process the primer sequence fasta files.

    e.g.
    python fetch_primers.py \
        --dependencies-file ../../docker/sars_cov_2/ncov.deps \
        --dest-schemes ../../docker/primer_schemes \
        --pathogen sars-cov-2
    """
    dest_schemes_path = Path(dest_schemes)
    shutil.rmtree(dest_schemes_path, ignore_errors=True)

    repo_commit = get_version(Path(dependencies_file), primer_repo_commit_key)

    # use a temporary directory for cloning the repository
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Cloning repo {primer_repo} on {repo_commit}")
        tmpdir_path = Path(tmpdir)
        repo = Repo.clone_from(primer_repo, tmpdir_path)
        repo.git.checkout(repo_commit)

        source_schemes_path = tmpdir_path / source_schemes / pathogen

        if not source_schemes_path.is_dir():
            raise FileNotFoundError("Path to the source scheme was not found. Check input arguments")

        ORGANISE_PRIMERS[pathogen](source_schemes_path, dest_schemes_path)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    fetch_primers()
