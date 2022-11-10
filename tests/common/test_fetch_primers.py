from pathlib import Path
import shutil
import pytest
from Bio import SeqIO

from scripts.common.fetch_primers import (
    get_version,
    extract_primer_sequences,
    ORGANISE_PRIMERS,
    SARS_COV_2,
)
from tests.utils_tests import assert_files_are_equal

FETCH_PRIMERS_DIR = "fetch_primers"
PRIMER_SCHEMES = "primer_schemes"


def rmfile(f: Path, fname: str) -> None:
    """
    remove a file called fname in f path
    """
    if f.is_file():
        if f.name == fname:
            f.unlink()
    else:
        for child in f.iterdir():
            rmfile(child, fname)


@pytest.mark.parametrize(
    "dependency_file,key,expected_value",
    [
        ("dependencies", "primer_schemes_commit", "a60a1e1e73bde1971c680bd3d53076127dd63fc6"),
    ],
)
def test_get_version(test_data_path, dependency_file, key, expected_value):
    dependency_path = test_data_path / FETCH_PRIMERS_DIR / dependency_file
    assert expected_value == get_version(dependency_path, key)


@pytest.mark.parametrize("pathogen", [SARS_COV_2])
@pytest.mark.parametrize(
    "primer,version",
    [
        ("ARTIC", "V4"),
        ("Midnight-ONT", "V2"),
    ],
)
@pytest.mark.parametrize(
    "ref_fasta,scheme_bed,scheme_fasta,expected_primer_counter",
    [
        (".reference.fasta", ".scheme.bed", ".scheme.fasta", 8),
    ],
)
def test_extract_primer_sequences(
    tmp_path, test_data_path, pathogen, primer, version, expected_primer_counter, ref_fasta, scheme_bed, scheme_fasta
):
    primer_dir = test_data_path / FETCH_PRIMERS_DIR / PRIMER_SCHEMES / pathogen / primer / version
    ref_fasta_path = primer_dir / f"{pathogen}{ref_fasta}"
    scheme_bed_path = primer_dir / f"{pathogen}{scheme_bed}"
    expected_scheme_fasta_path = primer_dir / f"{pathogen}{scheme_fasta}"
    scheme_fasta_path = tmp_path / f"{pathogen}{scheme_fasta}"

    extract_primer_sequences(ref_fasta_path, scheme_bed_path, scheme_fasta_path, f"{primer}_{version}")

    primer_counter = len(list(SeqIO.parse(scheme_fasta_path, "fasta")))
    assert primer_counter == expected_primer_counter
    assert_files_are_equal(scheme_fasta_path, expected_scheme_fasta_path)


@pytest.mark.parametrize("pathogen", [SARS_COV_2])
@pytest.mark.parametrize("primer", ["NO_STRAND"])
@pytest.mark.parametrize(
    "ref_fasta,scheme_bed",
    [
        (".reference.fasta", ".scheme.bed"),
    ],
)
def test_extract_primer_sequences_unknown_strand_error(
    tmp_path, test_data_path, pathogen, primer, ref_fasta, scheme_bed
):
    ref_fasta_path = test_data_path / FETCH_PRIMERS_DIR / primer / f"{pathogen}{ref_fasta}"
    scheme_bed_path = test_data_path / FETCH_PRIMERS_DIR / primer / f"{pathogen}{scheme_bed}"
    scheme_fasta_path = tmp_path / f"{pathogen}.fasta"

    with pytest.raises(ValueError, match="Unknown strand"):
        extract_primer_sequences(ref_fasta_path, scheme_bed_path, scheme_fasta_path, primer)


@pytest.mark.parametrize("pathogen", [SARS_COV_2])
@pytest.mark.parametrize(
    "primer,version",
    [
        ("ARTIC", "V4"),
        ("Midnight-ONT", "V2"),
    ],
)
@pytest.mark.parametrize(
    "ref_fasta,scheme_bed",
    [
        (".reference.fasta", None),
        (None, ".scheme.bed"),
    ],
)
def test_extract_primer_sequences_file_not_found_error(
    tmp_path, test_data_path, pathogen, primer, version, ref_fasta, scheme_bed
):
    primer_dir = test_data_path / FETCH_PRIMERS_DIR / PRIMER_SCHEMES / pathogen / primer / version
    ref_fasta_path = primer_dir / f"{pathogen}{ref_fasta}" if ref_fasta else None
    scheme_bed_path = primer_dir / f"{pathogen}{scheme_bed}" if scheme_bed else None
    scheme_fasta_path = tmp_path / f"{pathogen}.fasta"

    # assert that the other one exists, making sure both FileNotFoundError exceptions are raised
    if ref_fasta_path:
        assert ref_fasta_path.is_file()
    if scheme_bed_path:
        assert scheme_bed_path.is_file()

    with pytest.raises(FileNotFoundError, match="was not found"):
        extract_primer_sequences(ref_fasta_path, scheme_bed_path, scheme_fasta_path, primer)


@pytest.mark.parametrize(
    "pathogen,index,scheme_fasta",
    [
        (
            SARS_COV_2,
            f"{SARS_COV_2}_primer_fasta_index.txt",
            f"{SARS_COV_2}.scheme.fasta",
        ),
    ],
)
def test_organise_primers(tmp_path, test_data_path, pathogen, index, scheme_fasta):
    source_scheme_path = test_data_path / FETCH_PRIMERS_DIR / PRIMER_SCHEMES
    source_scheme_tmp_path = tmp_path / f"source_{PRIMER_SCHEMES}"

    # copy the input primers and delete the expected scheme.fasta
    shutil.copytree(source_scheme_path, source_scheme_tmp_path)
    rmfile(source_scheme_tmp_path, scheme_fasta)

    source_scheme_tmp_path = Path(source_scheme_tmp_path) / pathogen
    dest_scheme_path = tmp_path / PRIMER_SCHEMES

    ORGANISE_PRIMERS[pathogen](source_scheme_tmp_path, dest_scheme_path)

    # assert index file
    assert_files_are_equal(
        source_scheme_path / index,
        dest_scheme_path / index,
    )

    # assert scheme fasta files
    def get_scheme_fasta(p: Path):
        return sorted(list(p.rglob(scheme_fasta)))

    calc_scheme_fasta = get_scheme_fasta(dest_scheme_path)
    exp_scheme_fasta = get_scheme_fasta(source_scheme_path)
    assert len(calc_scheme_fasta) > 0
    for calc, exp in zip(calc_scheme_fasta, exp_scheme_fasta):
        assert_files_are_equal(calc, exp)
