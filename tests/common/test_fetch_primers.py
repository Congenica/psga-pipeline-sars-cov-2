from pathlib import Path
import shutil
import pytest
from Bio import SeqIO

from app.scripts.fetch_primers import (
    extract_primer_sequences,
    generate_primer_index_file,
    ORGANISE_PRIMERS,
    SARS_COV_2,
    EPI2ME_LABS,
    QUICK_LAB,
    BED,
    FASTA,
    PICKLE,
    REFERENCE,
    SCHEME,
)
from app.scripts.primer_autodetection import load_pickle
from tests.utils_tests import assert_files_are_equal

PATH = "path"
SOURCE_SCHEMES = "source_schemes"
INDEX = "index"
SOURCE_REFERENCE = f"source_{REFERENCE}"
DEST_SCHEME_FASTA = f"dest_{SCHEME}_{FASTA}"
DEST_SCHEME_BED = f"dest_{SCHEME}_{BED}"
DEST_SCHEME_PICKLE = f"dest_{SCHEME}_{PICKLE}"
PRIMER_COUNT = "primer_count"
DEST_SCHEME_BED_FILENAME = f"{SARS_COV_2}.{SCHEME}.{BED}"
DEST_SCHEME_FASTA_FILENAME = f"{SARS_COV_2}.{SCHEME}.{FASTA}"
DEST_SCHEME_PICKLE_FILENAME = f"{SARS_COV_2}.{SCHEME}.{PICKLE}"
INDEX_FILENAME = f"{SARS_COV_2}_primer_{INDEX}.csv"


@pytest.fixture
def primer_data(fetch_primers_data_path: Path) -> dict:
    return {
        EPI2ME_LABS: {
            PATH: fetch_primers_data_path / EPI2ME_LABS,
            SOURCE_SCHEMES: f"data/primer_schemes/{SARS_COV_2}",
            SOURCE_REFERENCE: f"{SARS_COV_2}.{REFERENCE}.{FASTA}",
            PRIMER_COUNT: 4,
            DEST_SCHEME_BED: DEST_SCHEME_BED_FILENAME,
            DEST_SCHEME_FASTA: DEST_SCHEME_FASTA_FILENAME,
            DEST_SCHEME_PICKLE: DEST_SCHEME_PICKLE_FILENAME,
            INDEX: INDEX_FILENAME,
        },
        QUICK_LAB: {
            PATH: fetch_primers_data_path / QUICK_LAB,
            SOURCE_SCHEMES: ".",
            SOURCE_REFERENCE: f"MN908947.3.{FASTA}",
            PRIMER_COUNT: 4,
            DEST_SCHEME_BED: DEST_SCHEME_BED_FILENAME,
            DEST_SCHEME_FASTA: DEST_SCHEME_FASTA_FILENAME,
            DEST_SCHEME_PICKLE: DEST_SCHEME_PICKLE_FILENAME,
            INDEX: INDEX_FILENAME,
        },
    }


def get_file_paths(my_path: Path, pattern: str) -> list[Path]:
    return sorted(list(my_path.rglob(pattern)))


def get_pickle_content(pickle_path: Path) -> list[tuple]:
    return [(index, sequence) for index, sequence in enumerate(load_pickle(pickle_path))]


def assert_pickle_files_are_equal(pickle_path_1: Path, pickle_path_2: Path) -> None:
    assert get_pickle_content(pickle_path_1) == get_pickle_content(pickle_path_2)


def assert_primer_files(dest_schemes_path: Path, test_data: dict, primer_file_key: str) -> None:
    calc_scheme_files = get_file_paths(dest_schemes_path, test_data[primer_file_key])
    expected_scheme_files = get_file_paths(test_data[PATH], test_data[primer_file_key])

    assert len(calc_scheme_files) == len(expected_scheme_files)

    for calc, exp in zip(sorted(calc_scheme_files), sorted(expected_scheme_files)):
        if primer_file_key == DEST_SCHEME_PICKLE:
            assert_pickle_files_are_equal(calc, exp)
        elif primer_file_key in [DEST_SCHEME_FASTA, DEST_SCHEME_BED, INDEX]:
            assert_files_are_equal(calc, exp)
        else:
            raise KeyError(f"assert_primer_files() does not support key: {primer_file_key}")


def rmfile(filepath: Path, filename: str) -> None:
    """
    remove a file called fname in f path
    """
    if filepath.is_file():
        if filepath.name == filename:
            filepath.unlink()
    else:
        for child in filepath.iterdir():
            rmfile(child, filename)


@pytest.mark.parametrize(
    "data_source,primer,primer_subdir,source_bed",
    [
        (EPI2ME_LABS, "ARTIC_V4", "ARTIC/V4", f"{SARS_COV_2}.{SCHEME}.{BED}"),
        (EPI2ME_LABS, "Midnight-ONT_V2", "Midnight-ONT/V2", f"{SARS_COV_2}.{SCHEME}.{BED}"),
        (QUICK_LAB, "ARTIC_v5-2-0-1200", "1200/v5.2.0_1200", f"{SARS_COV_2}_v5.2.0_1200.primer.{BED}"),
        (QUICK_LAB, "ARTIC_v5-3-2-400", "400/v5.3.2_400", f"SARs-CoV-2_v5.3.2_400.primer.{BED}"),
    ],
    ids=[
        "epi2me-labs, artic v4",
        "epi2me-labs, midnight-ont v2",
        "quick-lab, artic v5.2.0-1200",
        "quick-lab, artic v5.3.2-400 (source primer bed with mispelled name)",
    ],
)
@pytest.mark.jira(identifier="337f1ac9-8b1c-4ccd-ae61-b7a7fa482b4a", confirms="PSG-3621")
def test_extract_primer_sequences(
    tmp_path: Path,
    primer_data: dict,
    data_source: str,
    primer: str,
    primer_subdir: str,
    source_bed: str,
):
    # NOTE: for quick-lab the source bed file contains the version in the file name
    # It also shows name inconsistency.
    # The script handles that, but here an internal function is called, so the original name is specified
    data = primer_data[data_source]
    primer_dir = data[PATH] / data[SOURCE_SCHEMES] / primer_subdir
    ref_fasta_path = primer_dir / data[SOURCE_REFERENCE]
    scheme_bed_path = primer_dir / source_bed

    expected_scheme_fasta_path = primer_dir / data[DEST_SCHEME_FASTA]
    scheme_fasta_path = tmp_path / data[DEST_SCHEME_FASTA]

    extract_primer_sequences(ref_fasta_path, scheme_bed_path, scheme_fasta_path, primer)

    primer_counter = len(list(SeqIO.parse(scheme_fasta_path, FASTA)))
    assert primer_counter == data[PRIMER_COUNT]
    assert_files_are_equal(scheme_fasta_path, expected_scheme_fasta_path)


@pytest.mark.jira(identifier="ee8f4c43-c5d4-4904-8ef3-ef72035b5c5d", confirms="PSG-3621")
def test_extract_primer_sequences_unknown_strand_error(
    tmp_path: Path,
    fetch_primers_data_path: Path,
):
    primer = "NO_STRAND"
    ref_fasta_path = fetch_primers_data_path / primer / f"{SARS_COV_2}.{REFERENCE}.{FASTA}"
    scheme_bed_path = fetch_primers_data_path / primer / f"{SARS_COV_2}.{SCHEME}.{BED}"
    scheme_fasta_path = tmp_path / f"{SARS_COV_2}.{FASTA}"

    with pytest.raises(ValueError, match="Unknown strand"):
        extract_primer_sequences(ref_fasta_path, scheme_bed_path, scheme_fasta_path, primer)


@pytest.mark.parametrize(
    "data_source,primer,version",
    [
        (EPI2ME_LABS, "ARTIC", "V4"),
        (EPI2ME_LABS, "Midnight-ONT", "V2"),
    ],
)
@pytest.mark.parametrize(
    "reference,scheme_bed",
    [
        (True, False),
        (False, True),
    ],
    ids=[
        f"Missing {SCHEME} {BED} file",
        f"Missing {REFERENCE} {FASTA} file",
    ],
)
@pytest.mark.jira(identifier="b048656c-163a-4aec-bc76-83ad39491d12", confirms="PSG-3621")
def test_extract_primer_sequences_file_not_found_error(
    tmp_path: Path,
    primer_data: dict,
    data_source: str,
    primer: str,
    version: str,
    reference: bool,
    scheme_bed: bool,
):
    data = primer_data[data_source]
    primer_dir = data[PATH] / data[SOURCE_SCHEMES] / primer / version
    ref_fasta_path = primer_dir / data[SOURCE_REFERENCE] if reference else None
    scheme_bed_path = primer_dir / data[DEST_SCHEME_BED] if scheme_bed else None
    scheme_fasta_path = tmp_path / data[DEST_SCHEME_FASTA]

    # assert that the other one exists, making sure both FileNotFoundError exceptions are raised
    if ref_fasta_path:
        assert ref_fasta_path.is_file()
    if scheme_bed_path:
        assert scheme_bed_path.is_file()

    with pytest.raises(FileNotFoundError, match="was not found"):
        extract_primer_sequences(ref_fasta_path, scheme_bed_path, scheme_fasta_path, primer)


@pytest.mark.parametrize(
    "data_source",
    [
        EPI2ME_LABS,
        QUICK_LAB,
    ],
)
@pytest.mark.jira(identifier="7ba24142-060d-4b50-9100-0243d92c4231", confirms="PSG-3621")
def test_organise_primers(
    tmp_path: Path,
    primer_schemes_dir: str,
    primer_data: dict,
    data_source: str,
):
    data = primer_data[data_source]
    dest_schemes_path = tmp_path / primer_schemes_dir

    # copy the input primers, mocking repository cloning
    shutil.copytree(data[PATH], tmp_path, dirs_exist_ok=True)

    schemes = ORGANISE_PRIMERS[SARS_COV_2][data_source](
        tmp_path,
        data[SOURCE_SCHEMES],
        dest_schemes_path,
    )

    # assert scheme files are as expected
    assert_primer_files(dest_schemes_path, data, DEST_SCHEME_FASTA)
    assert_primer_files(dest_schemes_path, data, DEST_SCHEME_BED)
    assert_primer_files(dest_schemes_path, data, DEST_SCHEME_PICKLE)

    # assert generated index file is as expected
    generate_primer_index_file(SARS_COV_2, dest_schemes_path, schemes)
    assert_primer_files(dest_schemes_path, data, INDEX)
