import pytest

from click.testing import CliRunner

from scripts.db.models import Sample
from scripts.load_nextstrain_nt_muts_to_db import load_nextstrain_nt_muts_data


@pytest.mark.parametrize(
    "nt_muts_json,tree_nwk,expected_data",
    [
        (
            "nt_muts.json",
            "tree.nwk",
            {
                "NC_045512.2": "",
                "12704502": "C241T,C1059T,C3037T,C4423T,A6105G,T8041C,C11967T,C14408T,"
                "G22225T,A23403G,G25563T,C28005T",
                "12704501": "C241T,C1059T,C3037T,C3168A,T3654A,C12885T,C14408T,C16260T,"
                "C22987T,A23403G,C24370T,T24721C,G25563T,C27476T,A28113C,C28821A",
            },
        )
    ],
)
def test_load_nextstrain_nt_muts_data(db_session, test_data_path, nt_muts_json, tree_nwk, expected_data):
    rv = CliRunner().invoke(
        load_nextstrain_nt_muts_data,
        [
            "--nt-muts-json",
            test_data_path / nt_muts_json,
            "--tree-nwk",
            test_data_path / tree_nwk,
        ],
    )

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()

    assert len(samples) == 2

    for sample in samples:
        # elements in strings are sorted in the script
        assert sample.nucleotide_muts == expected_data[sample.lab_id]
