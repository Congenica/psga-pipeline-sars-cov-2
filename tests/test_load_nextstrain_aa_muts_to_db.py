import pytest

from click.testing import CliRunner

from scripts.db.models import Sample
from scripts.load_nextstrain_aa_muts_to_db import load_nextstrain_aa_muts_data


@pytest.mark.parametrize(
    "aa_muts_json,tree_nwk,expected_data",
    [
        (
            "aa_muts.json",
            "tree.nwk",
            {
                "NC_045512.2": "",
                "12704502": "ORF1a:T265I,Y1947C,A3901V;ORF1b:P314L;ORF3a:Q57H;ORF8:P38S;S:D614G",
                "12704501": "N:S183Y;ORF1a:T265I,A968D,L1130H,T4207I;ORF1b:P314L;"
                "ORF3a:Q57H;ORF7a:T28I;ORF8:I74L;ORF14:L30I;S:D614G",
            },
        )
    ],
)
def test_load_nextstrain_aa_muts_data(db_session, test_data_path, aa_muts_json, tree_nwk, expected_data):
    rv = CliRunner().invoke(
        load_nextstrain_aa_muts_data,
        [
            "--aa-muts-json",
            test_data_path / aa_muts_json,
            "--tree-nwk",
            test_data_path / tree_nwk,
        ],
    )

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()

    assert len(samples) == 3

    for sample in samples:
        # elements in strings are sorted in the script
        assert sample.amino_acid_muts == expected_data[sample.lab_id]
