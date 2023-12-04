from csv import DictWriter
import pathlib
import pytest
from click.testing import CliRunner

from scripts.common.reference_data import ReferenceDataNotFoundError, get_location


class TestReferenceData:

    @pytest.fixture(autouse=True)
    def setup_cli_runner(self):
        self.cli_runner = CliRunner()

    @pytest.fixture(autouse=True)
    def setup_valid_csv(self, tmpdir):
        headers = ["NAME", "LOCATION"]
        rows = [
            {"NAME": "ref-A", "LOCATION": "path/ref-A"},
            {"NAME": "ref-B", "LOCATION": "path/ref-B"},
        ]
        self.reference_data_csv = tmpdir.join("reference_data.csv")
        self.create_reference_data_csv(
            headers=headers,
            rows=rows
        )

    def create_reference_data_csv(self, headers: list[str], rows: list[dict[str, str]]) -> None:
        with self.reference_data_csv.open("w") as fd:
            writer = DictWriter(fd, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)

    @pytest.mark.jira(identifier="4bbf75c9-f54a-4585-9b86-61931eba22da", confirms="PSG-4792")
    def test_get_location_with_valid_csv_and_name_check_output(self):
        result = self.cli_runner.invoke(get_location, [
            str(self.reference_data_csv),
            "ref-A"
        ])
        assert result.exit_code == 0
        assert result.output == f"{pathlib.Path('path/ref-A').absolute()}\n"

    @pytest.mark.jira(identifier="f3076624-1740-438e-9622-b5f420668617", confirms="PSG-4792")
    def test_get_location_with_csv_mixed_case_check_reference_found(self, tmpdir):
        headers = ["nAmE", "lOcAtIoN"]
        rows = [{"nAmE": "ref-A", "lOcAtIoN": "path/ref-A"}]
        self.create_reference_data_csv(headers=headers, rows=rows)

        result = self.cli_runner.invoke(get_location, [
            str(self.reference_data_csv),
            "ref-A"
        ])
        assert result.exit_code == 0
        assert result.output == f"{pathlib.Path('path/ref-A').absolute()}\n"

    @pytest.mark.jira(identifier="149dee28-ac1f-4922-9b24-363db78c3a6d", confirms="PSG-4792")
    def test_get_location_with_empty_csv_check_exception(self):
        self.create_reference_data_csv(headers=["NAME", "LOCATION"], rows=[])

        result = self.cli_runner.invoke(get_location, [
            str(self.reference_data_csv),
            "ref-A"
        ])
        assert result.exit_code == 1
        assert isinstance(result.exception, ReferenceDataNotFoundError)

    @pytest.mark.jira(identifier="dd667211-7cc2-4b70-9840-093e5ad687ba", confirms="PSG-4792")
    def test_get_location_with_name_not_found_check_exception(self):
        result = self.cli_runner.invoke(get_location, [
            str(self.reference_data_csv),
            "NOT-PRESENT"
        ])
        assert result.exit_code == 1
        assert isinstance(result.exception, ReferenceDataNotFoundError)
