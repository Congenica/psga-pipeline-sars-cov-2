import pytest

from integration_tests.compare import compare_output_files_set, ValidationError


@pytest.mark.parametrize(
    "calc_output_files,exp_output_files,exc",
    [
        (
            set(),
            set(),
            None,
        ),
        (
            {"a", "b", "c"},
            {"a", "b", "c"},
            None,
        ),
        (
            {"a/b/c1", "a/b/c2", "d/e", "f"},
            {"a/b/c1", "a/b/c2", "d/e", "f"},
            None,
        ),
        (
            {"a/b/c1", "a/b/c2", "d/e", "f"},
            {"a/b/c1", "a/b/c2", "d/e"},
            None,
        ),
        (
            {"a/b/c1", "a/b/c2", "d/e"},
            {"a/b/c1", "a/b/c2", "d/e", "f"},
            "Validation FAILED. See above for details.",
        ),
    ],
)
@pytest.mark.jira(identifier="73dcd85e-163a-4e2b-8644-cfbc4d263abb", confirms="PSG-3621")
def test_compare_output_files_set(calc_output_files: set[str], exp_output_files: set[str], exc: str):
    if exc:
        with pytest.raises(ValidationError, match=exc):
            compare_output_files_set(calc_output_files, exp_output_files)
    else:
        compare_output_files_set(calc_output_files, exp_output_files)
