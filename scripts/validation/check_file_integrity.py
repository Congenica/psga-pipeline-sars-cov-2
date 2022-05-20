#!/usr/bin/env python

from pathlib import Path
import hashlib
import click

BLOCKSIZE = 65536


class FileIntegrityError(Exception):
    def __eq__(self, other):
        if isinstance(other, FileIntegrityError):
            return str(self) == str(other)
        return False


def md5sum(input_path: Path, blocksize: int = BLOCKSIZE) -> str:
    hashfun = hashlib.md5()
    with open(input_path, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hashfun.update(block)
    return hashfun.hexdigest()


def validate(input_path: Path, expected_md5: str, blocksize: int) -> None:
    computed_md5 = md5sum(input_path, blocksize)
    if computed_md5 != expected_md5:
        raise FileIntegrityError(
            f"Integrity check for file: {input_path.name} FAILED. Expected: {expected_md5}, computed: {computed_md5}"
        )
    print(f"File {input_path.name} is OK")


@click.command()
@click.option(
    "--input-path", required=True, type=click.Path(exists=True, readable=True), help="The input file to validate"
)
@click.option("--expected-md5", required=True, type=str, help="The expected md5 for the input file")
@click.option("--blocksize", default=BLOCKSIZE, type=int, help="The block size used for generating the file md5")
def check_file_integrity(
    input_path,
    expected_md5,
    blocksize,
):
    """
    Check the file integrity using the expected md5.
    """
    validate(Path(input_path), expected_md5, blocksize)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    check_file_integrity()
