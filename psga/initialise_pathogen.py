#!/usr/bin/env python3

from pathlib import Path
import click


def setup_pathogen_config(pathogen_name: str) -> None:
    """
    Create a config file for the pathogen
    """
    config_path = Path(f"{pathogen_name}.config")
    with open(config_path, "w") as fd:
        fd.write(
            "// load the main config\n"
            "includeConfig 'nextflow.config'\n"
            "\n"
            "params {\n"
            "    // mandatory. Directory containing the Nextflow configs and workflows for this pathogen\n"
            f'    pathogen_dir = "{pathogen_name}"\n'
            "}"
        )
    print(f"Generated config file: {config_path}")


def setup_pathogen_dir(pathogen_name: str) -> None:
    """
    Initialise the directory for this pathogen
    """
    pathogen_dir = Path(pathogen_name)

    # create pathogen directory
    pathogen_dir.mkdir(parents=True, exist_ok=True)

    # create pathogen psga.nf file
    with open(pathogen_dir / "psga.nf", "w") as fd:
        fd.write(
            "/*\n"
            f" * Main workflow for the pathogen: {pathogen_name}.\n"
            " */\n"
            "workflow psga {\n"
            "\n"
            "    main:\n"
            "\n"
            "        // to be updated\n"
            "        ch_qc_passed_fasta = Channel.empty()\n"
            "        ch_analysis_run_results_submitted = Channel.empty()\n"
            "\n"
            "    emit:\n"
            "        ch_qc_passed_fasta\n"
            "        ch_analysis_run_results_submitted\n"
            "}"
        )

    # create pathogen help.nf file
    with open(pathogen_dir / "help.nf", "w") as fd:
        fd.write(
            "def printPathogenConfig() {\n"
            "\n"
            '    log.info"""\n'
            f"        {pathogen_name} pathogen\n"
            "        ===================\n"
            "        Environment variables:\n"
            "        * <list pathogen-specific env vars here>\n"
            "\n"
            "        Parameters:\n"
            "        * <list pathogen-specific parameters here>\n"
            '    """.stripIndent()\n'
            "}\n"
            "\n"
            "def printPathogenHelp() {\n"
            '    log.info"""\n'
            f"        {pathogen_name} pathogen\n"
            "          Description:\n"
            "            <describe the pathogen analysis here>\n"
            "\n"
            "          Usage:\n"
            f"            nextflow run . -c {pathogen_name}.config --run [analysis_run] --metadata [csv] [parameters]\n"
            "\n"
            "          Mandatory environment variables:\n"
            "            <any mandatory env vars here>\n"
            "                                    <description>\n"
            "\n"
            "          Mandatory parameters:\n"
            "            <any mandatory parameters here>    <description>\n"
            "\n"
            "          Optional parameters:\n"
            "            <any optional parameters here>    <description>\n"
            '    """.stripIndent()\n'
            "}"
        )
    print(f"Initialised pathogen directory: {pathogen_dir}")


@click.command()
@click.option(
    "--pathogen-name",
    type=str,
    required=True,
    help="The name of the new pathogen to support",
)
def initialise_pathogen(pathogen_name):
    """
    This script initialises the files for analysing a new pathogen using the PSGA pipeline.
    """
    setup_pathogen_config(pathogen_name)
    setup_pathogen_dir(pathogen_name)
    print(f"Pathogen {pathogen_name} has been initialised")
    print("Customise the workflow for this pathogen adding / modifying files in the following paths:")
    print(f" * nextflow configs and workflows: ./{pathogen_name}")
    print(f" * Python scripts: ../scripts/{pathogen_name}")
    print(" * db schema: ../scripts/db")
    print(f" * Python tests: ../tests/{pathogen_name}")


if __name__ == "__main__":
    initialise_pathogen()  # pylint: disable=no-value-for-parameter,unexpected-keyword-arg