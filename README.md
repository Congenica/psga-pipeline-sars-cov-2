# Covid-pipeline documentation

## Operation

### Environment variables

Environment variables required to run the pipeline.
These environment variables are set in `${COVID_PIPELINE_ROOTDIR}/covid-pipeline/nextflow.config` .

| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| DB_HOST | Postgres database host address (e.g. 192.168.0.33) |
| DB_NAME | Postgres database name (e.g. bahrain_sars_cov_2) |
| DB_USER | Postgres database user name (e.g. postgres) |
| DB_PASSWORD | Postgres database user password (e.g. postgres) |
| COVID_PIPELINE_ROOTDIR | Path to the pipeline code (e.g. git checkout). Default: ${HOME}/ps-bahrain-covid |
| COVID_PIPELINE_FASTQ_PATH | Path to the input FASTQ files and TSV metadata file. Default: ${HOME}/Bahrain_COVID_s3_data_lite/sample_data |
| COVID_PIPELINE_WORKDIR | Path to the whole pipeline output. Default: ${HOME}/covid-pipeline |
| COVID_PIPELINE_REPORTS_PATH | Path to the pipeline reports. Default: ${HOME}/reports |


The following environment variables are set internally and should not be changed
| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| COVID_PIPELINE_MISSING_METADATA_PATH | Path to the missing metadata files. Set to: ${COVID_PIPELINE_WORKDIR}/no-metadata-found-fastq |
| COVID_PIPELINE_FASTA_PATH | Path to the re-headered ncov FASTA files. Set to: ${COVID_PIPELINE_WORKDIR}/reheadered-fasta |
| COVID_PIPELINE_FASTA_PATH_QC_FAILED | Path to the re-headered ncov QC_FAILED FASTA files. Set to: ${COVID_PIPELINE_WORKDIR}/reheadered-fasta-qc-failed |


### Dependencies

```commandline
# install nextflow (copy nextflow executable to a directory in your PATH environment variable):
wget -qO- https://get.nextflow.io | bash

# build the covid-pipeline image
docker build -t covid-pipeline:1.0.0 .

# add project submodules
git submodule init
git submodule update

# build ncov docker image
docker build -t ncov2019_artic_nf:1.0.0 -f Dockerfile.ncov2019-artic-nf .

# build pangolin docker image
docker build -t pangolin:1.0.0 -f Dockerfile.pangolin .

# build nextstrain docker image
docker build -t nextstrain:1.0.0 -f Dockerfile.nextstrain .

# build auspice docker image (for visualising the results with Auspice web-service):
docker build -t auspice:1.0.0 -f Dockerfile.auspice .
```

### Run covid-pipeline

The following command runs the covid-pipeline
```commandline
cd ${COVID_PIPELINE_ROOTDIR}/covid-pipeline
nextflow run .
```

The following command cleans up the previous run's work directories and cache, but retains the content of ${COVID_PIPELINE_WORKDIR}:
```commandline
nextflow clean -f
```

### Analyse the results via Auspice web-service
Once the Nextstrain process of the Covid-Pipeline has completed, the results can be visualised as follows.

Run the auspice web-service:
```commandline
# the port 4000 is already exposed in the Dockerfile
docker run -it --rm \
  -v ${COVID_PIPELINE_WORKDIR}/nextstrain_output/bahrain/ncov_with_accessions.json:/ncov_with_accessions.json \
  -p 4000:4000 \
  auspice:1.0.0 \
  auspice view --datasetDir=/
```

Start your local browser from the same machine running the auspice docker image using the address:
```commandline
google-chrome http://localhost:4000
```
In the section `Available datasets`, click on `ncov/with/accessions`


## Development

### Install dependency packages using Python Poetry tool
`Poetry` manages Python dependencies. Dependencies are declared in `pyproject.toml` and exact versions of both dependencies and sub-dependencies are stored in `poetry.lock`. Both are committed to the git repo.

#### Setting up your local environment
Follow the [instructions](https://python-poetry.org/docs/) to install Poetry depending on your local environment. Then run:
```commandline
poetry install
```
in the git checkout to create a virtual environment and install all dependencies. A [plugin](https://plugins.jetbrains.com/plugin/14307-poetry) is available to help with using the Poetry virtual environments in PyCharm.

We use git pre-commit hooks to automatically run code formatting and linting on commit. These hooks are managed using a utility called `pre-commit`. After `poetry install` run the following command in this directory to install the hooks:

```commandline
pre-commit install --install-hooks
```

#### Install additional libraries if needed
If you need to install additional python libraries (e.g. boto3), you can run the command:
```commandline
poetry add boto3
```
This will install the latest version of `boto3` and all it's dependencies and update the `poetry.lock` file. You can also add a version constraint to this if you don't want the latest version, e.g.
```commandline
poetry add "boto3>=1.14,<1.16"
```

#### Upgrading libraries
To update all libraries to their latest versions, while still matching the constraints in `${COVID_PIPELINE_ROOTDIR}/pyproject.toml`, run:
```commandline
poetry update
```
If you want to update a package to it's latest version, run
```commandline
poetry add package@latest
```
To update to a specific version that is not the latest version, re-run the add command specifying a different version constraint.


### Database
All database schema migrations are managed using `sqitch` tool.

#### Dependencies
Prerequisites for running the `sqitch`:

* `sqitch` is installed in the machine. See [sqitch downloads page](https://sqitch.org/download/)
* `Postgres` with `psql` installed in the environment. Has postgres user set up
* `Postgres` has password exported to env variable
```commandline
export PGPASSWORD=some_secret_password
```
* A dedicated database is created for the project. It may be created using the following cmd:
```commandline
createdb -h localhost -U postgres bahrain_sars_cov_2
```

#### Working with sqitch
The work is done in `${COVID_PIPELINE_ROOTDIR}/sqitch/` directory

To add the new migration, use the following command:
```commandline
sqitch add 0x-my_migration_file_name -n 'My changes described here'
```
An `.sql` file will be created in `deploy`, `revert` and `verify`. Populate the files.
Important - `verify` scripts only fail, if `.sql` query raises an exception. Create verify scripts
accordingly to throw exceptions in case of failed verification

To check the migration status (are we missing any migrations?):
```commandline
sqitch status db:pg:bahrain_sars_cov_2
```

Migrations can be made using the following:
```commandline
sqitch deploy db:pg:bahrain_sars_cov_2
```

To verify migrations, which were made:
```commandline
sqitch verify db:pg:bahrain_sars_cov_2
```

To revert the changes:
```commandline
sqitch revert db:pg:bahrain_sars_cov_2
```

## Troubleshootings

### Sqitch
If authentication fails, try adding connection info to the command-line. For example:
```commandline
sqitch --db-user postgres --db-host localhost --db-port 5432 deploy db:pg:bahrain_sars_cov_2
```

### Nextstrain
To manually test nextstrain:
1. Start a nextstrain container
2. copy the fasta and metadata to: `<nextstrain_container>:/nextstrain/data/nextstrain.fasta` and `<nextstrain_container>:/nextstrain/data/nextstrain_metadata.tsv`
  without changing the destination names.
3. Exec nextstrain container
4. Inside the contaienr, run `cd /nextstrain ; snakemake --profile /custom_profile`

Note: the custom configuration is in `<nextstrain_container>:/custom_profile`
