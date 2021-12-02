# Covid-pipeline documentation

## Operation

### Environment variables

Environment variables required to run the pipeline.
These environment variables must be available in the system.

| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| DB_HOST | Postgres database host address (e.g. 192.168.0.33) |
| DB_PORT | Postgres database port (e.g. 5432) |
| DB_NAME | Postgres database name (e.g. covid_pipeline_db) |
| DB_USER | Postgres database user name (e.g. postgres) |
| DB_PASSWORD | Postgres database user password (e.g. postgres) |
| COVID_PIPELINE_ROOTDIR | Path to the pipeline code (e.g. git checkout). Default: ${HOME}/covid-pipeline |
| COVID_PIPELINE_FASTQ_PATH | Path to the input FASTQ files and TSV metadata file. Default: ${HOME}/COVID_s3_data_lite/sample_data |
| COVID_PIPELINE_WORKDIR | Path to the whole pipeline output. Default: ${HOME}/covid-pipeline-workdir |
| COVID_PIPELINE_REPORTS_PATH | Path to the pipeline reports. Default: ${HOME}/covid-pipeline-reports |


The following environment variables are set internally and should not be changed
| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| COVID_PIPELINE_MISSING_METADATA_PATH | Path to the missing metadata files. Set to: ${COVID_PIPELINE_WORKDIR}/no-metadata-found-fastq |
| COVID_PIPELINE_NCOV_OUTPUT_PATH | Path to store all ncov2019-artic result files. Each run will be published to unique folder |
| COVID_PIPELINE_QC_PLOTS_PATH | Path to store all ncov2019-artic qc_plots graphs in single folder |
| COVID_PIPELINE_FASTA_PATH | Path to the re-headered ncov FASTA files. Set to: ${COVID_PIPELINE_WORKDIR}/reheadered-fasta |
| COVID_PIPELINE_FASTA_PATH_QC_FAILED | Path to the re-headered ncov QC_FAILED FASTA files. Set to: ${COVID_PIPELINE_WORKDIR}/reheadered-fasta-qc-failed |
| COVID_PIPELINE_PANGOLIN_PATH | Path to the results of pangolin pipeline with lineage reports. ach run will be published to unique folder |
| COVID_PIPELINE_GENBANK_PATH | Path to submission files, which were used to submit samples to GenBank programmatic interface |
| COVID_PIPELINE_NEXTSTRAIN_PATH | Path to store all nextstrain result files. Each run will be published to unique folder |
| COVID_PIPELINE_MICROREACT_PATH | Path to store microreact tsv file, generated from samples, found in database |
| COVID_PIPELINE_NOTIFICATIONS_PATH | Path to the pipeline notifications. Unexpected events regarding missing samples, files are reported here in text files |


### Requirements

#### Set up the required environment variables (the following configuration is just an example)
```commandline
export DB_HOST=127.0.0.1
export DB_PORT=5432
export DB_NAME=covid_pipeline_db
export DB_USER=postgres
export DB_PASSWORD=postgres

export COVID_PIPELINE_ROOTDIR="${HOME}/covid-pipeline"
export COVID_PIPELINE_FASTQ_PATH="${HOME}/COVID_s3_data_lite/sample_data"
export COVID_PIPELINE_WORKDIR="${HOME}/covid-pipeline-workdir"
export COVID_PIPELINE_REPORTS_PATH="${HOME}/covid-pipeline-reports"
```

#### Set up a local postgres database

A local database can be set up using a postgres docker image:
```commandline
docker pull postgres:11.14
docker run -d -p ${DB_PORT}:${DB_PORT} --name my-postgres-server -e POSTGRES_PASSWORD=${DB_PASSWORD} postgres:11.14

# test the connection from your local machine
psql -h ${DB_HOST} -p ${DB_PORT} -U ${DB_USER} -W

# once finished testing
docker stop my-postgres-server
docker rm my-postgres-server
```

All database schema migrations are managed using `sqitch` tool. Prerequisites for running the `sqitch`:

* `sqitch` is installed in the machine. See [sqitch downloads page](https://sqitch.org/download/). Choose docker installation.
* `Postgres` with `psql` installed in the environment. Has postgres user set up
* `Postgres` has password exported to env variable

Create a dedicated database for the project:
```commandline
export PGPASSWORD=${DB_PASSWORD}
createdb -h ${DB_HOST} -U ${DB_USER} ${DB_NAME}
```

Finally, deploy the required DB migrations:
```commandline
# all migrations are in the project sqitch dir
cd ${COVID_PIPELINE_ROOTDIR}/sqitch/
SQITCH_URI=db:pg://${DB_USER}:${DB_PASSWORD}@${DB_HOST}:${DB_PORT}/${DB_NAME}
sqitch deploy ${SQITCH_URI}
```

The following sqitch commands can be helpful for diagnosis:
```commandline
# check the migration status (are we missing any migrations?):
sqitch status ${SQITCH_URI}

# verify migrations, which were made:
sqitch verify ${SQITCH_URI}

# revert the changes:
sqitch revert ${SQITCH_URI}
```


#### Build the docker containers
```commandline
export VERSION=1.0.0

# install nextflow (copy nextflow executable to a directory in your PATH environment variable):
wget -qO- https://get.nextflow.io | bash

# build the covid-pipeline image
docker build -t covid-pipeline:${VERSION} .

# add project submodules
git submodule init
git submodule update

# update pangolin, ncov2019_artic_nf, and nextstrain to their latest commits
git submodule update --remote --merge

# build ncov docker image
docker build -t ncov2019_artic_nf:${VERSION} -f Dockerfile.ncov2019-artic-nf .

# build pangolin docker image
docker build -t pangolin:${VERSION} -f Dockerfile.pangolin .

# build nextstrain docker image
docker build -t nextstrain:${VERSION} -f Dockerfile.nextstrain .

# build auspice docker image (for visualising the results with Auspice web-service):
docker build -t auspice:${VERSION} -f Dockerfile.auspice .
```

### GenBank submission

To submit files to GenBank, appropriate files and values need to be prepared for submission:
* Submission template. Navigate to https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ and fill out the
  form. A file .sbt will be available for download. Add this file path to nextflow parameter
  `genbank_submission_template` in configuration file  `covid-pipeline/nextflow.config`
* Add Center/account abbreviation provided during account creation in MyNCBI to parameter
  `genbank_submitter_account_namespace` and `genbank_submitter_name`
  in configuration file  `covid-pipeline/nextflow.config`
* Provide GenBank FTP connection details in `covid-pipeline/nextflow.config`:
  * Username `genbank_storage_remote_username`
  * Password `genbank_storage_remote_password`
* Set the upload directory to `Production` for `genbank_storage_remote_directory` in `covid-pipeline/nextflow.config`

Samples to GenBank are submitted once. When Sample is uploaded to GenBank, it is marked with session id. After that,
sample won't be submitted to GenBank

If any of credentials are missing for GenBank, submission is skipped

#### GenBank Submission schema

File `scripts/genbank/submission.py` was generated using submission schema
[submission.xsd](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/common/).
A tool [generateDS](https://pypi.org/project/generateDS/) was used to generate a file.
A command line used to generate a file:
```commandline
generateDS -o "scripts/genbank/genbank_submission.py" submission.xsd
```

### Run covid-pipeline

The pipeline can be executed in two modes:
```commandline
cd ${COVID_PIPELINE_ROOTDIR}/covid-pipeline

# MODE 1: Retain the previous computations (e.g. samples processed previously)
nextflow run . -resume

# MODE 2: Overriding the output from the previous computations
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
  -v ${COVID_PIPELINE_WORKDIR}/nextstrain/latest/nextstrain_output/bahrain/ncov_with_accessions.json:/ncov_with_accessions.json \
  -p 4000:4000 \
  auspice:${VERSION} \
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

#### Testing the pipeline
Unit tests are implemented with pytest and stored in the project tests dir
```commandline
cd ${COVID_PIPELINE_ROOTDIR}/tests
pytest
```

To run the full pipeline:
```commandline
cd ${COVID_PIPELINE_ROOTDIR}/covid-pipeline
nextflow run .
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

#### Working with sqitch
The work is done in `${COVID_PIPELINE_ROOTDIR}/sqitch/` directory

To add the new migration, use the following command:
```commandline
sqitch add 0x-my_migration_file_name -n 'My changes described here'
```
An `.sql` file will be created in `deploy`, `revert` and `verify`. Populate the files.
Important - `verify` scripts only fail, if `.sql` query raises an exception. Create verify scripts
accordingly to throw exceptions in case of failed verification


## Troubleshootings

### Sqitch
If authentication fails, try adding connection info to the command-line. For example:
```commandline
sqitch --db-user ${DB_USER} --db-host ${DB_HOST} --db-port ${DB_PORT} deploy db:pg:${DB_NAME}
```

### Nextstrain
To manually test nextstrain:
1. Start a nextstrain container
2. copy the fasta and metadata to: `<nextstrain_container>:/nextstrain/data/nextstrain.fasta` and `<nextstrain_container>:/nextstrain/data/nextstrain_metadata.tsv`
  without changing the destination names.
3. Exec nextstrain container
4. Inside the contaienr, run `cd /nextstrain ; snakemake --profile /custom_profile`

Note: the custom configuration is in `<nextstrain_container>:/custom_profile`
