# ps-bahrain-covid

## Poetry

`Poetry` manages Python dependencies. Dependencies are declared in `pyproject.toml` and exact versions of both dependencies and sub-dependencies are stored in `poetry.lock` both of which must be committed to the git repo.

### Setting up your local environment

Follow the [instructions](https://python-poetry.org/docs/) to install Poetry however suits your local environment. Once installed simply run

```commandline
poetry install
```

in the git checkout to create a virtual environment and install all dependencies. A [plugin](https://plugins.jetbrains.com/plugin/14307-poetry) is available to help with using the Poetry virtual environments in PyCharm.

We use git pre-commit hooks to automatically run code formatting and linting on commit. These hooks are managed using a utility called `pre-commit`. After poetry install run the following command in this directory to install the hooks:

```commandline
pre-commit install --install-hooks
```

### Adding libraries

To add a new library, run (for example):

```commandline
poetry add boto3
```

This will install the latest version of `boto3` and all it's dependencies and update the `poetry.lock` file. You can also add a version constraint to this if you don't want the latest version, e.g.

```commandline
poetry add "boto3>=1.14,<1.16"
```

### Upgrading libraries

Run

```commandline
poetry update
```

to update all libraries to their latest versions, while still matching the constraints in `pyproject.toml`.

If you want to update a package to it's latest version, run

```commandline
poetry add package@latest
```

To update to a specific version that is not the latest version, re-run the add command specifying a different version constraint.

## Sqitch

`Sqitch` manages all schema migrations. Prerequisites for running the `sqitch`:

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

### Working with sqitch

The work is done in `sqitch/` directory

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

Manage migrations, deployment via docker image. Commands following postgres: status, deploy, verify, revert:
```commandline
docker run -v ${HOME}/ps-bahrain-covid/sqitch:/sqitch -w /sqitch -e PGPASSWORD=very_secret_password sqitch/sqitch --db-host sap-11281-bahrain-covid-pipeline.cbutwa8dg5pz.eu-west-1.rds.amazonaws.com --db-user postgres status db:pg:bahrain_sars_cov_2
```


### Sqitch troubleshoot

If authentication fails, try adding connection info to the command-line. For example:
```commandline
sqitch --db-user postgres --db-host localhost --db-port 5432 deploy db:pg:bahrain_sars_cov_2
```


### COVID pipeline

##### Environment variables

Environment variables required to run the pipeline. These environment variables are set in covid-pipeline/nextflow.config

| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| COVID_PIPELINE_ROOTDIR | Directory path, where all the pipeline git code is stored |
| COVID_PIPELINE_FASTQ_PATH | Directory path, where the input fastq files are stored |
| COVID_PIPELINE_WORKDIR | Directory path, where all the pipeline output is stored |
| COVID_PIPELINE_FASTA_PATH | Directory path, where re-headered ncov fasta files will be stored |
| COVID_PIPELINE_FASTA_PATH_QC_FAILED | Directory path, where re-headered ncov QC_FAILED fasta files will be stored  |

##### Install dependencies

```commandline
# get this s3 folder - it contains only a couple of FASTQ files for now
aws s3 cp s3://congenica-development-data-share/Bahrain_COVID_s3_data_lite ~/Bahrain_COVID_s3_data_lite --recursive

# install nextflow (copy nextflow cmd to a dir in your PATH):
wget -qO- https://get.nextflow.io | bash

# build the covid-pipeline image
docker build -t covid-pipeline:1.0.0 .

# set up the projects: ncov, pangolin
git submodule init
git submodule update

# build ncov docker image
cd ncov2019-artic-nf; docker build -f environments/illumina/Dockerfile -t ncov2019_artic_nf_base:1.0.0 . ; cd ..
docker build --build-arg NCOV_BASE_IMAGE_TAG=1.0.0 -f Dockerfile.ncov -t ncov2019_artic_nf:1.0.0 .

# build pangolin docker image
docker build -f Dockerfile.pangolin -t pangolin:1.0.0 .
```

##### Run covid-pipeline

```commandline
cd covid-pipeline
nextflow run .
```

