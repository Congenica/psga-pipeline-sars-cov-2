# ps-bahrain-covid

## Poetry

`Poetry` manages Python dependencies. Dependencies are declared in `pyproject.toml` and exact versions of both dependencies and sub-dependencies are stored in `poetry.lock` both of which must be committed to the git repo.

### Setting up your local environment

Follow the [instructions](https://python-poetry.org/docs/) to install Poetry however suits your local environment. Once installed simply run

```commandline
poetry install
```

in the git checkout to create a virtual environment and install all dependencies. A [plugin](https://plugins.jetbrains.com/plugin/14307-poetry) is available to help with using the Poetry virtual environments in PyCharm.

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

### Sqitch troubleshoot

If authentication fails, try adding connection info to the command-line. For example:
```commandline
sqitch --db-user postgres --db-host localhost --db-port 5432 deploy db:pg:bahrain_sars_cov_2
```


### COVID pipeline

To run this covid pipeline, the ncov illumina docker image must be edited and built first:

```
# get this s3 folder - it contains only a couple of FASTQ files for now
aws s3 cp s3://congenica-development-data-share/SAP-18211_Bahrain_COVID ~/Bahrain_COVID_s3_data_lite --recursive

# install nextflow: https://www.nextflow.io/docs/latest/getstarted.html

# set up the ncov project
git submodule init
git submodule update

cd ncov2019-artic-nf
# this is not required for running the pipeline using dockers, but for some reason it is required for running
# the pipeline with docker as part of our workflow. I haven't yet figured out why it works in one case and not in the other..
echo "COPY bin/qc.py /opt/conda/envs/artic-ncov2019-illumina/bin" >> environments/illumina/Dockerfile

# build the ncov docker image
docker build -f environments/illumina/Dockerfile -t ncov2019_edited:latest .
cd -


mkdir ~/ncov_results


# run our covid pipeline which executes ncov
cd covid-pipeline
nextflow run .



# Not sure whether it is best calling ncov as a workflow directly or via nested nextflow pipeline..


```

