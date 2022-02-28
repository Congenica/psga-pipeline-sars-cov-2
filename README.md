# Covid-pipeline documentation

## Operation

This pipeline runs on a Kubernetes environment. For the time being, there are two k8s deployments:
* `covid-pipeline`, which allows for the execution of the nextflow pipeline within a k8s pod
* `psql`, which is a postgresql database accessible by the pods within the environment
Ideally, the postgresql database should be stored in an RDS aurora system outside the cluster, but for development, the current setting is fine.

The following diagram offers an overview of the pipeline execution in k8s. Each nextflow process is executed on a dedicated pod spun up by the main `covid-pipeline` pod.
![Alt text](img/UKHSA_covid_project.png?raw=true "Covid pipeline in k8s environment")


### Environment variables

Environment variables required to run the pipeline. They are set up in the covid-pipeline k8s deployment.

| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| DB_HOST | Postgres database host address (e.g. 192.168.0.33) |
| DB_PORT | Postgres database port (e.g. 5432) |
| DB_NAME | Postgres database name (e.g. covid_pipeline_db) |
| DB_USER | Postgres database user name (e.g. postgres) |
| DB_PASSWORD | Postgres database user password (e.g. postgres) |
| COVID_PIPELINE_ROOT_PATH | Path to the pipeline code (e.g. git checkout). (e.g. /app) |
| COVID_PIPELINE_INPUT_PATH | Path to the required input BAM/FASTQ files and TSV metadata file. (e.g. /data/input, s3://synthetic-data-dev/UKHSA/piero-test-data/illumina_fastq ) |
| COVID_PIPELINE_OUTPUT_PATH | Path to the whole pipeline output. (e.g. /data/output) |


The following environment variables are set internally and should not be changed
| Variable | Description |
| :---------------- | :---------------------------------------------------------------- |
| COVID_PIPELINE_MISSING_METADATA_PATH | Path to the missing metadata files. Set to: ${COVID_PIPELINE_OUTPUT_PATH}/no-metadata-found-bam |
| COVID_PIPELINE_NCOV_OUTPUT_PATH | Path to store all ncov2019-artic result files. Each run will be published to unique folder |
| COVID_PIPELINE_QC_PLOTS_PATH | Path to store all ncov2019-artic qc_plots graphs in single folder |
| COVID_PIPELINE_FASTA_PATH | Path to the re-headered ncov FASTA files. Set to: ${COVID_PIPELINE_OUTPUT_PATH}/reheadered-fasta |
| COVID_PIPELINE_FASTA_PATH_QC_FAILED | Path to the re-headered ncov QC_FAILED FASTA files. Set to: ${COVID_PIPELINE_OUTPUT_PATH}/reheadered-fasta-qc-failed |
| COVID_PIPELINE_PANGOLIN_PATH | Path to the results of pangolin pipeline with lineage reports. Each run will be published to unique folder |
| COVID_PIPELINE_GENBANK_PATH | Path to submission files, which were used to submit samples to GenBank programmatic interface |
| COVID_PIPELINE_NOTIFICATIONS_PATH | Path to the pipeline notifications. Unexpected events regarding missing samples, files are reported here in text files |
| K8S_PULL_POLICY | The Kubernetes docker image pull policy (e.g. Always, Never) |

### Pipeline input parameters

Input parameters to run the pipeline.

| Argument | Value |
| :---------------- | :---------------------------------------------------------------- |
| --workflow | illumina_artic (default; input file extension: .fq.gz or .bam), medaka_artic (nanopore workflow; input file extension: .fastq.gz). |
| --filetype | fastq (default), bam . The type of input file. Currently bam is only supported by the illumina workflow |
| --run | The name for this analysis run |
| --scheme_repo_url | Repo to download your primer scheme from (e.g. 'https://github.com/artic-network/artic-ncov2019'). For efficiency, this repo was checked out and made available to the pipeline in the ncov docker images |
| --scheme_dir | Directory within schemeRepoURL that contains primer schemes (Default: 'primer_schemes') |
| --scheme | Scheme name (Default: 'nCoV-2019') |
| --scheme_version | ARTIC scheme version (Default: 'V3') |


Example of execution with parameter: `nextflow run . --workflow medaka_artic`

### Input files stored in aws s3
If you plan to read input files from an aws s3 bucket you will need to:

- copy your aws credentials to the /root dir in the covid-pipeline pod
- export the env var `COVID_PIPELINE_INPUT_PATH` (see section: `Environment variables`) to point to the s3 dir containing your data. For quick tests we have these three paths:
s3://synthetic-data-dev/UKHSA/piero-test-data/illumina_fastq
s3://synthetic-data-dev/UKHSA/piero-test-data/illumina_bams
s3://synthetic-data-dev/UKHSA/piero-test-data/medaka_fastq_fail/20200311_1427_X1_FAK72834_a3787181
s3://synthetic-data-dev/UKHSA/piero-test-data/medaka_fastq_pass                                     (TO BE REVIEWED)

If you intend to use your own path, do not forget to store a metadata.tsv file as well.


### Running the pipeline using K8s Minikube (local testing)

Download and install Minikube using the instructions provided here: https://minikube.sigs.k8s.io/docs/start/

Once minikube is active, we need to activate the minikube registry of docker images so that Minikube can find these locally.
See: https://medium.com/swlh/how-to-run-locally-built-docker-images-in-kubernetes-b28fbc32cc1d for additional ideas
```commandline
eval $(minikube -p minikube docker-env)
```

The next step is to build the pipeline docker images in the minikube docker environment. For simplicity, the database is stored on a pod. This is not ideal as this can be lost if the pod crashes or is deleted. However, as a proof of concept, this is fine. In the future, the database will be stored in an RDS aurora cluster, therefore outside the k8s environment.
```commandline
export DOCKER_IMAGE_PREFIX=144563655722.dkr.ecr.eu-west-1.amazonaws.com/congenica/dev
export VERSION_BASE=1.0.0
export VERSION=1.0.0

# create base images
docker build -t ${DOCKER_IMAGE_PREFIX}/covid-pipeline-base:${VERSION_BASE} -f docker/Dockerfile.covid-pipeline-base .
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-illumina-base:${VERSION_BASE} -f docker/Dockerfile.ncov2019-artic-nf-illumina-base .
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-nanopore-base:${VERSION_BASE} -f docker/Dockerfile.ncov2019-artic-nf-nanopore-base .
docker build -t ${DOCKER_IMAGE_PREFIX}/pangolin-base:${VERSION_BASE} -f docker/Dockerfile.pangolin-base .


# build main images
docker build -t ${DOCKER_IMAGE_PREFIX}/covid-pipeline:${VERSION} -f docker/Dockerfile.covid-pipeline .
docker build -t ${DOCKER_IMAGE_PREFIX}/covid-pipeline-db:${VERSION} -f docker/Dockerfile.postgres .

# add project submodules
git submodule init
git submodule update

# update pangolin, ncov2019_artic_nf to their latest commits
git submodule update --remote --merge

# build ncov docker images
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-illumina:${VERSION} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-nanopore:${VERSION} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .

# build pangolin docker image
docker build -t ${DOCKER_IMAGE_PREFIX}/pangolin:${VERSION} -f docker/Dockerfile.pangolin .
```

Once all the required images are generated, the deployments can be created:
```commandline
cd minikube
./startup.sh

# input files can either be copied to the covid-pipeline pod: /data/input
# or fetched from S3 (see examples: https://jira.congenica.net/browse/PSG-183).
# The env var: COVID_PIPELINE_INPUT_PATH must be set, accordingly.

# exec the covid-pipeline pod
kubectl exec -it covid-pipeline-XXXX -- bash

# ------------------
# WITHIN THE POD
# run the pipeline within the pod (processes are spun up as pod workers by this pipeline)
# the results will be stored in covid-pipeline pod: /data/output
# MODE 1: Fresh run, overriding the output from the previous computations
nextflow run . <input_parameters>

# MODE 2: run from the last successful process
nextflow run . <input_parameters> -resume

# The following command cleans up the previous run's work directories and cache, but retains the content of ${COVID_PIPELINE_OUTPUT_PATH}:
nextflow clean -f

# once finished
exit
# ------------------

# when finished:
./cleanup.sh
```


### Running the pipeline using K8s (currently in Congenica saas-dev cluster)
Start a new shell to make sure that the standard docker environment is used and not the one dedicated to minikube.
All the docker images mentioned in the Minikube section, except for `docker/Dockerfile.postgres` must be built and pushed to Congenica ECR.
To redeploy the covid pipeline components:
```
cd k8s
./startup.sh
```
The DB is stored in RDS Aurora and can be accessed via the covid-pipeline pod.



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


#### Set up the required environment variables (the following configuration is just an example)
```commandline
export DB_HOST=127.0.0.1
export DB_PORT=5432
export DB_NAME=covid_pipeline_db
export DB_USER=postgres
export DB_PASSWORD=postgres

export COVID_PIPELINE_ROOT_PATH="${HOME}/covid-pipeline"
export COVID_PIPELINE_INPUT_PATH="${HOME}/COVID_s3_data_lite/sample_data_0"
export COVID_PIPELINE_OUTPUT_PATH="${HOME}/covid-pipeline-output"
```

#### Set up a local postgres database

A local database must be available to run the tests
```commandline
export DOCKER_IMAGE_PREFIX=144563655722.dkr.ecr.eu-west-1.amazonaws.com/congenica/dev
export VERSION=1.0.0

docker build -t ${DOCKER_IMAGE_PREFIX}/covid-pipeline-db:${VERSION} -f docker/Dockerfile.postgres .

docker run -d -p ${DB_PORT}:${DB_PORT} --name my-postgres-server -e POSTGRES_PASSWORD=${DB_PASSWORD} ${DOCKER_IMAGE_PREFIX}/covid-pipeline-db:${VERSION}

# test the connection from your local machine
psql -h ${DB_HOST} -p ${DB_PORT} -U ${DB_USER} -W

# once finished testing
docker stop my-postgres-server
docker rm my-postgres-server
```

All database schema migrations are managed using `sqitch` tool. Prerequisites for running the `sqitch`. The docker image covid-pipeline-db
already has sqitch installed. The following commands are executed within the covid-pipeline-db docker container:

Create a dedicated database for the project:
```commandline
export PGPASSWORD=${DB_PASSWORD}
createdb -h ${DB_HOST} -U ${DB_USER} ${DB_NAME}
```

Finally, deploy the required DB migrations:
```commandline
# all migrations are in the project sqitch dir
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

#### Run the unit tests
Unit tests are implemented with pytest and stored in the project tests dir
```commandline
cd ${COVID_PIPELINE_ROOT_PATH}/tests
pytest
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
To update all libraries to their latest versions, while still matching the constraints in `${COVID_PIPELINE_ROOT_PATH}/pyproject.toml`, run:
```commandline
poetry update
```
If you want to update a package to it's latest version, run
```commandline
poetry add package@latest
```
To update to a specific version that is not the latest version, re-run the add command specifying a different version constraint.


#### Working with sqitch
The work is done in `${COVID_PIPELINE_ROOT_PATH}/sqitch/` directory

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
