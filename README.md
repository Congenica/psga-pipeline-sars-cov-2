# Pathogen SequencinG Analysis (PSGA) pipeline
The PSGA project is developed for the UK Health and Security Agency (UKHSA).

Currently, the only supported pathogen is: SARS-CoV-2. For this pathogen, Congenica sequencing protocol is based on the ARTIC consortium and works as follows:
- produce viral genome assemblies from sequence data (ncov2019-artic-nf);
- assign epidemiological lineages (Pangolin-Scorpio)

Support for other pathogens will be added.
A diagram for this pipeline is shown below:
![Alt text](img/psga.png?raw=true "PSGA pipeline")

## Operation
This pipeline runs on a Kubernetes environment. For the time being, there are two k8s deployments:
* `psga`, which allows for the execution of the nextflow pipeline within a k8s pod
* `psql`, which is a postgresql database accessible by the pods within the environment
Ideally, the postgresql database should be stored in an RDS aurora system outside the cluster, but for development, the current setting is fine.

The following diagram offers an overview of the pipeline execution in k8s. Each nextflow process is executed on a dedicated pod spun up by the main `psga` pod.
![Alt text](img/psga_k8s.png?raw=true "PSGA pipeline in k8s environment")


### Environment variables and input parameters
See: psga/modules/help.nf .

The help can also be printed with the command: `nextflow run . --help`.
The current configuration can be printed with the command: `nextflow run . --print_config`.


### Input files stored in aws s3
To process samples stored in s3, you need to export the env var `PSGA_INPUT_PATH` (see section: `Environment variables and input parameters`) to point to the s3 dir containing your data.
Input paths containing test datasets can be found below:

Small size datasets (processing time: few minutes):
- s3://synthetic-data-dev/UKHSA/small_tests/illumina_fastq          (2 samples)
- s3://synthetic-data-dev/UKHSA/small_tests/illumina_bams           (2 samples)
- s3://synthetic-data-dev/UKHSA/small_tests/medaka_fastq            (2 samples)

Medium size datasets (processing time: 1-2 hours). These datasets are used in our Jenkins CI validation (see jenkins/ dir):
- s3://synthetic-data-dev/UKHSA/validation_ci/illumina_artic_fastq    (30 samples: 10 alpha, 10 delta, 10 omicron)
- s3://synthetic-data-dev/UKHSA/validation_ci/illumina_artic_bam      (30 samples: 10 alpha, 10 delta, 10 omicron)
- s3://synthetic-data-dev/UKHSA/validation_ci/ont_artic_fastq         (30 samples: mix variants)

Large size datasets (processing time: 4-6 hours). These datasets were used as part of the validation in the UKHSA tender:
- s3://synthetic-data-dev/UKHSA/validation/illumina_ARTIC_fastq_COG_MARCH     (100 samples, alpha variants)
- s3://synthetic-data-dev/UKHSA/validation/illumina_ARTIC_fastq_COG_OCT2021   (100 samples, mostly delta variants)
- s3://synthetic-data-dev/UKHSA/validation/illumina_ARTIC_fastq_COG_FEB2022   (100 samples, omicron variants)
- s3://synthetic-data-dev/UKHSA/validation/illumina_ARTIC_bam_COG_MAR2021   (100 samples, alpha variants)
- s3://synthetic-data-dev/UKHSA/validation/illumina_ARTIC_bam_COG_OCT2021   (100 samples, mostly delta variants)
- s3://synthetic-data-dev/UKHSA/validation/illumina_ARTIC_bam_COG_FEB2022   (100 samples, omicron variants)
- s3://synthetic-data-dev/UKHSA/validation/ONT_ARTIC_fastq_COG   (300 samples, mix variants)


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
export PSGA_DOCKER_IMAGE_TAG_BASE=1.0.0
export NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG_BASE=1.0.0
export NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG_BASE=1.0.0
export PANGOLIN_DOCKER_IMAGE_TAG_BASE=1.0.0
export PSGA_DOCKER_IMAGE_TAG=1.0.0
export NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG=1.0.0
export NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG=1.0.0
export PANGOLIN_DOCKER_IMAGE_TAG=1.0.0

# create base images
docker build -t ${DOCKER_IMAGE_PREFIX}/psga-base:${PSGA_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.psga-base .
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-illumina-base:${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.ncov2019-artic-nf-illumina-base .
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-nanopore-base:${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.ncov2019-artic-nf-nanopore-base .
docker build -t ${DOCKER_IMAGE_PREFIX}/pangolin-base:${PANGOLIN_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.pangolin-base .


# build main images
docker build -t ${DOCKER_IMAGE_PREFIX}/psga:${PSGA_DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga .
docker build -t ${DOCKER_IMAGE_PREFIX}/psga-db:${PSGA_DOCKER_IMAGE_TAG} -f docker/Dockerfile.postgres .

# add project submodules
git submodule init
git submodule update

# update pangolin, ncov2019_artic_nf to their latest commits
# If you run this command, you need to regenerate the base images for pangolin and ncov2019
git submodule update --remote --merge

# build ncov docker images
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-illumina:${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-nanopore:${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .

# build pangolin docker image
docker build -t ${DOCKER_IMAGE_PREFIX}/pangolin:${PANGOLIN_DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .
```

Once all the required images are generated, the deployments can be created:
```commandline
cd minikube
./startup.sh

# input files can either be copied to the psga pod: /data/input
# or fetched from S3 (see examples: https://jira.congenica.net/browse/PSG-183).
# The env var: PSGA_INPUT_PATH must be set, accordingly.

# exec the psga pod
kubectl exec -it psga-XXXX -- bash

# ------------------
# WITHIN THE POD
# run the pipeline within the pod (processes are spun up as pod workers by this pipeline)
# the results will be stored in psga pod: /data/output
# MODE 1: Fresh run, overriding the output from the previous computations
nextflow run . <input_parameters>

# MODE 2: run from the last successful process
nextflow run . <input_parameters> -resume

# The following command cleans up the previous run's work directories and cache, but retains the content of ${PSGA_OUTPUT_PATH}:
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
To redeploy PSGA pipeline components:
```
cd k8s
./startup.sh
```
The DB is stored in RDS Aurora and can be accessed via the psga pod.



### GenBank submission

To submit files to GenBank, appropriate files and values need to be prepared for submission:
* Submission template. Navigate to https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ and fill out the
  form. A file .sbt will be available for download. Add this file path to nextflow parameter
  `genbank_submission_template` in configuration file  `psga/nextflow.config`
* Add Center/account abbreviation provided during account creation in MyNCBI to parameter
  `genbank_submitter_account_namespace` and `genbank_submitter_name`
  in configuration file  `psga/nextflow.config`
* Provide GenBank FTP connection details in `psga/nextflow.config`:
  * Username `genbank_storage_remote_username`
  * Password `genbank_storage_remote_password`
* Set the upload directory to `Production` for `genbank_storage_remote_directory` in `psga/nextflow.config`

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
For a description of these environment variables, see section `Environment variables and input parameters`.

```commandline
export DB_HOST=127.0.0.1
export DB_PORT=5432
export DB_NAME=psga_db
export DB_USER=postgres
export DB_PASSWORD=postgres

export PSGA_ROOT_PATH="/app"
export PSGA_INPUT_PATH="/data/input"
export PSGA_OUTPUT_PATH="/data/output"
```

#### Set up a local postgres database

A local database must be available to run the tests
```commandline
export DOCKER_IMAGE_PREFIX=144563655722.dkr.ecr.eu-west-1.amazonaws.com/congenica/dev
export PSGA_DOCKER_IMAGE_TAG=1.0.0

docker build -t ${DOCKER_IMAGE_PREFIX}/psga-db:${PSGA_DOCKER_IMAGE_TAG} -f docker/Dockerfile.postgres .

docker run -d -p ${DB_PORT}:${DB_PORT} --name my-postgres-server -e POSTGRES_PASSWORD=${DB_PASSWORD} ${DOCKER_IMAGE_PREFIX}/psga-db:${PSGA_DOCKER_IMAGE_TAG}

# test the connection from your local machine
psql -h ${DB_HOST} -p ${DB_PORT} -U ${DB_USER} -W

# once finished testing
docker stop my-postgres-server
docker rm my-postgres-server
```

All database schema migrations are managed using `sqitch` tool. Prerequisites for running the `sqitch`. The docker image psga-db
already has sqitch installed. The following commands are executed within the psga-db docker container:

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
cd ${PSGA_ROOT_PATH}/tests
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
To update all libraries to their latest versions, while still matching the constraints in `${PSGA_ROOT_PATH}/pyproject.toml`, run:
```commandline
poetry update
```
If you want to update a package to it's latest version, run
```commandline
poetry add package@latest
```
To update to a specific version that is not the latest version, re-run the add command specifying a different version constraint.


#### Working with sqitch
The work is done in `${PSGA_ROOT_PATH}/sqitch/` directory

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
