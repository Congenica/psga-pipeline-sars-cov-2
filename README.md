# Pathogen Sequence Genome Analysis (PSGA) pipeline

This repository contains a collection of pipelines for pathogen sequence genome analysis.

## Supported pathogens

### SARS-CoV-2

Congenica sequencing protocol is based on the ARTIC consortium and works as follows:

- produce viral genome assemblies from sequence data (ncov2019-artic-nf);
- assign epidemiological lineages (Pangolin-Scorpio)

## Operation

This pipeline runs on a Kubernetes environment. The main workflow coordinates the execution of processes within Kubernetes jobs.

### Environment variables and input parameters

See: app/modules/help.nf .

The help can also be printed with the command: `nextflow run . --help`.
The current configuration can be printed with the command: `nextflow run . --print_config`.

### Pipeline metadata.csv

To process samples, set up a metadata CSV file (e.g. tests/test_data/check_metadata/good_metadata_illumina.csv) which stores the sample id and the sample files. Two files are required for running illumina fastq samples. 1 file is required for running illumina bam / nanopore medaka fastq / fasta samples.

### Running the pipeline using K8s Minikube (local testing)

Download and install Minikube using the instructions provided here: https://minikube.sigs.k8s.io/docs/start/ .
Once minikube is active, we need to activate the minikube registry of docker images so that Minikube can find these locally.
See: https://medium.com/swlh/how-to-run-locally-built-docker-images-in-kubernetes-b28fbc32cc1d for additional ideas

```commandline
eval $(minikube -p minikube docker-env)
```

Note on Apple Silicon ensure you DO NOT have DOCKER_DEFAULT_PLATFORM set.

If you get errors around kubelet not being available try the following:

```
minikube delete --all --purge
unset DOCKER_DEFAULT_PLATFORM

# Mounting is required for reference data to work
# Recommended to throw more than the default resources at it
minikube start --mount-string="/$PWD/ref-data:/app/resources" --mount --cpus=5 --memory=7000
```

Docker Desktop 4.27 broke minikube which breaks docker running amd64 images. Set the base image to this in the meantime.

```
minikube start --base-image gcr.io/k8s-minikube/kicbase-builds:v0.0.42-1703092832-17830
```

Make sure that the submodules are installed:

```commandline
git submodule update --init --recursive
```

Build the pipeline docker images in the minikube docker environment:

```commandline
make sars-cov-2-images
```

Once all the required images are generated, the deployments can be created:

```commandline
cd minikube
./startup.sh

# exec the <pathogen> pod
kubectl exec -it <pathogen>-pipeline-XXXX -- bash
kubectl exec -it sars-cov-2-pipeline-minikube -- bash

# ------------------
# WITHIN THE pathogen POD
# run the pipeline within the pod (processes are spun up as pod workers by this pipeline). The results will be stored in: /data/output
# use `-resume` flag to resume the previous pipeline execution
nextflow run . --metadata <metadata_path> --run <analysis_run> --sequencing_technology <sequencing_technology> --kit <kit> --output_path <output_path>
nextflow run . --metadata <metadata_path> --run <analysis_run> --sequencing_technology <sequencing_technology> --kit <kit> --output_path <output_path>

Running sars-cov-2/fasta test
nextflow -log unknown.log run /app/main.nf --run unknown --sequencing_technology unknown --kit none --config-path s3://psga-test-static-data/pipeline_ci_tests/sars_cov_2/unknown/none/ --output_path /data/output/pipeline_ci_unknown
Running sars-cov-2/illumina_unknown test
nextflow -log illumina.log run /app/main.nf --run illumina --sequencing_technology illumina --kit unknown --config-path s3://psga-test-static-data/pipeline_ci_tests/sars_cov_2/illumina/unknown/ --output_path /data/output/pipeline_ci_illumina
Running sars-cov-2/ont_unknown test
nextflow -log ont.log run /app/main.nf --run ont --sequencing_technology ont --kit unknown --config-path s3://psga-test-static-data/pipeline_ci_tests/sars_cov_2/ont/unknown/ --output_path /data/output/pipeline_ci_ont

# The following command cleans up the previous run's work directories and cache, but retains the published output:
nextflow clean -f

# once finished
exit
# ------------------

# when finished:
./cleanup.sh
```

N.B.
This requires AWS Access and Secret keys set in the env, you you can procure via:

(Fill in your own profile name)

```
CREDS=`aws configure export-credentials --profile <profile-name>`
echo export AWS_SECRET_ACCESS_KEY=$(echo ${CREDS} | jq .SecretAccessKey | tr -d \")
echo export AWS_ACCESS_KEY_ID=$(echo ${CREDS} | jq .AccessKeyId | tr -d \")
echo export AWS_SESSION_TOKEN=\"$(echo ${CREDS} | jq .SessionToken | tr -d \")\"
```

## Development

### Install dependency packages using Python Poetry tool

`Poetry` manages Python dependencies. Dependencies are declared in `pyproject.toml` and exact versions of both dependencies and sub-dependencies are stored in `poetry.lock`. Both are committed to the git repo.

#### Setting up your local environment

The pipeline scripts and tests are code for Python 3.10.
Follow the [instructions](https://python-poetry.org/docs/) to install Poetry depending on your local environment. Then run:

```commandline
poetry install
```

in the git checkout to create a virtual environment and install all dependencies. A [plugin](https://plugins.jetbrains.com/plugin/14307-poetry) is available to help with using the Poetry virtual environments in PyCharm.

We use git pre-commit hooks to automatically run code formatting and linting on commit. These hooks are managed using a utility called `pre-commit`. After `poetry install` run the following command in this directory to install the hooks:

```commandline
pre-commit install --install-hooks
```

#### Set up the required environment variables

For a description of these environment variables, see section `Environment variables and input parameters`.

#### Run the unit tests

Unit tests are implemented with pytest and stored in the project tests dir

```commandline
cd tests
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

## TODO

- extract scripts to use template
- install + use results writer
- stop using tee in read it and keep to get errors (test with two files of different sizes)
- Integrate ncov pipelines
- Stop running generate files script
- Remove python from pyproject or environment files
- bump mypy commit hook
- add **init** to at least scripts.util
- Stop copying integration tests into the docker image

### Notes on installing packages

- Pangolin v4.3 requires exactly v1.74 of BioPython.
  Github has a release of 4.3.1 which is more generous, however this hasn't been released to conda.
  When this gets updated, it should be possible to bring the pangolin code into the main docker image.
  An alternative is to install the pangolin release directly from github and then install dependencies in the docker image.
  ```
  pangolin = {git = "https://github.com/cov-lineages/pangolin.git", rev = "v4.3.1"}
  ```

## notes

- This repo has tests and scripts outside the app dir
- Scripts are copied in at runtime. The tests that run on them assume they do not live inside app/.
- Integration tests are copied into the docker image.

## UPDATED - Running and testing locally

Platform is now set in all the docker images but to ensure you are building for the correct platform.

```
export DOCKER_DEFAULT_PLATFORM=linux/amd64
```

### Fastq processing

Sars-CoV-2 can be built using

```
make build_sars_cov_2_local
```

There is a command to enter a shell with app mounted to /app.

```
make mounted_shell_local
```

Once inside the shell the updated pipeline can be run using

```
seqt="ont"
nextflow run ./fastq.nf --run 61c06b0a-e5e8-4dbf-8bb0-729cce46a223 -params-file /app/local_test/$seqt/settings.json --config-path /app/local_test/$seqt/ --output_path /app/output/$seqt/
```

Note, set seqt to "illumina" to run illumina.

### NCOV

```
make build_ncov_local

cd app
seqt="ncov_ont"
nextflow run ./workflows/ncov.nf --run 61c06b0a-e5e8-4dbf-8bb0-729cce46a223 -params-file /app/local_test/$seqt/settings.json --config-path /app/local_test/$seqt/ --output_path /app/output/$seqt/
```
