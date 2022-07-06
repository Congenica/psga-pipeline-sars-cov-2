# Pathogen Sequence Genome Analysis (PSGA) pipeline
The architecture of the PSGA pipeline is described [here](https://confluence.congenica.net/display/PSG/PSGA+Pipeline+Architecture).

Currently, the only supported pathogen is: SARS-CoV-2. For this pathogen, Congenica sequencing protocol is based on the ARTIC consortium and works as follows:
- produce viral genome assemblies from sequence data (ncov2019-artic-nf);
- assign epidemiological lineages (Pangolin-Scorpio)

Support for other pathogens will be added.

## Operation
This pipeline runs on a Kubernetes environment. The main workflow coordinates the execution of processes within Kubernetes jobs.


### Environment variables and input parameters
See: psga/modules/help.nf .

The help can also be printed with the command: `nextflow run . --help`.
The current configuration can be printed with the command: `nextflow run . --print_config`.


### Input files stored in aws s3
To process samples stored in s3, set up a metadata CSV file (see tests/test_data/good_metadata.csv for reference) including the paths to the sample input files. Two files are required for running illumina fastq samples. 1 file is required for running illumina bam / nanopore medaka fastq / no_ncov fasta samples.

For each pathogen, test datasets can be found in `jenkins/files`.


### Running the pipeline using K8s Minikube (local testing)
Install git submodules
```commandline
# initialise the submodules and fetch the code from origin of this repository
make submodule-setup
# update the submodules from origin of their repositories
make submodule-update
```

Download and install Minikube using the instructions provided here: https://minikube.sigs.k8s.io/docs/start/ .
Once minikube is active, we need to activate the minikube registry of docker images so that Minikube can find these locally.
See: https://medium.com/swlh/how-to-run-locally-built-docker-images-in-kubernetes-b28fbc32cc1d for additional ideas
```commandline
eval $(minikube -p minikube docker-env)
```

Build the pipeline docker images in the minikube docker environment:
```commandline
make build-base-images
make build-images
```

Once all the required images are generated, the deployments can be created:
```commandline
cd minikube
./startup.sh

# exec the psga pod
kubectl exec -it psga-pipeline-XXXX -- bash

# ------------------
# WITHIN THE psga POD
# run the pipeline within the pod (processes are spun up as pod workers by this pipeline). The results will be stored in psga-pipeline pod: /data/output
# use `-resume` flag to resume the previous pipeline execution
nextflow run . -c <pathogen>.config <input_parameters>

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
All the docker images mentioned in the Minikube section must be built and pushed to Congenica ECR.
To redeploy PSGA pipeline components:
```
cd k8s
./startup.sh
```


## Development

### Adding new pathogens
In order to add the pathogen `pathogenX` to the pipeline, change dir to `psga` and follow the instructions below:
1. run the script: `python initialise_pathogen.py --pathogen-name pathogenX`
2. add nextflow configs and workflows to the following files: `pathogenX.config`, `pathogenX/psga.nf`, `pathogenX/help.nf`.
3. add Python scripts to: `../scripts/pathogenX/`
4. edit the DB schema as necessary
5. run the pipeline from the psga directory using the command: `nextflow run . -c pathogenX.config <pathogenX-specific parameters>`

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


#### Set up the required environment variables
For a description of these environment variables, see section `Environment variables and input parameters`.


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
