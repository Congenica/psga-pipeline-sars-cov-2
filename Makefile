DOCKER_IMAGE_URI_PATH=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-nonprod
DOCKER_IMAGE_TAG=dev_latest
SARS_COV_2=sars_cov_2
DOCKER_IMAGE_NAME=sars-cov-2-pipeline

# build images per pathogen
sars-cov-2-images:
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/sars-cov-2-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/htstools:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.htstools .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/read-it-and-keep:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.read-it-and-keep .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/fastqc:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.fastqc .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/primer-autodetection:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.primer-autodetection .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-illumina:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-nanopore:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/pangolin:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .

build-sars-cov-2-local:
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .


shell_local: build-sars-cov-2-local
	docker run \
	--rm \
	-it \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash

test:
	poetry run pytest tests/
