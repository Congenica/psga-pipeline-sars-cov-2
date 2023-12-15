DOCKER_IMAGE_URI_PATH=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-nonprod
DOCKER_IMAGE_TAG=dev_latest
SARS_COV_2=sars_cov_2
SYNTHETIC=synthetic
S_AUREUS=s_aureus

# build base images
bactopia-base-image:
	docker build -t ${DOCKER_IMAGE_URI_PATH}/bactopia:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.bactopia .

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

synthetic-images:
	docker build --build-arg pathogen=${SYNTHETIC} -t ${DOCKER_IMAGE_URI_PATH}/synthetic-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .

s-aureus-images:
	docker build --build-arg pathogen=${S_AUREUS} -t ${DOCKER_IMAGE_URI_PATH}/s-aureus-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .
	docker build --build-arg pathogen=${S_AUREUS} -t ${DOCKER_IMAGE_URI_PATH}/s-aureus:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.s-aureus .


build-sars-cov-2-local:
		docker build --build-arg pathogen=${SARS_COV_2} -t sars-cov-2-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .

test:
	poetry run pytest tests/
