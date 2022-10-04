DOCKER_IMAGE_URI_PATH=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-dev
SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG=1.0.0
NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG=1.0.0
NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG=1.0.0
PANGOLIN_DOCKER_IMAGE_TAG=1.0.0


# build images per pathogen
sars-cov-2-images:
	docker build --build-arg pathogen=sars_cov_2 -t ${DOCKER_IMAGE_URI_PATH}/sars-cov-2-pipeline:${SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .
	docker build -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-illumina:${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
	docker build -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-nanopore:${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .
	docker build -t ${DOCKER_IMAGE_URI_PATH}/pangolin:${PANGOLIN_DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .
