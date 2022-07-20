DOCKER_IMAGE_URI_PATH=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-dev
DOCKER_IMAGE_TAG=1.0.0


# build images per pathogen
sars-cov-2-images:
	docker build --build-arg pathogen=sars_cov_2 -t ${DOCKER_IMAGE_URI_PATH}/sars-cov-2-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .
	docker build -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-illumina:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
	docker build -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-nanopore:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .
	docker build -t ${DOCKER_IMAGE_URI_PATH}/pangolin:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .

dummy-pathogen-images:
	docker build --build-arg pathogen=dummy_pathogen -t ${DOCKER_IMAGE_URI_PATH}/dummy-pathogen-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .
