DOCKER_IMAGE_URI_PATH=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-nonprod
DOCKER_IMAGE_TAG=dev_latest
SARS_COV_2=sars_cov_2
DOCKER_IMAGE_NAME=sars-cov-2-pipeline

TEST_RUN_ID=10a0d649-045d-4419-a4c3-90892c0aa583
TEST_OUTPUT_LOCAL=/app/output/${TEST_RUN_ID}

CONTAINER_TEST_DATA_PATH=/app/local_test
ONT_TEST_DATA_PATH=${CONTAINER_TEST_DATA_PATH}/ont/
ILLUMINA_TEST_DATA_PATH=${CONTAINER_TEST_DATA_PATH}/illumina/

# build images per pathogen
sars-cov-2-images:
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/sars-cov-2-pipeline:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline-sars-cov-2 .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-illumina:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/ncov2019-artic-nf-nanopore:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/pangolin:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .

build_sars_cov_2_local:
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline-sars-cov-2 .


build_local_images: build_sars_cov_2_local
	docker build --build-arg pathogen=${SARS_COV_2} -t ncov2019-artic-nf-illumina:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
	# docker build --build-arg pathogen=${SARS_COV_2} -t ncov2019-artic-nf-nanopore:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .
	# docker build --build-arg pathogen=${SARS_COV_2} -t pangolin:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .

shell_local: build_sars_cov_2_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app:/app \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash


mounted_shell_local:# build_sars_cov_2_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:/app/local_test \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash


test_fastq_local_ont:# build_sars_cov_2_local
	docker run \
	--rm \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./fastq.nf \
		-params-file ${ONT_TEST_DATA_PATH}settings.json \
		--config-path ${ONT_TEST_DATA_PATH}
		--output_path ${TEST_OUTPUT_LOCAL}/ont

test_fastq_local_illumina:# build_sars_cov_2_local
	docker run \
	--rm \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./fastq.nf \
		-params-file ${ILLUMINA_TEST_DATA_PATH_TEST_DATA_PATH}settings.json \
		--config-path ${ILLUMINA_TEST_DATA_PATH}
		--output_path ${TEST_OUTPUT_LOCAL}/illumina

test_local: build_sars_cov_2_local
	docker run \
	--rm \
	-it \
	--env-file app/test/env.test \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run . \
		-params-file /app/test/settings.json \
		--console.echo \
		--process.executor=local \
		--run ${TEST_RUN_ID} \
		--config-path /app/test \
		--output_path /app/output/${TEST_RUN_ID}

test:
	poetry run pytest tests/
