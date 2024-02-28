DOCKER_IMAGE_URI_PATH=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-nonprod
DOCKER_IMAGE_TAG=dev_latest
SARS_COV_2=sars_cov_2
DOCKER_IMAGE_NAME=sars-cov-2-pipeline
NCOV_DOCKER_IMAGE_NAME=ncov2019-artic-nf
NCOV_ONT_DOCKER_IMAGE_NAME=${NCOV_DOCKER_IMAGE_NAME}-nanopore
NCOV_ILLUMINA_DOCKER_IMAGE_NAME=${NCOV_DOCKER_IMAGE_NAME}-illumina
PANGOLIN_DOCKER_IMAGE_NAME=pangolin

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

build_ncov_local:
	docker build --build-arg pathogen=${SARS_COV_2} -t ${NCOV_ONT_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.${NCOV_ONT_DOCKER_IMAGE_NAME} .
	docker build --build-arg pathogen=${SARS_COV_2} -t ${NCOV_ILLUMINA_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.${NCOV_ILLUMINA_DOCKER_IMAGE_NAME} .

build_pangolin_local:
	docker build --build-arg pathogen=${SARS_COV_2} -t ${PANGOLIN_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .

build_minikube_local:
	minikube image build --build-opt=build-arg=pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline-sars-cov-2 .
	# minikube image build --build-opt=build-arg=pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${NCOV_ONT_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.${NCOV_ONT_DOCKER_IMAGE_NAME} .
	# minikube image build --build-opt=build-arg=pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${NCOV_ILLUMINA_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.${NCOV_ILLUMINA_DOCKER_IMAGE_NAME} .
	minikube image build --build-opt=build-arg=pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${PANGOLIN_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .

	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${NCOV_ONT_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.${NCOV_ONT_DOCKER_IMAGE_NAME} .
	minikube image load ${DOCKER_IMAGE_URI_PATH}/${NCOV_ONT_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}
	docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${NCOV_ILLUMINA_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.${NCOV_ILLUMINA_DOCKER_IMAGE_NAME} .
	minikube image load ${DOCKER_IMAGE_URI_PATH}/${NCOV_ILLUMINA_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}

	# docker build --build-arg pathogen=${SARS_COV_2} -t ${DOCKER_IMAGE_URI_PATH}/${PANGOLIN_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .
	# minikube image load ${DOCKER_IMAGE_URI_PATH}/${PANGOLIN_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}

reload_minikube:
	kubectl delete --ignore-not-found=true -f minikube/pipelines/sars-cov-2.yaml
	sleep 35
	minikube image rm ${DOCKER_IMAGE_URI_PATH}/${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}
	make build_minikube_local
	kubectl apply -f minikube/pipelines/sars-cov-2.yaml

shell_local: build_sars_cov_2_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app:/app \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash


mounted_shell_local: build_sars_cov_2_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash


# Running this requires commenting/disabling nextflow.config
test_fastq_local_ont: build_sars_cov_2_local
	docker run \
	--rm \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./main.nf \
		-c /app/nextflow.local.config \
		--run 61c06b0a-e5e8-4dbf-8bb0-729cce46a224 \
		-params-file ${ONT_TEST_DATA_PATH}settings.json \
		--config-path ${ONT_TEST_DATA_PATH}
		--output_path ${TEST_OUTPUT_LOCAL}/ont

# Running this requires commenting/disabling nextflow.config
test_fastq_local_illumina: build_sars_cov_2_local
	docker run \
	--rm \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./main.nf \
		-c /app/nextflow.local.config \
		--run 61c06b0a-e5e8-4dbf-8bb0-729cce46a224 \
		-params-file ${ILLUMINA_TEST_DATA_PATH}settings.json \
		--config-path ${ILLUMINA_TEST_DATA_PATH}
		--output_path ${TEST_OUTPUT_LOCAL}/illumina


mounted_ncov_ont_shell_local: build_ncov_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app/modules:/app/modules \
	--volume ${PWD}/app/workflows:/app/workflows \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${NCOV_ONT_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash


test_ncov_ont_local: build_ncov_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app/modules:/modules \
	--volume ${PWD}/app/workflows:/workflows \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	--volume ${PWD}/app/output/ncov_ont/:/app/output/ncov_ont/ \
	${NCOV_ONT_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./workflows/ncov.nf \
		-c /app/local_test/ncov_ont/nextflow.config \
		--run 61c06b0a-e5e8-4dbf-8bb0-729cce46a224 \
		-params-file /app/local_test/ncov_ont/settings.json \
		--config-path /app/local_test/ncov_ont/ \
		--output_path /app/output/ncov_ont/

mounted_ncov_illumina_shell_local: build_ncov_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app/modules:/app/modules \
	--volume ${PWD}/app/workflows:/app/workflows \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${NCOV_ILLUMINA_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash

test_ncov_illumina_local: build_ncov_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app/modules:/modules \
	--volume ${PWD}/app/workflows:/workflows \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	--volume ${PWD}/app/output/ncov_illumina/:/app/output/ncov_illumina/ \
	${NCOV_ILLUMINA_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./workflows/ncov.nf \
		-c /app/local_test/ncov_illumina/nextflow.config \
		--run 61c06b0a-e5e8-4dbf-8bb0-729cce46a223 \
		-params-file /app/local_test/ncov_illumina/settings.json \
		--config-path /app/local_test/ncov_illumina/ \
		--output_path /app/output/ncov_illumina/


mounted_pangolin_shell_local: build_pangolin_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app/modules:/app/modules \
	--volume ${PWD}/app/workflows:/app/workflows \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${PANGOLIN_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	bash

# N.B. for this to work, you need to add
#   - nextflow=23.10.1
# To docker/sars_cov_2/pangolin.yml
test_pangolin_local: build_pangolin_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app/modules:/modules \
	--volume ${PWD}/app/workflows:/workflows \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	--volume ${PWD}/app/output/pangolin/:/app/output/pangolin/ \
	${PANGOLIN_DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./workflows/pangolin.nf \
		--run 61c06b0a-e5e8-4dbf-8bb0-729cce46a223 \
		--config-path /app/local_test/fasta/ \
		--output_path /app/output/pangolin/

# N.B. for this to work, you need to add
#   - nextflow=23.10.1
# To docker/sars_cov_2/pangolin.yml
test_fasta_local: build_sars_cov_2_local
	docker run \
	--rm \
	-it \
	--volume ${PWD}/app:/app \
	--volume ${PWD}/local_test/:${CONTAINER_TEST_DATA_PATH} \
	${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG} \
	nextflow \
		run ./main.nf \
		-c /app/nextflow.local.config \
		--run 61c06b0a-e5e8-4dbf-8bb0-729cce46a223 \
		--config-path /app/local_test/fasta/ \
		-params-file /app/local_test/fasta/settings.json \
		--output_path /app/output/fasta/

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
