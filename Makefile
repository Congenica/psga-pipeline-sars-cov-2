DOCKER_IMAGE_PREFIX=566277102435.dkr.ecr.eu-west-2.amazonaws.com/congenica/psga-dev
DOCKER_IMAGE_PREFIX_LEGACY=144563655722.dkr.ecr.eu-west-1.amazonaws.com/congenica/dev
PSGA_PIPELINE_DOCKER_IMAGE_TAG_BASE=1.0.6
NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG_BASE=1.0.6
NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG_BASE=1.0.6
PANGOLIN_DOCKER_IMAGE_TAG_BASE=1.0.3
SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG=1.0.0
NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG=1.0.0
NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG=1.0.0
PANGOLIN_DOCKER_IMAGE_TAG=1.0.0

# create base images
base-images:
	docker build -t ${DOCKER_IMAGE_PREFIX_LEGACY}/psga-pipeline-base:${PSGA_PIPELINE_DOCKER_IMAGE_TAG_BASE} -t ${DOCKER_IMAGE_PREFIX}/psga-pipeline-base:${PSGA_PIPELINE_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.psga-pipeline-base .
	docker build -t ${DOCKER_IMAGE_PREFIX_LEGACY}/ncov2019-artic-nf-illumina-base:${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG_BASE} -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-illumina-base:${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.ncov2019-artic-nf-illumina-base .
	docker build -t ${DOCKER_IMAGE_PREFIX_LEGACY}/ncov2019-artic-nf-nanopore-base:${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG_BASE} -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-nanopore-base:${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.ncov2019-artic-nf-nanopore-base .
	docker build -t ${DOCKER_IMAGE_PREFIX_LEGACY}/pangolin-base:${PANGOLIN_DOCKER_IMAGE_TAG_BASE} -t ${DOCKER_IMAGE_PREFIX}/pangolin-base:${PANGOLIN_DOCKER_IMAGE_TAG_BASE} -f docker/Dockerfile.pangolin-base .

# build images per pathogen
sars-cov-2-images:
	docker build --build-arg pathogen=sars_cov_2 -t ${DOCKER_IMAGE_PREFIX}/sars-cov-2-pipeline:${SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG} -f docker/Dockerfile.psga-pipeline .
	docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-illumina:${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-illumina .
	docker build -t ${DOCKER_IMAGE_PREFIX}/ncov2019-artic-nf-nanopore:${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG} -f docker/Dockerfile.ncov2019-artic-nf-nanopore .
	docker build -t ${DOCKER_IMAGE_PREFIX}/pangolin:${PANGOLIN_DOCKER_IMAGE_TAG} -f docker/Dockerfile.pangolin .
