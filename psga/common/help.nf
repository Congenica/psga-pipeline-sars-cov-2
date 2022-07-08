def printMainConfig() {
    log.info"""
        =====================
        ${workflow.manifest.name} v ${workflow.manifest.version}
        =====================
        Global environment variables:
        * PSGA_ROOT_PATH                              : ${PSGA_ROOT_PATH}
        * PSGA_OUTPUT_PATH                            : ${PSGA_OUTPUT_PATH}
        * PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH          : ${PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH}
        * PSGA_MAX_ATTEMPTS                           : ${PSGA_MAX_ATTEMPTS}
        * PSGA_SLEEP_TIME_BETWEEN_ATTEMPTS            : ${PSGA_SLEEP_TIME_BETWEEN_ATTEMPTS}
        * DOCKER_IMAGE_PREFIX                         : ${DOCKER_IMAGE_PREFIX}
        * PSGA_PIPELINE_DOCKER_IMAGE_TAG              : ${PSGA_PIPELINE_DOCKER_IMAGE_TAG}
        * K8S_NODE                                    : ${K8S_NODE}
        * K8S_PULL_POLICY                             : ${K8S_PULL_POLICY}
        * K8S_SERVICE_ACCOUNT                         : ${K8S_SERVICE_ACCOUNT}
        * K8S_QUEUE_SIZE                              : ${K8S_QUEUE_SIZE}
        * K8S_STORAGE_CLAIM_NAME                      : ${K8S_STORAGE_CLAIM_NAME}
        * K8S_STORAGE_MOUNT_PATH                      : ${K8S_STORAGE_MOUNT_PATH}
        * NXF_WORK                                    : ${NXF_WORK}
        * NXF_EXECUTOR                                : ${NXF_EXECUTOR}
        * NXF_ANSI_LOG                                : ${NXF_ANSI_LOG}
        * NXF_OPTS                                    : ${NXF_OPTS}

        Global parameters:
        * metadata                                    : ${params.metadata}
        * run                                         : ${params.run}
    """.stripIndent()
}

def printMainHelp() {
    log.info"""
    Pathogen Sequence Genome Analysis pipeline.

    Generic configuration for all pathogens:
      Mandatory environment variables:
        PSGA_ROOT_PATH
                                Path to the pipeline code (e.g. git checkout). (e.g. /app) |
        PSGA_OUTPUT_PATH
                                Path to the whole pipeline output. (e.g. /data/output, s3://data/output)
        PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH
                                Path containing the analysis runs which are in progress or incomplete. (e.g. /data/incomplete_analysis_run)
        PSGA_MAX_ATTEMPTS
                                The maximum number of attempts to resume the an interrupted pipeline run
        PSGA_SLEEP_TIME_BETWEEN_ATTEMPTS
                                The sleep time between attempts in seconds
        DOCKER_IMAGE_PREFIX     The prefix of the docker image, excluded the image name
        PSGA_PIPELINE_DOCKER_IMAGE_TAG
                                The tag of the psga docker image
        K8S_NODE                The Kubernetes node for nodeAffinity
        K8S_PULL_POLICY         The Kubernetes docker image pull policy (e.g. Always, Never)
        K8S_SERVICE_ACCOUNT     The Kubernetes service account
        K8S_QUEUE_SIZE          The maximum number of processes to run at the same time (default: 20)
        K8S_STORAGE_CLAIM_NAME  The Kubernetes PVC claim
        K8S_STORAGE_MOUNT_PATH  The Kubernetes mount path (default: /data)
        K8S_PROCESS_MAX_RETRIES The maximum number that a process can be retried if a resource-based exit code (137-143) is raised (default: 3)
        K8S_PROCESS_CPU_LOW     Value for a process using little CPU. There is no need to change this as the pipeline was designed for high scalability (default: 1)
        K8S_PROCESS_CPU_HIGH    Value for a process using a lot of CPU. There is no need to change this as the pipeline was designed for high scalability (default: 2)
        K8S_PROCESS_MEMORY_VERY_LOW
                                Value for a process using very low memory in MB (default: 250)
        K8S_PROCESS_MEMORY_LOW  Value for a process using low memory in MB (default: 500)
        K8S_PROCESS_MEMORY_MEDIUM
                                Value for a process using medium memory in MB (default: 1500)
        K8S_PROCESS_MEMORY_HIGH Value for a process using high memory in MB (default: 3000)
        K8S_PROCESS_MEMORY_VERY_HIGH
                                Value for a process using very high memory in MB (default: 6000)
        NXF_WORK                Set Nextflow work directory (e.g. /data/work)
        NXF_EXECUTOR            Set Nextflow executor (default: k8s)
        NXF_ANSI_LOG            Enable Nextflow ANSI log (default: false)
        NXF_OPTS                Pass JVM options to Nextflow (default: -Xms1g -Xmx4g)

      Mandatory parameters:
        --metadata              The path to the metadata file. This can be an s3 path.
        --run                   A (unique) string identifying the analysis run (batch).

      Optional parameters:
        --help                  Print this help
        --print_config          Print the pipeline configuration
    """.stripIndent()
}
