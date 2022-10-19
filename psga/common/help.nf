def printMainConfig() {
    log.info"""
        =====================
        ${workflow.manifest.name} v ${workflow.manifest.version}
        =====================
        Global environment variables:
        * DOCKER_IMAGE_URI_PATH                       : ${DOCKER_IMAGE_URI_PATH}
        * DOCKER_IMAGE_TAG                            : ${DOCKER_IMAGE_TAG}
        * AWS_ACCESS_KEY                              : ${AWS_ACCESS_KEY}
        * AWS_SECRET_ACCESS_KEY                       : ${AWS_SECRET_ACCESS_KEY}
        * AWS_REGION                                  : ${AWS_REGION}
        * AWS_ROLE_ARN                                : ${AWS_ROLE_ARN}
        * AWS_CONNECTION_TIMEOUT                      : ${AWS_CONNECTION_TIMEOUT}
        * AWS_MAX_CONNECTIONS                         : ${AWS_MAX_CONNECTIONS}
        * AWS_MAX_PARALLEL_TRANSFERS                  : ${AWS_MAX_PARALLEL_TRANSFERS}
        * K8S_NODE                                    : ${K8S_NODE}
        * K8S_PULL_POLICY                             : ${K8S_PULL_POLICY}
        * K8S_SERVICE_ACCOUNT                         : ${K8S_SERVICE_ACCOUNT}
        * K8S_STORAGE_CLAIM_NAME                      : ${K8S_STORAGE_CLAIM_NAME}
        * K8S_STORAGE_MOUNT_PATH                      : ${K8S_STORAGE_MOUNT_PATH}
        * QUEUE                                       : ${QUEUE}
        * QUEUE_SIZE                                  : ${QUEUE_SIZE}
        * NXF_WORK                                    : ${NXF_WORK}
        * NXF_EXECUTOR                                : ${NXF_EXECUTOR}
        * NXF_OPTS                                    : ${NXF_OPTS}

        Global parameters:
        * metadata                                    : ${params.metadata}
        * run                                         : ${params.run}
        * sequencing_technology                       : ${params.sequencing_technology}
        * kit                                         : ${params.kit}
        * output_path                                 : ${params.output_path}

    """.stripIndent()
}

def printMainHelp() {
    log.info"""
    Pathogen Sequence Genome Analysis pipeline.

    Generic configuration for all pathogens:
      Mandatory environment variables:
        DOCKER_IMAGE_URI_PATH   The prefix of the docker image, excluded the image name
        DOCKER_IMAGE_TAG
                                The tag of the docker images docker image
        AWS_ACCESS_KEY
                                The AWS access key
        AWS_SECRET_ACCESS_KEY
                                The AWS secret access key
        AWS_REGION
                                The AWS region
        AWS_ROLE_ARN
                                The AWS role ARN
        AWS_CONNECTION_TIMEOUT
                                The amount of time to wait (in milliseconds) when initially establishing a connection before giving up and timing out
        AWS_MAX_CONNECTIONS
                                The maximum number of allowed open HTTP connections
        AWS_MAX_PARALLEL_TRANSFERS
                                Max parallel upload/download transfer operations per job
        K8S_NODE                The Kubernetes node for nodeAffinity
        K8S_PULL_POLICY         The Kubernetes docker image pull policy (e.g. Always, Never)
        K8S_SERVICE_ACCOUNT     The Kubernetes service account
        K8S_STORAGE_CLAIM_NAME  The Kubernetes PVC claim
        K8S_STORAGE_MOUNT_PATH  The Kubernetes mount path (default: /data)
        QUEUE                   The name of the queue where jobs are scheduled when using a grid based executor in your pipeline
        QUEUE_SIZE              The maximum number of processes to run at the same time (default: 20)
        PROCESS_MAX_RETRIES     The maximum number that a process can be retried if a non-zero exit status is returned (default: 3)
        PROCESS_CPU_LOW         Value for a process using little CPU. There is no need to change this as the pipeline was designed for high scalability (default: 1)
        PROCESS_CPU_HIGH        Value for a process using a lot of CPU. There is no need to change this as the pipeline was designed for high scalability (default: 2)
        PROCESS_MEMORY_VERY_LOW
                                Value for a process using very low memory in MB (default: 250)
        PROCESS_MEMORY_LOW  Value for a process using low memory in MB (default: 500)
        PROCESS_MEMORY_MEDIUM
                                Value for a process using medium memory in MB (default: 1500)
        PROCESS_MEMORY_HIGH Value for a process using high memory in MB (default: 3000)
        PROCESS_MEMORY_VERY_HIGH
                                Value for a process using very high memory in MB (default: 6000)
        NXF_WORK                Set Nextflow work directory (e.g. /data/work; s3://work/dir)
        NXF_EXECUTOR            Set Nextflow executor (e.g. k8s, awsbatch)
        NXF_OPTS                Pass JVM options to Nextflow (default: -Xms1g -Xmx4g)

      Mandatory parameters:
        --metadata              The path to the metadata file. This can be an s3 path
        --run                   A (unique) string identifying the analysis run (batch)
        --sequencing_technology The technology used for sequencing the samples. Values: 'illumina', 'ont', 'unknown'
        --kit                   The kit used for sequencing the samples (e.g. the scheme version of the primers)
        --output_path
                                Path to the whole pipeline output. (e.g. /data/output, s3://data/output)

      Optional parameters:
        --help                  Print this help
        --print_config          Print the pipeline configuration
    """.stripIndent()
}
