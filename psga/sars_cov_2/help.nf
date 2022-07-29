def printPathogenConfig() {

    log.info"""
        SARS-CoV-2 pathogen
        ===================
        Environment variables:
        * SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG        : ${SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG}
        * NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG : ${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG}
        * NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG : ${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG}
        * PANGOLIN_DOCKER_IMAGE_TAG                   : ${PANGOLIN_DOCKER_IMAGE_TAG}

        Parameters:
        * sequencing_technology                 : ${params.sequencing_technology}
        * scheme_repo_url                       : ${params.scheme_repo_url}
        * scheme_dir                            : ${params.scheme_dir}
        * scheme                                : ${params.scheme}
        * scheme_version                        : ${params.scheme_version}
    """.stripIndent()
}

def printPathogenHelp() {
    log.info"""
    SARS-CoV-2 pathogen
      Description:
        Map sequencing reads to consensus sequences to phylogenetic lineages.
          - Nanopore: ARTIC (https://github.com/artic-network/fieldbioinformatics)
          - Illumina: iVar (https://github.com/andersen-lab/ivar)
          - Pangolin: pangolin (https://github.com/cov-lineages/pangolin)

      Usage:
        nextflow run . --run [analysis_run] --sequencing_technology [sequencing_technology] [workflow-options]

      Mandatory environment variables:
        SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG
                                The tag of the sars-cov-2 docker image
        NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG
                                The tag of the ncov2019-artic-nf-illumina docker image
        NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG
                                The tag of the ncov2019-artic-nf-nanopore docker image
        PANGOLIN_DOCKER_IMAGE_TAG
                                The tag of the pangolin docker image

      Mandatory parameters:
        --sequencing_technology The sequencing technology for to use for ncov2019_artic pipeline: 'illumina' (input file extension: .fastq.gz or .bam), 'ont' (input file extension: .fastq), 'unknown' (input file extension: .fasta). If unknown, ncov is not executed.

      Optional parameters:
        --scheme_version        ARTIC scheme version (Default: 'V3')
        --scheme_repo_url       Repo to download your primer scheme from (e.g. 'https://github.com/artic-network/artic-ncov2019'). For efficiency, this repo was checked out and made available to the pipeline in the ncov docker images.
        --scheme_dir            Directory within scheme_repo_url that contains primer schemes (Default: 'primer_schemes')
        --scheme                Scheme name (Default: 'nCoV-2019')
    """.stripIndent()
}

