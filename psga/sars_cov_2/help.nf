def printPathogenConfig() {

    log.info"""
        SARS-CoV-2 pathogen
        ===================
        Environment variables:
        * SARS_COV_2_PIPELINE_DOCKER_IMAGE            : ${SARS_COV_2_PIPELINE_DOCKER_IMAGE}
        * SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG        : ${SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG}
        * NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG : ${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG}
        * NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG : ${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG}
        * PANGOLIN_DOCKER_IMAGE_TAG                   : ${PANGOLIN_DOCKER_IMAGE_TAG}
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
        SARS_COV_2_PIPELINE_DOCKER_IMAGE
                                The name of the sars-cov-2 docker image
        SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG
                                The tag of the sars-cov-2 docker image
        NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG
                                The tag of the ncov2019-artic-nf-illumina docker image
        NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG
                                The tag of the ncov2019-artic-nf-nanopore docker image
        PANGOLIN_DOCKER_IMAGE_TAG
                                The tag of the pangolin docker image
    """.stripIndent()
}

