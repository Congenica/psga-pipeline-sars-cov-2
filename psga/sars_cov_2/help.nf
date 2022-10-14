def printPathogenConfig() {

    log.info"""
        SARS-CoV-2 pathogen
        ===================
        Environment variables:
        * SARS_COV_2_PIPELINE_DOCKER_IMAGE            : ${SARS_COV_2_PIPELINE_DOCKER_IMAGE}
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
    """.stripIndent()
}

