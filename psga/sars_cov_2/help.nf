def printPathogenConfig() {

    log.info"""
        SARS-CoV-2 pathogen
        ===================
        Environment variables:
        * NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG : ${NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG}
        * NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG : ${NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG}
        * PANGOLIN_DOCKER_IMAGE_TAG                   : ${PANGOLIN_DOCKER_IMAGE_TAG}

        Parameters:
        * ncov_workflow                         : ${params.ncov_workflow}
        * filetype                              : ${params.filetype}
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
        nextflow run . -c sars_cov_2.config --run [analysis_run] --ncov_workflow [ncov_workflow] --filetype [filetype] [workflow-options]

      Mandatory environment variables:
        NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG
                                The tag of the ncov2019-artic-nf-illumina docker image
        NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG
                                The tag of the ncov2019-artic-nf-nanopore docker image
        PANGOLIN_DOCKER_IMAGE_TAG
                                The tag of the pangolin docker image

      Mandatory parameters:
        --ncov_workflow         The ncov2019artic workflow to run: 'illumina_artic' (input file extension: .fastq.gz or .bam), 'medaka_artic' (input file extension: .fastq), 'no_ncov' (input file extension: .fasta).
        --filetype              The type of input file: 'fasta', 'fastq' or 'bam'.

      Optional parameters:
        --scheme_version        ARTIC scheme version (Default: 'V3')
        --scheme_repo_url       Repo to download your primer scheme from (e.g. 'https://github.com/artic-network/artic-ncov2019'). For efficiency, this repo was checked out and made available to the pipeline in the ncov docker images.
        --scheme_dir            Directory within scheme_repo_url that contains primer schemes (Default: 'primer_schemes')
        --scheme                Scheme name (Default: 'nCoV-2019')
    """.stripIndent()
}

