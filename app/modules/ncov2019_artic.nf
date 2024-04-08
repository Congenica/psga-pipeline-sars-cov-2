/*
 * Run: ncov2019-artic-nf nextflow pipeline (illumina workflow)
 * see: https://github.com/Congenica/ncov2019-artic-nf
 *
 * Run: ncov2019-artic-nf nextflow pipeline (nanopore/medaka workflow)
 * see: https://github.com/Congenica/ncov2019-artic-nf
 * Note: This runs as a shell block
 */
process NCOV2019_ARTIC_NF_PIPELINE {
  publishDir "${params.output_path}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: 'output_{bam,fasta,plots,typing,variants}/*'
  publishDir "${params.output_path}", mode: 'copy', overwrite: true, pattern: 'reheadered_fasta/*.fasta'

  input:
    tuple val(meta), path(read_paths)

  output:
    tuple val(meta), path("reheadered_fasta/*.fasta"), emit: ch_ncov_sample_fasta
    path "ncov_output/*.qc.csv", emit: ch_ncov_qc_csv
    path "output_bam/*.bam*"
    path "output_fasta/*.fa"
    path "output_plots/*.png"
    path "output_typing/*"
    path "output_variants/*.vcf.gz*", optional: { params.sequencing_technology != "ont"} // ont
    path "output_variants/*.tsv", optional: { params.sequencing_technology != "illumina"} // illumina

  shell:

    sample_id = meta.SAMPLE_ID
    //   # extract scheme and version from primer
    (scheme, scheme_version) = meta.PRIMER.tokenize( '_' )

    ncov2019_artic_nf_typing_gff = "/MN908947.3.gff"
    ncov2019_artic_nf_typing_yaml = "/SARS-CoV-2.types.yaml"

    ncov_out_dir="ncov_output"
    untrimmed_bam_file_ext="sorted.bam"
    output_bam="output_bam"
    output_fasta="output_fasta"
    output_plots="output_plots"
    output_typing="output_typing"
    output_variants="output_variants"
    reheadered_fasta="reheadered_fasta"

    if( params.sequencing_technology == 'illumina' ) {
      '''
      # set nextflow env vars
      export NXF_ANSI_LOG="false"

      ncov_untrimmed_bam_out_dir="!{ncov_out_dir}/ncovIllumina_sequenceAnalysis_readMapping"
      ncov_trimmed_bam_out_dir="!{ncov_out_dir}/ncovIllumina_sequenceAnalysis_trimPrimerSequences"
      ncov_fasta_out_dir="!{ncov_out_dir}/ncovIllumina_sequenceAnalysis_makeConsensus"
      ncov_qc_plots_dir="!{ncov_out_dir}/qc_plots"
      ncov_tsv_out_dir="!{ncov_out_dir}/ncovIllumina_sequenceAnalysis_callVariants"
      ncov_typing_out_dir="!{ncov_out_dir}/ncovIllumina_Genotyping_typeVariants"

      trimmed_bam_file_ext="mapped.primertrimmed.sorted.bam"

      mkdir -p !{sample_id}
      for fastq_file in "!{read_paths}"; do
        echo $fastq_file
        mv -f $fastq_file !{sample_id}
      done

      # note: we inject our configuration into ncov to override parameters
      # note: `pwd` is the workdir for this nextflow process
      # note: use /tmp as work dir, so that the intermediate files for this pipeline remain local
      #       instead of being shared with our pipeline

      nextflow run /ncov2019-artic-nf \
          --illumina \
          --prefix !{params.run} \
          --directory !{sample_id} \
          --outdir !{ncov_out_dir} \
          --schemeDir /primer_schemes \
          --scheme !{scheme} \
          --schemeVersion !{scheme_version} \
          --gff !{ncov2019_artic_nf_typing_gff} \
          --yaml !{ncov2019_artic_nf_typing_yaml} \
          -work-dir /tmp \
          -c /ncov-illumina.config

      mkdir -p !{output_bam} !{output_fasta} !{output_plots} !{output_typing} !{output_variants} !{reheadered_fasta}
      mv ${ncov_untrimmed_bam_out_dir}/*!{sample_id}.!{untrimmed_bam_file_ext} !{output_bam}
      mv ${ncov_trimmed_bam_out_dir}/*!{sample_id}.${trimmed_bam_file_ext} !{output_bam}

      mv ${ncov_fasta_out_dir}/*.primertrimmed.consensus.fa !{output_fasta}
      mv ${ncov_qc_plots_dir}/*.png !{output_plots}
      mv ${ncov_typing_out_dir}/**/* !{output_typing}
      mv ${ncov_tsv_out_dir}/*.variants.tsv !{output_variants}

      # rename the qc csv
      mv !{ncov_out_dir}/!{params.run}.qc.csv !{ncov_out_dir}/!{sample_id}.qc.csv

      # index bam files
      samtools index !{output_bam}/!{sample_id}.!{untrimmed_bam_file_ext} !{output_bam}/!{sample_id}.!{untrimmed_bam_file_ext}.bai
      samtools index !{output_bam}/!{sample_id}.${trimmed_bam_file_ext} !{output_bam}/!{sample_id}.${trimmed_bam_file_ext}.bai

      # reheader fasta file
      python ${PSGA_ROOT_PATH}/scripts/ncov/reheader_fasta.py --input-dir !{output_fasta} --output-dir !{reheadered_fasta}
      '''
    } else if( params.sequencing_technology == 'ont' ) {
        // We will need to gunzip the file
        gz_fastq_file = read_paths
        fastq_file = gz_fastq_file.getBaseName()

        // ncov-artic pipeline parameters
        // ONT/medaka parameters used by
        // artic minion
        // {pore}_{device}_{caller variant}_{caller version}
        default_medaka_model = 'r941_min_hac_variant_g507'
        normalise = 200
        // //guppyplex
        min_len = null
        max_len = null

        '''
        # set nextflow env vars
        export NXF_ANSI_LOG="false"

        ncov_out_dir="ncov_output"
        ncov_minion_out_dir="!{ncov_out_dir}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"
        ncov_qc_plots_dir="!{ncov_out_dir}/qc_plots"
        ncov_typing_out_dir="!{ncov_out_dir}/articNcovNanopore_Genotyping_typeVariants"

        trimmed_bam_file_ext="primertrimmed.rg.sorted.bam"
        vcf_file_ext="pass.vcf.gz"

        # ncov/ONT expects decompressed fastq in artic 1.1.3
        # This is always going to be compressed because contamination removal gzips it
        gunzip -c !{gz_fastq_file} > !{fastq_file}

        # We try to infer the medaka model, and fall back if it fails
        medaka_model=$(medaka tools resolve_model --auto_model consensus !{fastq_file} || echo !{default_medaka_model})

        # move fastq file to a specific directory so that the output files will have the filename pattern:
        # <analysis_run>_<sample_id>
        # where analysis_run is ncov_prefix, and
        # sample_id is the sample UUID (=file name of the fastq file)

        # there is only one input fastq file here as ncov is executed per sample
        mkdir -p !{sample_id}
        mv -f !{fastq_file} !{sample_id}

        # note: we inject our configuration into ncov to override parameters
        # note: --basecalled_fastq is the directory containing the barcodes or the fastq files
        # note: use /tmp as work dir, so that the intermediate files for this pipeline remain local
        #       instead of being shared with our pipeline
        nextflow run /ncov2019-artic-nf \
            --medaka \
            --medakaModel ${medaka_model} \
            --normalise !{normalise} \
            --min_len !{min_len} \
            --max_len !{max_len} \
            --prefix !{params.run} \
            --basecalled_fastq !{sample_id} \
            --outdir !{ncov_out_dir} \
            --schemeDir /primer_schemes \
            --scheme !{scheme} \
            --schemeVersion !{scheme_version} \
            --gff !{ncov2019_artic_nf_typing_gff} \
            --yaml !{ncov2019_artic_nf_typing_yaml} \
            -work-dir /tmp \
            -c /ncov-nanopore.config

        mkdir -p !{output_bam} !{output_fasta} !{output_plots} !{output_typing} !{output_variants} !{reheadered_fasta}
        mv ${ncov_minion_out_dir}/*!{sample_id}.!{untrimmed_bam_file_ext} !{output_bam}
        mv ${ncov_minion_out_dir}/*!{sample_id}.${trimmed_bam_file_ext} !{output_bam}
        mv ${ncov_minion_out_dir}/*.fasta !{output_fasta}
        mv ${ncov_qc_plots_dir}/*.png !{output_plots}
        mv ${ncov_typing_out_dir}/**/* !{output_typing}
        mv ${ncov_minion_out_dir}/*.${vcf_file_ext} !{output_variants}

        # this is a code correction to the nanopore medaka workflow in order to restore the correct sample names
        # in file content and file name
        # UPDATE sample name in file content
        for file_to_update in `ls !{output_bam}/!{params.run}_!{sample_id}*.bam`; do
            samtools view -H ${file_to_update}  | sed "s/!{params.run}_!{sample_id}/!{sample_id}/g" | samtools reheader - ${file_to_update} > ${file_to_update}.new
            mv ${file_to_update}.new ${file_to_update}
            samtools index ${file_to_update} ${file_to_update}.bai
        done

        for file_to_update in `ls !{output_fasta}/!{params.run}_!{sample_id}*.fasta`; do
            sed -i "s/!{params.run}_!{sample_id}/!{sample_id}/" ${file_to_update}
            # use .fa extension so that this is the same as for the ncov illumina workflow
            mv ${file_to_update} ${file_to_update%.*}.fa
        done

        for file_to_update in `ls !{output_typing}/!{params.run}_!{sample_id}*.*`; do
            sed -i "s/!{params.run}_!{sample_id}/!{sample_id}/" ${file_to_update}
        done

        for file_to_update in `ls !{output_variants}/!{params.run}_!{sample_id}*.vcf*`; do
            zcat ${file_to_update} | sed "s/!{params.run}_!{sample_id}/!{sample_id}/" | bgzip > "${file_to_update}.new"
            mv ${file_to_update}.new ${file_to_update}
            bcftools index -t "${file_to_update}"
        done

        sed -i "s/!{params.run}_!{sample_id}/!{sample_id}/g" !{ncov_out_dir}/!{params.run}.qc.csv
        mv !{ncov_out_dir}/!{params.run}.qc.csv !{ncov_out_dir}/!{sample_id}.qc.csv

        # UPDATE sample name in file name
        for file_to_rename in `find !{output_bam} !{output_fasta} !{output_plots} !{output_typing} !{output_variants} -name !{params.run}_!{sample_id}*`; do
            file_dir=`dirname ${file_to_rename}`
            file_name=`basename ${file_to_rename}`
            cd ${file_dir}
            rename "s/!{params.run}_!{sample_id}/!{sample_id}/" ${file_name}
            cd - > /dev/null
        done

        # TODO: FASTA Channel also does this
        # reheader fasta file
        python ${PSGA_ROOT_PATH}/scripts/ncov/reheader_fasta.py --input-dir !{output_fasta} --output-dir !{reheadered_fasta}
        '''
    }
}

// Processes run by ONT NCOV
// ```
// articNcovNanopore:sequenceAnalysisMedaka:articGuppyPlex (61c06b0a-e5e8-4dbf-8bb0-729cce46a223-4a32dffd-584c-4f2a-8b17-e1e27b630f66)
// articNcovNanopore:sequenceAnalysisMedaka:articMinIONMedaka (61c06b0a-e5e8-4dbf-8bb0-729cce46a223_4a32dffd-584c-4f2a-8b17-e1e27b630f66)
// articNcovNanopore:sequenceAnalysisMedaka:articRemoveUnmappedReads (61c06b0a-e5e8-4dbf-8bb0-729cce46a223_4a32dffd-584c-4f2a-8b17-e1e27b630f66)
// articNcovNanopore:sequenceAnalysisMedaka:makeQCCSV (61c06b0a-e5e8-4dbf-8bb0-729cce46a223_4a32dffd-584c-4f2a-8b17-e1e27b630f66)
// articNcovNanopore:Genotyping:typeVariants (61c06b0a-e5e8-4dbf-8bb0-729cce46a223_4a32dffd-584c-4f2a-8b17-e1e27b630f66)
// articNcovNanopore:Genotyping:mergeTypingCSVs (61c06b0a-e5e8-4dbf-8bb0-729cce46a223)
// articNcovNanopore:sequenceAnalysisMedaka:collateSamples (61c06b0a-e5e8-4dbf-8bb0-729cce46a223_4a32dffd-584c-4f2a-8b17-e1e27b630f66)
// articNcovNanopore:sequenceAnalysisMedaka:writeQCSummaryCSV (61c06b0a-e5e8-4dbf-8bb0-729cce46a223)
// ```
// Processes run by Illumina NCOV
// ncovIllumina:prepareReferenceFiles:indexReference (SARS-CoV-2.reference.fasta)
// ncovIllumina:sequenceAnalysis:readTrimming (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:readMapping (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:trimPrimerSequences (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:callVariants (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:makeConsensus (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:makeQCCSV (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:collateSamples (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:sequenceAnalysis:writeQCSummaryCSV (61c06b0a-e5e8-4dbf-8bb0-729cce46a223)
// ncovIllumina:Genotyping:typeVariants (15d3a409-d84e-4e7f-a832-5980956f8ee3)
// ncovIllumina:Genotyping:mergeTypingCSVs (61c06b0a-e5e8-4dbf-8bb0-729cce46a223)