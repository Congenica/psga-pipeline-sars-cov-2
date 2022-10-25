/*
 * Run: primer-autodetection (use first read only for illumina)
 */
process primer_autodetection {
  publishDir "${params.output_path}/primer_autodetection", mode: 'copy', overwrite: true, pattern: '{*_primer_data.csv,*_primer_detection.tsv,trimmomatic.out,*.bowtie2*}'

  tag "${task.index} - ${fastq}"

  input:
    path fastq

  output:
    tuple path("*_primer_*.txt"), path(fastq), emit: ch_files
    path "*_primer_data.csv", emit: ch_primer_data
    path "*_primer_detection.tsv", emit: ch_primer_coverage
    path "trimmomatic.out"
    path "*.bowtie2*"

  shell:
  '''
  input_primer_val=!{params.kit}

  # select the input fastq file based on illumina or ont.
  # For illumina, 1 read is sufficient for primer autodetection
  if [[ "!{params.sequencing_technology}" == "illumina" ]]; then
    file_1="$(ls *_1.fastq.gz)"
  elif [[ "!{params.sequencing_technology}" == "ont" ]]; then
    file_1="!{fastq}"
  else
    echo "sequencing technology must be either 'illumina' or 'ont' for primer autodetection"
    exit 1
  fi

  # extract the sample id, whereas this is ID.fastq.gz or ID_1.fastq.gz
  sample_id=$( echo ${file_1} | cut -d '.' -f1 | cut -d '_' -f1)
  trimmed_sample="trimmed_sample.fastq.gz"
  crop="!{params.primer_minimum_length}"

  # trim sample reads so that they are long as the shortest primer
  trimmomatic SE -phred33 ${file_1} ${trimmed_sample} CROP:${crop} 2> trimmomatic.out

  # iterate over the available primer scheme fasta.
  # These fasta are generated by us and contain the sequences for each primer scheme name/version
  readarray -t primer_fasta_array < "/primer_schemes/SARS-CoV-2_primer_fasta_index.txt"
  for primer_fasta in "${primer_fasta_array[@]}"; do
      # extract the primer (e.g. ARTIC_V3)
      primer="$(echo "${primer_fasta}" | cut -d/ -f3,5 | sed 's#/#_#')"

      # build the index for the primer fasta
      bowtie2-build ${primer_fasta} ${primer}.ref_idx &> ${primer}.bowtie2-build

      # end-to-end align the trimmed sample to the primer scheme fasta
      bowtie2 --no-unal -x ${primer}.ref_idx -U ${trimmed_sample} 2> ${primer}.bowtie2 | samtools sort -o ${primer}.bam

      # extract the coverage and other metrics from the bam. `samtools coverage` complains if the bam is not sorted
      samtools coverage ${primer}.bam > ${primer}.coverage.tsv
  done

  ###############################
  # THE FOLLOWING CODE WILL BE REFACTORED TO PYTHON IN https://jira.congenica.net/browse/PSG-1275
  ###############################

  # concatenate the coverage files to 1 single file with 1 single header
  awk 'FNR==1 && NR!=1{next;}{print}' *.coverage.tsv > ${sample_id}_primer_detection.tsv

  # extract 1) the primer with the highest number of reads and 2) the number of reads
  detected_primer="$(awk 'BEGIN{nreads=0;primer="none"}{if(NR>1 && $4>nreads) {nreads=$4; primer=$1}} END {print primer,nreads}' ${sample_id}_primer_detection.tsv)"
  primer="$(echo "${detected_primer}" | awk '{print $1}')"
  primer_nreads="$(echo "${detected_primer}" | awk '{print $2}')"

  # store the data for the selected primer. Use CSV for convenience here
  awk -v name="${primer}" 'NR==1 || $1==name {$1=$1; print}' OFS=, ${sample_id}_primer_detection.tsv > ${sample_id}_primer_data.csv.tmp

  # add additional columns
  if [[ "${input_primer_val}" == "unknown" ]]; then
    if [[ ${primer_nreads} -gt 0 ]]; then
      qc_val="PASS"
    else
      qc_val="FAIL"
    fi
  elif [[ "${input_primer_val}" == "none" ]]; then
    if [[ ${primer_nreads} -gt 0 ]]; then
      qc_val="FAIL"
    else
      qc_val="PASS"
    fi
  else
    if [[ "${input_primer_val}" == "${primer}" ]]; then
      qc_val="PASS"
    else
      qc_val="FAIL"
    fi
  fi
  awk -v s="sample_id" -v qc="primer_qc" -v ip="primer_input" -v s_val="${sample_id}" -v qc_val="${qc_val}" -v ip_val="${input_primer_val}" \
    'BEGIN{FS=OFS=","} {print (NR==1 ? s OFS qc OFS ip : s_val OFS qc_val OFS ip_val), $0}' \
    ${sample_id}_primer_data.csv.tmp > ${sample_id}_primer_data.csv

  # reheader
  new_header="sample_id,primer_qc,primer_input,primer_detected,startpos,endpos,primer_numreads,primer_covbases,primer_coverage,meandepth,meanbaseq,meanmapq"
  sed -i "1s/.*/${new_header}/" ${sample_id}_primer_data.csv

  # store the primer scheme name/version
  echo "${primer}" > ${sample_id}_primer_${qc_val}.txt

  # add a default line if no primer was autodetected
  if [[ "${primer}" == "none" ]]; then
      echo "${sample_id},${qc_val},${input_primer_val},${primer},,,0,0,0,0,0,0" >> ${sample_id}_primer_data.csv

      # WARNING: this is a hack
      # certain samples could not have primers at all (e.g. bam files in which adapters and primers were trimmed in the past)
      # ncov and artic expect a primer scheme and version as input parameters, though.
      # here we mock ncov by passing a primer scheme/version to make it run.
      # This won't have consequences because none of the sequences in this primer was found in the sample,
      # so they won't be trimmed
      if [[ "${input_primer_val}" == "none" ]]; then
          echo "ARTIC_V3" > ${sample_id}_primer_${qc_val}.txt
      fi
  fi
  '''
}