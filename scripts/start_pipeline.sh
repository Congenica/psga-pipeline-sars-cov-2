#!/bin/bash

#
# This script can be run by cron to start the pipeline.
#
# Ensure COVID_PIPELINE_FASTQ_PATH is set pointing to a directory, the script assumes that there
# will be directories under here that will each contain the FASTQ files and metadata.tsv. If this
# script sees a file with a name starting "transfer_complete", it will create a "pipeline_started"
# file in the same directory, start a pipeline run for that directory and then exit. It will not
# start a pipeline if the "pipeline_started" file already exists.
#
# Note that the "transfer_complete" name is checked case-insensitively, the words can be separated
# by either a hyphen or an underscore, and the file can have any extension.
#

shopt -s nocaseglob extglob nullglob

PIPELINE_STARTED_FLAG_FILE="pipeline_started"
PIPELINE_COMPLETE_FLAG_FILE="pipeline_complete"

if [ -z "${COVID_PIPELINE_FASTQ_PATH}" ]
then
  echo "COVID_PIPELINE_FASTQ_PATH not set!"
  exit
fi

complete_files=("${COVID_PIPELINE_FASTQ_PATH}"/*/transfer[-_]complete*)

for complete_file in "${complete_files[@]}"
do
  dir=$(dirname "$complete_file")
  echo "Found a transfer_complete file in ${dir}"

  if [ -e "${dir}/${PIPELINE_STARTED_FLAG_FILE}" ]
  then
    echo "Pipeline already started for ${dir}"
    continue
  elif [ -e "${dir}/${PIPELINE_COMPLETE_FLAG_FILE}" ]
  then
    echo "Pipeline already run for for ${dir}"
    continue
  else
    outfile="nextflow_$(date +%s).out"
    echo "Starting pipeline for ${dir}, output to ${outfile}"
    touch "${dir}/${PIPELINE_STARTED_FLAG_FILE}"
    export COVID_PIPELINE_FASTQ_PATH=${dir}
    nextflow run -ansi-log false covid-pipeline > "${outfile}" 2>&1 &
    break
  fi
done
