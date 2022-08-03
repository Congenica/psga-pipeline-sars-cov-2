#!/bin/bash

# Copy reports, logs and workdir logs for this analysis run to the output directory
# To enable the main log and reports, invoke nextflow as follows:
# nextflow -log <ANALYSIS_RUN>.log run . -with-trace <PIPELINE PARAMETERS>

analysis_run=$1
output_path=$2

[ -z "${analysis_run}" ] && echo "Error: analysis run is unset" && exit 1
[ -z "${output_path}" ] && echo "Error: output path is unset" && exit 1

# clean nextflow work directory (env var: NXF_WORK), but retain the log files for each process
nextflow clean -force -keep-logs -quiet

# create a compressed archive
# NOTE: to speed up the tests, don't create / save this
#tar czf ${analysis_run}_workdir.tar.gz ${NXF_WORK}

if [[ "${output_path}" =~ ^s3.* ]]; then
    # copy to s3
    aws s3 cp . ${output_path} --exclude "*" --include "${analysis_run}*.*" --recursive
else
    cp ${analysis_run}*.* ${output_path}
fi

# cleanup
rm -f ${analysis_run}*.*
rm -rf ${NXF_WORK}
