#!/usr/bin/env bash
set -euo pipefail # exit on any failures

source config.sh

PSGA_ROOT_PATH="/app"

for name in $PIPELINES; do
  for config in $(ls test_configs/$name); do
    echo "Running $name/$config test"
    source test_configs/$name/$config
    pipeline_pod="$( kubectl get pods -l app=$name-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
    #kubectl -n psga-minikube exec $pipeline_pod -- bash -c "AWS_SECRET_ACCESS_KEY=${AWS_SECRET_ACCESS_KEY} AWS_ACCESS_KEY_ID=${AWS_ACCESS_KEY_ID} nextflow -log ${ANALYSIS_RUN}.log run ${PSGA_ROOT_PATH}/psga/main.nf --run ${ANALYSIS_RUN} --sequencing_technology ${SEQUENCING_TECHNOLOGY} --kit ${KIT} --metadata ${METADATA} --output_path ${OUTPUT_PATH}"
    kubectl -n psga-minikube exec $pipeline_pod -- bash -c "cd /app/jenkins && pytest test_validation.py --expected-results-csv ${PSGA_ROOT_PATH}/jenkins/files/sars_cov_2/expected_results/pipeline_ci_illumina/results.csv --results-csv ${OUTPUT_PATH}/results.csv --output-path ${OUTPUT_PATH} --pathogen sars_cov_2 --sequencing-technology ${SEQUENCING_TECHNOLOGY} && cd -"
  done
done

# string(name: 'METADATA', value: "s3://synthetic-data-dev/UKHSA/small_tests/ont/metadata.csv"),
# string(name: 'OUTPUT_PATH', value: "/data/output/ont_none"),
# string(name: 'NXF_WORK', value: "/data/work/ont_none"),
# string(name: 'ANALYSIS_RUN', value: "ont_none"),
# string(name: 'SEQUENCING_TECHNOLOGY', value: "ont"),
# string(name: 'KIT', value: "none"),
#
# string(name: 'METADATA', value: "s3://synthetic-data-dev/UKHSA/requirements/CL0015b/ont_v4_input/metadata.csv"),
# string(name: 'OUTPUT_PATH', value: "/data/output/ont_artic_v4-1"),
# string(name: 'NXF_WORK', value: "/data/work/ont_artic_v4-1"),
# string(name: 'ANALYSIS_RUN', value: "ont_artic_v4-1"),
# string(name: 'SEQUENCING_TECHNOLOGY', value: "ont"),
# string(name: 'KIT', value: "unknown"),
#
# string(name: 'METADATA', value: "s3://synthetic-data-dev/UKHSA/small_tests/ont_midnight_v2/metadata.csv"),
# string(name: 'OUTPUT_PATH', value: "/data/output/ont_midnight_v2"),
# string(name: 'NXF_WORK', value: "/data/work/ont_midnight_v2"),
# string(name: 'ANALYSIS_RUN', value: "ont_midnight_v2"),
# string(name: 'SEQUENCING_TECHNOLOGY', value: "ont"),
# string(name: 'KIT', value: "Midnight-ONT_V2"),
#
# string(name: 'METADATA', value: "s3://synthetic-data-dev/UKHSA/small_tests/unknown/metadata.csv"),
# string(name: 'OUTPUT_PATH', value: "/data/output/unknown_short"),
# string(name: 'NXF_WORK', value: "/data/work/unknown_short"),
# string(name: 'ANALYSIS_RUN', value: "unknown_short"),
# string(name: 'SEQUENCING_TECHNOLOGY', value: "unknown"),
# string(name: 'KIT', value: "none"),
#
# string(name: 'METADATA', value: "s3://synthetic-data-dev/UKHSA/from_them/gold_standard/illumina/v4/2022-08-18/pipeline_metadata.csv"),
# string(name: 'OUTPUT_PATH', value: "/data/output/empty_sample_and_sample_failing_ncov_qc"),
# string(name: 'NXF_WORK', value: "/data/work/empty_sample_and_sample_failing_ncov_qc"),
# string(name: 'ANALYSIS_RUN', value: "empty_sample_and_sample_failing_ncov_qc"),
# string(name: 'SEQUENCING_TECHNOLOGY', value: "illumina"),
# string(name: 'KIT', value: "ARTIC_V4-1"),
