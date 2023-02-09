#!/usr/bin/env bash
set -euo pipefail # exit on any failures

source config.sh

PSGA_ROOT_PATH="/app"

for name in $PIPELINES; do
  for config in $(ls test_configs/$name); do
    echo "Running $name/$config test"
    source test_configs/$name/$config
    pipeline_pod="$( kubectl get pods -l app=$name-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
    kubectl -n psga-minikube exec $pipeline_pod -- bash -c "AWS_SECRET_ACCESS_KEY=${AWS_SECRET_ACCESS_KEY} AWS_ACCESS_KEY_ID=${AWS_ACCESS_KEY_ID} nextflow -log ${ANALYSIS_RUN}.log run ${PSGA_ROOT_PATH}/psga/main.nf --run ${ANALYSIS_RUN} --sequencing_technology ${SEQUENCING_TECHNOLOGY} --kit ${KIT} --metadata ${METADATA} --output_path ${OUTPUT_PATH}"
    kubectl -n psga-minikube exec $pipeline_pod -- bash -c "cd ${PSGA_ROOT_PATH}/jenkins && pytest test_validation.py --expected-results-csv ${PSGA_ROOT_PATH}/jenkins/files/sars_cov_2/expected_results/${ANALYSIS_RUN}/results.csv --results-csv ${OUTPUT_PATH}/results.csv --output-path ${OUTPUT_PATH} --pathogen sars_cov_2 --sequencing-technology ${SEQUENCING_TECHNOLOGY} && cd -"
  done
done
