#!/bin/bash

set -e
set -o pipefail

source sars_cov_2.deps

git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
git checkout ${pangolin_commit}
cp environment.yml ../pangolin_environment.yml
cd ../

rm -rf pangolin
