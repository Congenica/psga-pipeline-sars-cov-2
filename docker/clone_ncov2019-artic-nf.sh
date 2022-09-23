#!/bin/bash

set -e
set -o pipefail

source sars_cov_2.deps

git clone --recursive https://github.com/Congenica/ncov2019-artic-nf.git
cd ncov2019-artic-nf
git checkout ${ncov_commit}

cp environments/extras.yml ../ncov2019-artic-nf_extras.yml
cp environments/illumina/environment.yml ../ncov2019-artic-nf_illumina_environment.yml
cp environments/nanopore/environment.yml ../ncov2019-artic-nf_nanopore_environment.yml
cd ../

rm -rf ncov2019-artic-nf
