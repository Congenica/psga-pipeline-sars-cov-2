#!/bin/bash

set -e
set -o pipefail

pathogen="sars_cov_2"

source ${pathogen}/ncov.deps

git clone --recursive https://github.com/Congenica/ncov2019-artic-nf.git
cd ncov2019-artic-nf
git checkout ${ncov_commit}

cp environments/extras.yml ../${pathogen}/ncov2019-artic-nf_extras.yml
cp environments/illumina/environment.yml ../${pathogen}/ncov2019-artic-nf_illumina.yml
cp environments/nanopore/environment.yml ../${pathogen}/ncov2019-artic-nf_nanopore.yml

cd ../

rm -rf ncov2019-artic-nf
