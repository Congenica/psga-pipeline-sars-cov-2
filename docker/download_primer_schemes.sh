#!/bin/bash

set -e
set -o pipefail

# Organise the primer schemes (e.g. ARTIC, Midnight-IDT, Midnight-ONT, NEB-VarSkip)
# This script is executed in a Dockerfile

# In the future, we might want to keep these primer-schemes in our own repository
# This will require some regular maintenance in order to make sure that the latest versions are available

# the commit to checkout
primer_schemes_commit="$1"
pathogen="$2"

git clone https://github.com/epi2me-labs/wf-artic.git wf-artic
cd wf-artic
ls -lart /wf-artic/data/primer_schemes/SARS-CoV-2
ls -lart /wf-artic/data/primer_schemes/${pathogen}
git checkout ${primer_schemes_commit}
cd -

# restructure the primer schemes in order to execute `artic minion` correctly to:
# /primer_schemes/<SCHEME_NAME>/<PATHOGEN>/<SCHEME_VERSION>
for i in wf-artic/data/primer_schemes/${pathogen}/*; do
    scheme=$(basename "$i")
    mkdir -p "primer_schemes/${scheme}"
    mv "$i" "primer_schemes/${scheme}"/"${pathogen}"
done

rm -rf wf-artic
