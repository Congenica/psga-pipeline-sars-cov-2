#!/bin/bash

set -e
set -o pipefail

source sars_cov_2.deps

# Organise the primer schemes (e.g. ARTIC, Midnight-IDT, Midnight-ONT, NEB-VarSkip) so that ncov can find the paths correctly

pathogen="SARS-CoV-2"
primers_dir="primer_schemes"

rm -rf ${primers_dir} wf-artic

git clone https://github.com/epi2me-labs/wf-artic.git
cd wf-artic
git checkout ${primer_schemes_commit}

# restructure the primer schemes in order to execute `artic minion` correctly to:
# /${primers_dir}/<SCHEME_NAME>/<PATHOGEN>/<SCHEME_VERSION>
for i in data/${primers_dir}/${pathogen}/*; do
    scheme=$(basename "$i")
    mkdir -p "../${primers_dir}/${scheme}"
    cp -R "$i" "../${primers_dir}/${scheme}"/"${pathogen}"

    # replace version containg a dot with hyphen (e.g. V4.1 => V4-1)
    for version_path in ../${primers_dir}/${scheme}/${pathogen}/*; do
        version=$(basename "${version_path}")
        if [[ "${version}" =~ \. ]]; then
            version_parent=$(dirname "${version_path}")
            version_hyphen="${version//\./-}"
            mv ${version_path} ${version_parent}/${version_hyphen}
        fi
    done
done
cd ../

rm -rf wf-artic
