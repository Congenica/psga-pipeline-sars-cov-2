#!/usr/bin/env bash


wait_for_pod() {
    local __POD="${1}"

    printf "waiting for pod ${__POD} to run "
    while [[ $(kubectl get pods ${__POD} -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
        printf "." && sleep 1
    done
    printf "\n${__POD} is running\n"
}

# set new namespace
kubectl apply -f create_namespace.yaml
kubectl config set-context $(kubectl config current-context) --namespace=ukhsa-covid

# set service account
kubectl apply -f service_account.yaml

# I think this is no longer required as we use the same service account config as in congenica k8s cluster
## create rolebinding. e.g. kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=<namespace>:default
#kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=ukhsa-covid:ukhsa-covid-admin
#kubectl create rolebinding default-view --clusterrole=view --serviceaccount=ukhsa-covid:ukhsa-covid-admin


## deploy DB. Currently, the DB is stored in a pod.
# This is fine for a proof of concept, but obviously not the long term solution.
# The DB will be stored in a RDS aurora system
kubectl create -f deploy_db.yaml
db_pod="$( kubectl get pods -l app=psql --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${db_pod}"
kubectl exec -it ${db_pod} -- bash -c './setup_db.sh'

## deploy covid-pipeline
kubectl create -f deploy_covid_pipeline.yaml
pipeline_pod="$( kubectl get pods -l app=covid-pipeline --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
# you will need to copy files to the /data/input or change COVID_PIPELINE_INPUT_PATH to point to an s3 location
kubectl exec -it ${pipeline_pod} -- bash -c 'mkdir -p /data/input'
# copy your aws credentials so that you can fetch files from s3 within the pod
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/
