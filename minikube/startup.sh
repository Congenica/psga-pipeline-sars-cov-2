#!/usr/bin/env bash


wait_for_pod() {
    local __POD="${1}"

    printf "waiting for pod ${__POD} to run "
    while [[ $(kubectl get pods ${__POD} -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
        printf "." && sleep 2
    done
    printf "\n${__POD} is running\n"
}


## create rolebinding. e.g. kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=<namespace>:default
kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=default:default
kubectl create rolebinding default-view --clusterrole=view --serviceaccount=default:default

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
# copy some basic input files so that we do not need to copy over every time
kubectl exec -it ${pipeline_pod} -- bash -c 'mkdir -p /data/input && cp /app/minikube/input_files/* /data/input/'
