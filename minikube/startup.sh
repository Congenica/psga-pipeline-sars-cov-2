#!/usr/bin/env bash

wait_for_pod() {
    local __POD="${1}"

    printf "waiting for pod ${__POD} to run "
    while [[ $(kubectl get pods ${__POD} -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
        printf "." && sleep 1
    done
    printf "\n${__POD} is running\n"
}

# create new namespace and set it as default
kubectl apply -f create_namespace.yaml
kubectl config set-context $(kubectl config current-context) --namespace=ukhsa-covid-minikube

# set service account
kubectl apply -f service_account.yaml

# set psga resources (e.g. pvc)
kubectl apply -f deploy_psga_resources.yaml

## deploy DB. Currently, the DB is stored in a pod.
# This is fine for a proof of concept, but obviously not the long term solution.
# The DB will be stored in a RDS aurora system
kubectl apply -f deploy_db.yaml
db_pod="$( kubectl get pods -l app=psql --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${db_pod}"
kubectl exec -it ${db_pod} -- bash -c './setup_db.sh'

## deploy psga
kubectl apply -f deploy_psga.yaml
pipeline_pod="$( kubectl get pods -l app=psga-minikube --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
# you will need to copy files to the /data/input or change PSGA_INPUT_PATH to point to an s3 location
kubectl exec -it ${pipeline_pod} -- bash -c 'mkdir -p /data/input'
# copy your aws credentials so that you can fetch files from s3 within the pod
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/
