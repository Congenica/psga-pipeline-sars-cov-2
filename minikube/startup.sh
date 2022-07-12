#!/usr/bin/env bash

wait_for_pod() {
    local __POD="${1}"

    printf "waiting for pod ${__POD} to run "
    while [[ $(kubectl get pods ${__POD} -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
        printf "." && sleep 1
    done
    printf "\n${__POD} is running\n"
}

# label minikube node so that it matches against the cluster
kubectl label --overwrite node minikube farmNode=true

# create new namespace and set it as default
kubectl apply -f create_namespace.yaml
kubectl config set-context $(kubectl config current-context) --namespace=psga-minikube

# set service account
kubectl apply -f service_account.yaml

# set psga resources (e.g. pvc)
kubectl apply -f deploy_psga_resources.yaml

## deploy psga
kubectl apply -f deploy_psga.yaml
pipeline_pod="$( kubectl get pods -l app=psga-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
# set up a metadata.csv file containing the sample file paths
# copy your aws credentials so that you can fetch files from s3 within the pod
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/
