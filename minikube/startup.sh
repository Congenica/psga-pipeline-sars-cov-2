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

## deploy sars-cov-2 pipeline
kubectl apply -f deploy_sars_cov_2.yaml
pipeline_pod="$( kubectl get pods -l app=sars-cov-2-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/
