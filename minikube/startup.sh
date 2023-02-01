#!/usr/bin/env bash

wait_for_pod() {
    local __POD="${1}"
    echo "Waiting for pod ${__POD} to run"
    while [[ $(kubectl get pods ${__POD} -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
      kubectl describe pod ${__POD}
      sleep 1
    done
    echo "${__POD} is running"
}

echo "Labeling minikube node so that it matches against the cluster"
kubectl label --overwrite node minikube farmNode=true

echo "Creating new namespace and setting it as default"
kubectl apply -f create_namespace.yaml
kubectl config set-context $(kubectl config current-context) --namespace=psga-minikube

echo "Setting service account"
kubectl apply -f service_account.yaml

echo "Setting psga resources (e.g. pvc)"
kubectl apply -f deploy_psga_resources.yaml

echo "Listing running pods"
kubectl get pods

echo "Deploying sars-cov-2 pipeline"
kubectl apply -f deploy_sars_cov_2_pipeline.yaml
echo "Waiting for the sars-cov-2-pipeline-minikube pod to be ready"
pipeline_pod="$( kubectl get pods -l app=sars-cov-2-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/

echo "Deploying synthetic pipeline"
kubectl apply -f deploy_synthetic_pipeline.yaml
echo "Waiting for the synthetic-pipeline-minikube pod to be ready"
pipeline_pod="$( kubectl get pods -l app=synthetic-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/

echo "Deploying s-aureus pipeline"
kubectl apply -f deploy_s_aureus_pipeline.yaml
echo "Waiting for the s-aureus-pipeline-minikube pod to be ready"
pipeline_pod="$( kubectl get pods -l app=s-aureus-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )"
wait_for_pod "${pipeline_pod}"
kubectl cp ${HOME}/.aws ${pipeline_pod}:/root/
