#!/usr/bin/env bash
set -euo pipefail # exit on any failures

source config.sh

wait_for_pod() {
    POD="$1"
    while [[ $(kubectl -n psga-minikube get pods $POD -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
      echo "Waiting for pod $POD to run"
      kubectl get pvc
      kubectl describe pod $POD
      sleep 5
    done
    echo "$POD is running"
}

echo "Labeling minikube node so that it matches against the cluster"
kubectl label --overwrite node minikube farmNode=true

echo "Creating new namespace and setting it as default"
kubectl apply -f create_namespace.yaml
kubectl config set-context minikube --namespace=psga-minikube

echo "Setting service account"
kubectl apply -f service_account.yaml

# this raises an error if already present
kubectl create secret generic regcred \
  --from-file=.dockerconfigjson=$HOME/.docker/config.json \
  --type=kubernetes.io/dockerconfigjson
kubectl patch serviceaccount psga-minikube-admin -p '{"imagePullSecrets": [{"name": "regcred"}]}'

echo "Setting psga resources (e.g. pvc)"
kubectl apply -f deploy_psga_resources.yaml

echo "Listing running pods"
kubectl get pods

for name in $PIPELINES; do
  echo "Deploying $name pipeline"
  kubectl apply -f pipelines/$name.yaml
done

for name in $PIPELINES; do
  app_name=$name-pipeline-minikube
  echo "Waiting for the $app_name pod to be ready"
  kubectl -n psga-minikube get pods -l app=$app_name --no-headers -o custom-columns=':metadata.name'
  pipeline_pod=$( kubectl -n psga-minikube get pods -l app=$name-pipeline-minikube --no-headers -o custom-columns=':metadata.name' )
  echo "pipeline_pod: $pipeline_pod"
  wait_for_pod "$pipeline_pod"
done
