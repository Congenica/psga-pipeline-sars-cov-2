#!/usr/bin/env bash
set -euo pipefail # exit on any failures

source config.sh

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
  echo "Waiting for the $app_name deployment to be ready"
  kubectl wait deployment -n psga-minikube $app_name --for condition=Available=True --timeout=300s
done
