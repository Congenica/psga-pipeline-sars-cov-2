#!/usr/bin/env bash
set -euxo pipefail # exit on any failures

source config.sh

echo "Labeling minikube node so that it matches against the cluster"
kubectl label --overwrite node minikube farmNode=true

echo "Making sure namespace $KUBE_NAMESPACE exists and it's set as default for kubectl context"
kubectl apply -f create_namespace.yaml
kubectl config set-context minikube --namespace=$KUBE_NAMESPACE

echo "Setting up service account"
kubectl apply -f service_account.yaml

if kubectl get secret | grep -q regcred; then
  echo "Docker secret exists, removing it"
  kubectl delete secret generic regcred
else
  echo "Creating docker secret for ECR pulls within minikube"
  kubectl create secret generic regcred \
    --from-file=.dockerconfigjson=$HOME/.docker/config.json \
    --type=kubernetes.io/dockerconfigjson
  kubectl patch serviceaccount psga-minikube-admin -p '{"imagePullSecrets": [{"name": "regcred"}]}'
fi

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
  kubectl wait deployment $app_name --for condition=Available=True --timeout=300s
done
