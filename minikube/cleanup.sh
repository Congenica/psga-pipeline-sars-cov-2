#!/usr/bin/env bash
set -euo pipefail # exit on any failures

source config.sh

## delete psga-minikube
for name in $PIPELINES; do
  kubectl delete -f pipelines/$name.yaml
done
kubectl delete pvc psga-minikube-pvc
kubectl delete rolebinding psga-minikube-admin
kubectl delete serviceaccount psga-minikube-admin
kubectl delete role psga-minikube-admin
kubectl delete namespace psga-minikube
