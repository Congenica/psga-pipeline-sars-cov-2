#!/usr/bin/env bash

## delete any nf pod
kubectl get pods -n psga-minikube --no-headers=true | awk '/nf/{print $1}'| xargs  kubectl delete -n psga-minikube pod

## delete psga-minikube
kubectl delete deployment psga-pipeline-minikube
kubectl delete pvc psga-minikube-pvc
kubectl delete rolebinding psga-minikube-admin
kubectl delete serviceaccount psga-minikube-admin
kubectl delete role psga-minikube-admin
kubectl delete namespace psga-minikube
