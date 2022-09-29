#!/usr/bin/env bash

## delete any nf pod
kubectl get pods -n psga-minikube --no-headers=true | awk '/nf/{print $1}'| xargs  kubectl delete -n psga-minikube pod
kubectl get jobs -n psga-minikube --no-headers=true | awk '/nf/{print $1}'| xargs  kubectl delete -n psga-minikube job

## delete psga-minikube
kubectl delete deployment sars-cov-2-pipeline-minikube
kubectl delete pvc psga-minikube-pvc
kubectl delete rolebinding psga-minikube-admin
kubectl delete serviceaccount psga-minikube-admin
kubectl delete role psga-minikube-admin
kubectl delete namespace psga-minikube
