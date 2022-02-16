#!/usr/bin/env bash


## delete DB
kubectl delete deployment psql-minikube
kubectl delete service psql-minikube-host
kubectl delete pvc psql-minikube-pvc

## delete covid-pipeline-minikube
kubectl delete deployment covid-pipeline-minikube
kubectl delete pvc covid-pipeline-minikube-pvc
kubectl delete rolebinding covid-pipeline-minikube-admin
kubectl delete serviceaccount covid-pipeline-minikube-admin
kubectl delete role covid-pipeline-minikube-admin
kubectl delete namespace ukhsa-covid-minikube


## delete any nf pod
#kubectl get pods -n default --no-headers=true | awk '/nf/{print $1}'| xargs  kubectl delete -n default pod
