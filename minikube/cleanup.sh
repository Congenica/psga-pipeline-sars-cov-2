#!/usr/bin/env bash


## delete DB
kubectl delete deployment psql-minikube
kubectl delete service psql-minikube-host
kubectl delete pvc psql-minikube-pvc

## delete psga-minikube
kubectl delete deployment psga-minikube
kubectl delete pvc psga-minikube-pvc
kubectl delete rolebinding psga-minikube-admin
kubectl delete serviceaccount psga-minikube-admin
kubectl delete role psga-minikube-admin
kubectl delete namespace ukhsa-covid-minikube


## delete any nf pod
#kubectl get pods -n default --no-headers=true | awk '/nf/{print $1}'| xargs  kubectl delete -n default pod
