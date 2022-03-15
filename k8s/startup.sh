#!/usr/bin/env bash

kubectl config set-context $(kubectl config current-context) --namespace=ukhsa-covid

kubectl apply -f k8s_service_account.yaml
kubectl apply -f k8s_deploy_psga_resources.yaml
kubectl apply -f k8s_deploy_psga.yaml
