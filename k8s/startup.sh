#!/usr/bin/env bash

kubectl config set-context $(kubectl config current-context) --namespace=psga

kubectl apply -f k8s_namespace.yaml
kubectl apply -f k8s_service_account.yaml
kubectl apply -f k8s_psga_resources.yaml
kubectl apply -f k8s_sars_cov_2.yaml
