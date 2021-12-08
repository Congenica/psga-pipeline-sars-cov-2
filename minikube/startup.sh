#!/bin/bash


## create rolebinding. e.g. kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=<namespace>:default
kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=default:default
kubectl create rolebinding default-view --clusterrole=view --serviceaccount=default:default

## create DB. Currently, the DB is stored in a pod. This is fine for a proof of concept, but obviously not the long term solution. The DB will be stored in a RDS aurora system
kubectl create -f deploy_db.yaml
db_pod="$( kubectl get pods -l app=psql --no-headers -o custom-columns=':metadata.name' )"

while [[ $(kubectl get pods ${db_pod} -o 'jsonpath={..status.conditions[?(@.type=="Ready")].status}') != "True" ]]; do
    echo "waiting for pod ${db_pod} to run" && sleep 2
done
kubectl exec -it ${db_pod} -- bash -c './setup_db.sh'


kubectl create -f deploy_covid_pipeline.yaml
