#!/usr/bin/env bash
set -euo pipefail # exit on any failures

source config.sh

kubectl delete namespace $KUBE_NAMESPACE
