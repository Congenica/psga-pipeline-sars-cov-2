#!/usr/bin/env bash

# this is executed within the k8s pod containing the DB
createdb covid_pipeline_db
sqitch deploy
