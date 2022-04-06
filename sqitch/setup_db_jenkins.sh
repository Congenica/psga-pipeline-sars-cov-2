#!/usr/bin/env bash

# This script is executed by Jenkins
# Deploy DB schema migrations

SQITCH_URI="db:pg://${DB_USER}:${DB_PASSWORD}@${DB_HOST}:${DB_PORT}/${DB_NAME}"

sqitch deploy ${SQITCH_URI}
sqitch verify ${SQITCH_URI}
