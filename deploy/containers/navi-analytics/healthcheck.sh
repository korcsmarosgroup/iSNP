#!/bin/bash
set -e

if [ -e /executor/executor_healthy ]
then
    echo "executor service started"
    exit 0
else
    echo "executor service not started yet"
    exit 1
fi
