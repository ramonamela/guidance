#!/bin/bash

source ./configure.sh

if gcloud compute instances list | grep -q -w ${1}; then
    gcloud compute instances delete ${1} -q &
fi

