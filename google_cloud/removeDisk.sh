#!/bin/bash

if gcloud compute disks list | grep -q ${1}; then
    gcloud compute disks delete ${1} -q
fi
