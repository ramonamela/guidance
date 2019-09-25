#!/bin/bash

source configure.sh

gcloud auth activate-service-account --key-file=${identificationJson}
gcloud config set project ${projectName}
availZone=$(gcloud compute zones list | grep ${bucketLocation} | grep UP | awk '{print $1 }' | head -n1)
gcloud config set compute/zone ${availZone}
