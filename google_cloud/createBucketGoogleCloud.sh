#!/bin/bash

## https://cloud.google.com/storage/docs/creating-buckets ##

source ./configure.sh
echo "${bucketProjectName}"
echo "${storageClass}"
echo "${bucketLocation}"
echo "gs://${bucketName}/"
if ! gsutil ls | grep -q -w ${bucketName}; then
    gsutil mb -p "${bucketProjectName}" -c "${storageClass}" -l "${bucketLocation}" "gs://${bucketName}/"
fi
## sudo apt-get install google-cloud-sdk
