#!/bin/bash

source ./configure.sh

echo "Init google cloud"
./initGCloud.sh

echo "Remove base instance"
./removeBaseInstance.sh

echo "Create new base instance"
gcloud compute instances create ${instanceName} --zone=${zone} --machine-type=${machineType} --image=ubuntu-minimal-1810-cosmic-v20190108 --image-project=ubuntu-os-cloud --boot-disk-type=pd-standard --boot-disk-size=10GB --boot-disk-device-name=${instanceName}

echo "Add public key google cloud"
./addPublicKeyGoogleCloud.sh ${instanceName}

echo "Wait until the image is running"
currentIP=$(gcloud compute instances list --filter="${instanceName}" | tail -n1 | awk '{print $5}')
ssh -q -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "exit"
while [ ! "$?" -eq "0" ]; do
    echo "Could not stablish connection, retry"
    sleep 1
    currentIP=$(gcloud compute instances list --filter="${instanceName}" | tail -n1 | awk '{print $5}')
    ssh -q -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "exit"
done

echo "The image is already running"
## Create and mount bucket
./createBucketGoogleCloud.sh
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "mkdir -p \$HOME/${bucketName}"
scp -o "StrictHostKeyChecking no" installGFuse.sh ${instanceUsername}@${currentIP}:\$HOME
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "sh \$HOME/installGFuse.sh"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "rm -f \$HOME/installGFuse.sh"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "gcsfuse bucket-guidance \$HOME/${bucketName}"

echo "Install guidance dependencies"
./installGuidanceDependenciesGoogle.sh

echo "Install compss"
./installCOMPSsGoogle.sh
