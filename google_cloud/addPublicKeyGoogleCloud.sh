#!/bin/bash

## CONFIGURE SSH KEYS INTO THE NEW INSTANCE ##

currentInstanceName=$1

#gcloud init
source ./configure.sh
gcloud compute instances stop ${currentInstanceName}
gcloud compute instances describe ${currentInstanceName} | grep ssh-rsa | xargs -i echo {} > /tmp/publicKeys.txt
currentUsername="$(cat ${publicSSHfile} | awk '{ print $3 }')" && currentPublicKey="$(cat ${publicSSHfile})" && echo "${instanceUsername}:${currentPublicKey}" >> /tmp/publicKeys.txt
gcloud compute instances add-metadata ${currentInstanceName} --metadata-from-file ssh-keys=/tmp/publicKeys.txt
gcloud compute instances start ${currentInstanceName}
rm /tmp/publicKeys.txt
currentIP=$(gcloud compute instances list --filter="${currentInstanceName}" | tail -n1 | awk '{print $5}')
echo $currentIP
