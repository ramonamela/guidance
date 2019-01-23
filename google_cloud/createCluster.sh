#!/bin/bash

source ./configure.sh

./removeCluster.sh

./createBaseInstance.sh

./createSnapshot.sh

serviceAccount=$(gcloud compute instances describe ${instanceName} | grep service | awk '{ print $3 }')
machineType=$(gcloud compute instances describe ${instanceName} | grep machine | awk '{ print $2 }' | tr "/" "\t" | awk '{ print $NF }')

for ((i=1;i<=amountOfNodes;++i)); do
    currentName="${baseInstanceName}$(printf %04d $i)"
    gcloud compute disks create ${currentName} --source-snapshot "${snapName}"
    gcloud compute instances create "${currentName}" --machine-type ${machineType} --service-account ${serviceAccount} --disk "name=${currentName},device-name=${currentName},mode=rw,boot=yes,auto-delete=yes"
    ./addPublicKeyGoogleCloud.sh ${currentName}
done

