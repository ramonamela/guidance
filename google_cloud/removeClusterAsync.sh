#!/bin/bash

source ./configure.sh

#serviceAccount=$(gcloud compute instances describe ${instanceName} | grep service | awk '{ print $3 }')
#machineType=$(gcloud compute instances describe ${instanceName} | grep machine | awk '{ print $2 }' | tr "/" "\t" | awk '{ print $NF }')

## Look into the gcloud compute instances create --image flag

#gcloud compute instances delete "${baseInstanceName}${index}" -q
#gcloud compute disks delete "${baseInstanceName}${index}" -q

#gcloud compute disks create "${baseInstanceName}${index}" --source-snapshot "${snapName}" 
#gcloud compute instances create "${baseInstanceName}${index}" --machine-type ${machineType} --service-account ${serviceAccount} --disk "name=${baseInstanceName}${index},device-name=${baseInstanceName}${index},mode=rw,boot=yes,auto-delete=yes"

for ((i=1;i<=amountOfNodes;++i)); do
    currentName="${baseInstanceName}$(printf %04d $i)"
    ./removeInstanceAsync.sh ${currentName}
    #gcloud compute instances delete ${currentName} -q
    #./removeDisk.sh ${currentName}
    #gcloud compute disks delete ${currentName} -q
done

#gcloud compute disks create ${allNames} --source-snapshot "${snapName}" 
#gcloud compute instances create "${baseInstanceName}${index}" --machine-type ${machineType} --service-account ${serviceAccount} --source-instance-template ${instanceName}
