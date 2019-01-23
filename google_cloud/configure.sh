#!/bin/bash

contains() {
    if [[ $1 =~ (^|[[:space:]])$2($|[[:space:]]) ]]
    then
        echo 0
    else
        echo 1
    fi
}

currentFolder=$(pwd)

## General information
instanceUsername="guidanceproject2018"
instanceName="guidancebase"
publicSSHfile="${HOME}/.ssh/id_rsa.pub"
projectName="guidance"
identificationJson="${currentFolder}/guidance-06986d1d6789.json"
project="guidance"


## Bucket information  
bucketLocation="us-east1"
bucketProjectName="${project}"
storageClass="regional"
bucketName="bucket-${bucketProjectName}"
zone=$(gcloud compute zones list | grep ${bucketLocation} | awk '{ print $1 }' | sort | head -n1)


## Instance information
baseInstanceName="guidancecluster"
## standard, highmem, highcpu
memoryInstanceConfig="highcpu"
#memoryInstanceConfig="standard"
## St: 1 2 4 8 16 32 64 96
## High{Mem-Cpu}: 2 4 8 16 32 64 96
cpuInstance="2"


## Snapshot name
snapName="snap${project}"

## Cluster information
amountOfNodes="3"

#echo "## MEMORY OPTION ##"
## Param verification
if [ ${memoryInstanceConfig} == "standard" ]
then
    #echo "Standard"
    correctOptions="1 2 4 8 16 32 64 96"
elif [ 0 == $(contains "highcpu highmem" "${memoryInstanceConfig}") ]
then
    #echo "High memory or high cpu"
    correctOptions="2 4 8 16 32 64 96"
else
    #echo "Memory not recognized: ${memoryInstanceConfig}"
    exit 1
fi

#echo "## AMOUNT OF CPUs PER NODE"
if [ 1 == $(contains "${correctOptions[@]}" "${cpuInstance}") ]
then
    #echo "Bad amount of CPUs: ${cpuInstance}"
    exit 1
else
    #echo ${cpuInstance}
    true
fi

machineType="n1-${memoryInstanceConfig}-${cpuInstance}"
