#!/bin/bash

source configure.sh
currentIP=$(gcloud compute instances list --filter="${instanceName}" | tail -n1 | awk '{print $5}')
scp -o "StrictHostKeyChecking no" installCOMPSsJava.sh ${instanceUsername}@${currentIP}:/home/${instanceUsername}
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "sudo rm -rf ~/2.4"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "sh /home/${instanceUsername}/installCOMPSsJava.sh"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "rm /home/${instanceUsername}/installCOMPSsJava.sh"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "sudo rm -rf ~/2.4"
