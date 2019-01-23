#!/bin/bash

source configure.sh
currentIP=$(gcloud compute instances list --filter="${instanceName}" | tail -n1 | awk '{print $5}')
scp -o "StrictHostKeyChecking no" installGuidanceDependencies.sh ${instanceUsername}@${currentIP}:/home/${instanceUsername}
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "rm -rf ~/TOOLS"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "sh /home/${instanceUsername}/installGuidanceDependencies.sh"
ssh -o "StrictHostKeyChecking no" ${instanceUsername}@${currentIP} "rm /home/${instanceUsername}/installGuidanceDependencies.sh"
