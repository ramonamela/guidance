#!/bin/bash

source configure.sh
currentIP=$(gcloud compute instances list --filter="${instanceName}" | tail -n1 | awk '{print $5}')
scp installGuidanceDependencies.sh ${instanceUsername}@${currentIP}:/home/${instanceUsername}
ssh ${instanceUsername}@${currentIP} "sh /home/${instanceUsername}/installGuidanceDependencies.sh"
ssh ${instanceUsername}@${currentIP} "rm /home/${instanceUsername}/installGuidanceDependencies.sh"
