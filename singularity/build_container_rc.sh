#!/bin/bash

rm -f singularity.img

#sudo docker run -d -p 5000:5000 --name registry registry:2
#sudo docker pull docker_guidance
#sudo docker image tag docker_guidance localhost:5000/docker_guidance_image
#sudo docker push localhost:5000/docker_guidance_image

#export SINGULARITY_NOHTTPS=true

## DOCKER IMAGE GENERATION ##

sudo docker ps -a | tail -n +2 | awk '{ print $1 }' | xargs -i sudo docker stop {}
sudo docker ps -a | tail -n +2 | awk '{ print $1 }' | xargs -i sudo docker rm {}

sudo docker image ls | grep "docker_guidance" | awk '{ print $3 }' | tail -n +2 | xargs -i sudo docker rmi -f {}

pushd docker
./build_docker_rc.sh
popd

#sudo docker stop registry
#sudo docker ps -a | grep registry:2 | xargs -i sudo docker rm {}
sudo docker rmi -f registry:2
sudo docker run -d -p 5000:5000 --restart=always --name registry registry:2
sudo docker tag docker_guidance_rc localhost:5000/docker_guidance_rc
sudo docker push localhost:5000/docker_guidance_rc

## SINGULARITY IMAGE GENERATION ##

# https://github.com/singularityware/singularity/issues/429 

sudo SINGULARITY_NOHTTPS=yes singularity build guidance_singularity.img singularity/guidance_rc.def
