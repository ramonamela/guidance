#!/bin/bash

rm -f singularity.img

#sudo docker run -d -p 5000:5000 --name registry registry:2
#sudo docker pull docker_guidance
#sudo docker image tag docker_guidance localhost:5000/docker_guidance_image
#sudo docker push localhost:5000/docker_guidance_image

#export SINGULARITY_NOHTTPS=true

## DOCKER IMAGE GENERATION ##

pushd docker
./build_docker.sh
popd

sudo docker run -d -p 5000:5000 --restart=always --name registry registry:2
sudo docker tag docker_guidance localhost:5000/docker_guidance
sudo docker push localhost:5000/docker_guidance

## SINGULARITY IMAGE GENERATION ##

# https://github.com/singularityware/singularity/issues/429 

sudo SINGULARITY_NOHTTPS=yes singularity build singularity.img singularity/guidance.def
