#!/bin/bash

rm -f guidance_singularity.img

guidance_previous_registry_name=registry
guidance_registry_name=guidance_registry
guidance_previous_image_name=docker_guidance
guidance_image_name=${guidance_previous_image_name}
registry_port=5000

#sudo docker run -d -p 5000:5000 --name registry registry:2
#sudo docker pull docker_guidance
#sudo docker image tag docker_guidance localhost:5000/docker_guidance_image
#sudo docker push localhost:5000/docker_guidance_image

#export SINGULARITY_NOHTTPS=true

## DOCKER IMAGE GENERATION ##

sudo docker ps -a | grep ${guidance_previous_registry_name} | awk '{ print $1 }' | xargs -i sudo docker stop {}
sudo docker ps -a | grep ${guidance_previous_registry_name} | awk '{ print $1 }' | xargs -i sudo docker rm {}

sudo docker image ls | grep ${guidance_previous_image_name} | awk '{ print $3 }' | tail -n +2 | xargs -i sudo docker rmi -f {}

pushd docker
./build_docker.sh ${guidance_image_name}
popd

#sudo docker stop registry
#sudo docker ps -a | grep registry:2 | xargs -i sudo docker rm {}
sudo docker rmi -f registry:2
sudo docker run -d -p ${registry_port}:${registry_port} --restart=always --name ${guidance_registry_name} registry:2
sudo docker tag ${guidance_image_name} localhost:${registry_port}/${guidance_image_name}
sudo docker push localhost:${registry_port}/${guidance_image_name}

## SINGULARITY IMAGE GENERATION ##

# https://github.com/singularityware/singularity/issues/429 

sed "s/\${docker_image_name}/${guidance_image_name}/g" singularity/base.def | tee singularity/guidance.def
sudo SINGULARITY_NOHTTPS=yes singularity build guidance_singularity.img singularity/guidance.def
