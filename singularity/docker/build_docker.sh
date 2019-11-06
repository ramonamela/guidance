#!/bin/bash

guidance_image_name=${1}

rm -rf ./TOOLS/R_scripts/
mkdir -p ./TOOLS/R_scripts/
cp -r ../../src/main/R/* ./TOOLS/R_scripts/

sudo docker build -f GuidanceDockerfile -t ${guidance_image_name} .

echo "[INFO] Docker build successfully executed."
#sudo docker run docker_guidance&
#sudo docker save --output=docker_singularity.tar docker_guidance
#sudo docker ps
#sudo sudo docker exec -i -t {name or id} /bin/bash

