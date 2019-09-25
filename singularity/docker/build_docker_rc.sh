#!/bin/bash

cp ../../src/main/R/* ./TOOLS/R_scripts/

sudo docker build -f GuidanceDockerfileRc -t docker_guidance_rc .

echo "[GUIDANCE] Docker build successfully executed."
echo "[GUIDANCE] Begin image generation."
#sudo docker run docker_guidance&
#sudo docker save --output=docker_singularity.tar docker_guidance
#sudo docker ps
#sudo sudo docker exec -i -t {name or id} /bin/bash

