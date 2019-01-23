#!/bin/bash

# Uninstall current java version
export DEBIAN_FRONTEND=noninteractive && \
sudo -E apt-get update && sudo apt-get install -y --no-install-recommends apt-utils && \
dpkg-query -W -f='${binary:Package}\n' | grep -E -e '^(ia32-)?(sun|oracle)-java' -e '^openjdk-' -e '^icedtea' -e '^(default|gcj)-j(re|dk)' -e '^gcj-(.*)-j(re|dk)' | xargs sudo apt-get -y remove && \
sudo -E apt-get -y --no-install-recommends install openjdk-8-jre openjdk-8-jdk && \
sudo -E apt-get -y --no-install-recommends install python && \
sudo -E apt-get -y --no-install-recommends install maven subversion && \
sudo -E apt-get -y --no-install-recommends install openjdk-8-jdk graphviz xdg-utils && \
sudo -E apt-get -y --no-install-recommends install libtool automake build-essential && \
sudo -E apt-get -y --no-install-recommends install openssh-server openssh-client && \
sudo -E apt-get -y --no-install-recommends install libxml2 libxml2-dev gfortran libpapi-dev papi-tools && \
sudo -E apt-get -y --no-install-recommends install openmpi-bin openmpi-doc libopenmpi-dev uuid-runtime curl bc git && \
sudo rm -rf ~/2.4 && cd ~ && \
git clone https://github.com/bsc-wdc/compss.git 2.4 && \
cd ~/2.4/ && ./submodules_get.sh && ./submodules_patch.sh && cd - && \
cat ~/.bashrc | grep -v JAVA_HOME > ~/newbashrc && mv ~/newbashrc ~/.bashrc && \
echo "export JAVA_HOME=\"/usr/lib/jvm/java-8-openjdk-amd64/\"" >> ~/.bashrc && \
cat ~/.bashrc && \
export JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/" && \
cd ~/2.4/builders && sudo -E ./buildlocal -M -A && cd -


