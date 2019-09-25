#!/bin/bash

# Uninstall current java version
sudo apt-get update && sudo apt-get install -y --no-install-recommends apt-utils && \ #apt utils and vi to make the life easier
dpkg-query -W -f='${binary:Package}\n' | grep -E -e '^(ia32-)?(sun|oracle)-java' -e '^openjdk-' -e '^icedtea' -e '^(default|gcj)-j(re|dk)' -e '^gcj-(.*)-j(re|dk)' -e '^java-common' | xargs sudo apt-get -y remove && \ # Uninstall current java version
sudo apt-get install -y --no-install-recommends openjdk-8-jre openjdk-8-jdk && \ # Java 8
sudo apt-get -y --no-install-recommends install python && \ # Now python is a dependency
sudo apt-get -y --no-install-recommends install maven subversion && \ # Build dependencies
sudo apt-get -y --no-install-recommends install openjdk-8-jdk graphviz xdg-utils && \ # Runtime dependencies
sudo apt-get -y --no-install-recommends install libtool automake build-essential && \ # Bindings-common-dependencies
sudo apt-get -y --no-install-recommends install openssh-server openssh-client && \ # SSH
sudo apt-get -y --no-install-recommends install libxml2 gfortran libpapi-dev papi-tools && \ # Extrae dependencies
sudo apt-get -y --no-install-recommends install openmpi-bin openmpi-doc libopenmpi-dev uuid-runtime curl bc git && \ # Misc dependencies
rm -rf ~/2.4 && cd ~ && \
git clone https://github.com/bsc-wdc/compss.git 2.4 && \
pushd ~/2.4/ && ./submodules_get.sh && ./submodules_patch.sh && popd && \
cat ~/.bashrc | grep -v JAVA_HOME > ~/newbashrc && mv ~/newbashrc ~/.bashrc && \
echo "export JAVA_HOME=\"/usr/lib/jvm/java-8-openjdk-amd64/\"" >> ~/.bashrc && \
source ~/.bashrc && \
echo "**************************************************************************************************${JAVA_HOME}"
pushd ~/2.4/builders && sudo -E ./buildlocal -M -B -A && popd


