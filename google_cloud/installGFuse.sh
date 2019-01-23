#!/bin/bash

sudo apt-get update && \
sudo apt-get install -y --no-install-recommends software-properties-common && \
sudo add-apt-repository ppa:longsleep/golang-backports && \
sudo apt-get update && \
sudo apt-get -y --no-install-recommends install git golang-go fuse && \
export GO15VENDOREXPERIMENT=1 && \
export GOPATH="$HOME/go" && \
go get -u github.com/googlecloudplatform/gcsfuse && \
sudo mkdir -p /opt/userBin && \
sudo mv $HOME/go/bin/gcsfuse /opt/userBin/ && \
rm -rf ~/go && \
sudo ln -s /opt/userBin/gcsfuse /bin/gcsfuse
