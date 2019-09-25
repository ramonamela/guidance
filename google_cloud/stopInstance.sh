#!/bin/bash

source ./configure.sh

gcloud compute instances stop "${instanceName}"
