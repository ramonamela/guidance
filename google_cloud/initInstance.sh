#!/bin/bash

source ./configure.sh

gcloud compute instances start "${instanceName}" --zone=${bucketLocation}
