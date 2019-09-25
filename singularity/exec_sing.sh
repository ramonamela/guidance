#!/bin/bash

module load intel/2018.1
module load singularity
module load COMPSs
module remove python

pushd matmul
./launch_sing.sh None 2 10 48 false 4 128 1
popd
