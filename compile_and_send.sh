#!/bin/bash -e

mvn clean package
scp guidance.jar pr1ees14@mn1.bsc.es:/gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/
