#!/bin/bash -e

mvn clean package
scp {guidance.jar,set_environment.sh,set_environment_original.sh,set_environment_separate.sh,set_environment_separate_constant.sh,set_environment_variable.sh} pr1ees14@mn1.bsc.es:/gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE
