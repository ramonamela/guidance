#!/bin/bash

#SBATCH --job-name=COMPSs_GERA_T2D
#SBATCH --exclusive
#SBATCH -t2:00:00
#SBATCH --workdir=/gpfs/scratch/pr1ees00/pr1ees14/GERA_T2D
#SBATCH -o compss-%J.out
#SBATCH -e compss-%J.err
#SBATCH -N1
#SBATCH -n48
#SBATCH --qos=debug

sample_file=$1
phasing_sample_file=$2
response_var=$3
covariables=$4

FILE_DIR=/gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/
CLASSPATH=${FILE_DIR}/guidance.jar
$JAVA_HOME/bin/java -cp $CLASSPATH guidance.TestFunction ${sample_file} ${phasing_sample_file} ${response_var} ${covariables}