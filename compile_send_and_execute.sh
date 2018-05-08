#!/bin/bash -e

mvn clean package
scp guidance.jar pr1ees14@mn1.bsc.es:/gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE
#ssh pr1ees14@mn1.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/launch_guidance_sub1000.batch'
#ssh pr1ees14@mn1.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/launch_guidance_T1D_chr23.batch'
ssh pr1ees14@mn1.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/launch_guidance_T1D_test.batch'
