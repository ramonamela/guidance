#!/bin/bash -e

mvn clean package
scp {guidance.jar,set_environment.sh,set_environment_original.sh,set_environment_separate.sh,set_environment_separate_constant.sh,set_environment_variable.sh} pr1ees03@mn1.bsc.es:/gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/

#ssh pr1ees03@mn1.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_sub1000.batch'

ssh pr1ees03@mn3.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_T1D_chr20_23_eagle_minimac.batch'

#ssh pr1ees03@mn2.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_T1D_chr20_23_eagle_impute.batch'

#ssh pr1ees03@mn2.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_T1D_chr20_23_shapeit_minimac.batch'

#ssh pr1ees03@mn2.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_T1D_chr20_23_eagle_minimac.batch'

#ssh pr1ees03@mn1.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_T1D_chr23.batch'
#ssh pr1ees03@mn1.bsc.es 'sbatch /gpfs/scratch/pr1ees00/pr1ees03/GUIDANCE/launch_guidance_T1D_test.batch'
