#!/bin/bash -e

main(){
  local snptest_binary=${1}
  local input_merged_gen_file=${2}
  local input_merged_sample_file=${3}
  local output_snptest_file=${4}
  local output_snptest_log_file=${5}
  local input_response_var=${6}
  local chr=${7}
  local input_covariables=${8}

  if [ "${chr}" -eq "23" ]; then
    if [ -z "${input_covariable}" ]; then
      ${snptest_binary} -data ${input_merged_gen_file} ${input_merged_sample_file} -o ${output_snptest_file} -pheno ${input_response_var} -hwe -log ${output_snptest_log_file} -method newml -assume_chromosome X -stratify_on sex -frequentist 1
    else
      ${snptest_binary} -data ${input_merged_gen_file} ${input_merged_sample_file} -o ${output_snptest_file} -pheno ${input_response_var} -cov_names ${input_covariables} -hwe -log ${output_snptest_log_file} -method newml -assume_chromosome X -stratify_on sex -frequentist 1
    fi
  else
    if [ -z "${input_covariable}" ]; then
      ${snptest_binary} -data ${input_merged_gen_file} ${input_merged_sample_file} -o ${output_snptest_file} -pheno ${input_response_var} -hwe -log ${output_snptest_log_file} -method em -frequentist 1 2 3 4 5
    else
      ${snptest_binary} -data ${input_merged_gen_file} ${input_merged_sample_file} -o ${output_snptest_file} -pheno ${input_response_var} -cov_names ${input_covariables} -hwe -log ${output_snptest_log_file} -method em -frequentist 1 2 3 4 5
    fi
  fi
}

main $@
