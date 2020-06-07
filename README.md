# Guidance

GUIDANCE is an integrated framework that is able to perform haplotype phasing, genotype imputation, association testing assuming different models of inheritance and phenome-wide association analysis (PheWAS) of large GWAS datasets. Moreover, this application allows performing all these steps in a single execution, as well as in a modular way with optional user intervention.

The GUIDANCE strength is based on the possibility of increasing the potential of GWAS by means of the analysis of the X chromosome, as well as the autosomes analysis, using one or multiple reference panels, do the assosiation testing using several models (additive, dominant, recessive, heterodominant and genotypic inheritance models) and performing a cross-phenotype association analysis when more than one disease is available in the cohort under study.

# Table of Contents
1. [Getting started](#getstarted)
    * [Installing on a Singularity image](#singularity)
    * [Installing on bare metal](#bare)
2. [Configuration file description](#config)
    * [Input description](#input)
    * [Configuration file examples](#examples_conf)
    * [Environment file example](#environ_exam)
3. [GUIDANCE execution](#execution)
    * [Running on a Singularity image](#singularityexec)
    * [Running on bare metal](#bareexec)
4. [Run guidance test example](#test)
5. [Authors](#authors)
6. [License](#license)
7. [Acknowledgments](#acknowledgments)

## Getting Started <a name="getstarted"></a>

GUIDANCE was developed to run on HPC infrastructures and cloud environments. The instructions given here consider this execution in HPC environment. In order to run it in the cloud, specific instructions are given in this [project](https://github.com/ramonamela/guidance_cloud).

These instructions will get you a copy of the project up and running on your machine. There are two main options to run GUIDANCE: directly on top of your OS and inside a [Singularity image](https://github.com/sylabs/singularity). We strongly recommend to use Singularity as it makes all the process much more easier.

### Installing on a Singularity image <a name="singularity"></a>
<details><summary>Show instructions</summary>

<p>

#### Prerequisites

In order to run GUIDANCE in this mode, you only need to install the guidance [COMPSs’](https://github.com/bsc-wdc/compss/tree/guidance) branch and generate the Singularity image containing all the depedencies. If you want to generate execution traces, you should install this [COMPSs’](https://github.com/ramonamela/compss/tree/guidance) branch.

#### Installation

It is possible to find the instructions to install COMPSs in their github repository. 

Regarding the Singularity image generation, it is enough to run [this script](https://github.com/ramonamela/guidance/blob/master/singularity/build_container.sh). It will automatically create a local docker registry and upload to it the base docker image necessary to create the Singularity image. Nevertheless, you should modify the file placed in [$PROJECT_PATH/singularity/singularity/base.def](https://github.com/ramonamela/guidance/blob/master/singularity/singularity/base.def) in order to include the folder in your local file system that should be visible into the image. More precisely, you should only touch this lines:
```
%setup
    mkdir -p $SINGULARITY_ROOTFS/gpfs/home/
    mkdir -p $SINGULARITY_ROOTFS/gpfs/scratch/
    mkdir -p $SINGULARITY_ROOTFS/gpfs/apps/MN4
    mkdir -p $SINGULARITY_ROOTFS/gpfs/projects/
    mkdir -p /opt/intel
    mkdir -p /scratch
```

In addition, a GUIDANCE binary should be generated from this repository. In order to do so, you should run the following command in the root of the repository:

```
mvn clean install
```
This will generate a binary guidance.jar.

</p>
</details>

### Installing on bare metal <a name="bare"></a>
<details><summary>Show instructions</summary>

<p>

#### Prerequisites

First of all, you need to install the guidance [COMPSs’](https://github.com/ramonamela/compss/tree/guidance) branch in your system.
Afterwards, all the dependencies should be installed in the machine in which GUIDANCE will run. More precisely, the following binaries should be installed and accessible from all the worker nodes:

* Minimac4
* Minimac3
* Plink 1.9
* Eagle 2.4
* Impute 2.3.2
* Qctool 1.4
* Snptest 2.5
* Tabix
* Bgzip
* Shapeit v2 r727
* Samtools
* Bcftools
* R (tested version 3.5.1)

Next, the following environment variables should be set in order to tell GUIDANCE where the binaries are located in each machine:

```
export BCFTOOLSBINARY=/PATH_TO/BCFTOOLS/1.8/INTEL/bin/bcftools
export SAMTOOLSBINARY=/PATH_TO/SAMTOOLS/1.5-DNANEXUS/INTEL/IMPI/bin
export PLINKBINARY=/PATH_TO/plink_1.9/plink
export EAGLEBINARY=/PATH_TO/Eagle_v2.4.1/eagle
export RSCRIPTDIR=/PATH_TO/R_SCRIPTS/
export QCTOOLBINARY=/PATH_TO/qctool_v1.4-linux-x86_64/qctool
export SHAPEITBINARY=/PATH_TO/shapeit.v2.r727.linux.x64
export IMPUTE2BINARY=/PATH_TO/impute_v2.3.2_x86_64_static/impute2
export SNPTESTBINARY=/PATH_TO/snptest_v2.5
export MINIMAC3BINARY=/PATH_TO/Minimac3/bin/Minimac3
export MINIMAC4BINARY=/PATH_TO/Minimac4/release-build/minimac4
export TABIXBINARY=/PATH_TO/SAMTOOLS/1.5-DNANEXUS/INTEL/IMPI/bin/tabix
export BGZIPBINARY=/PATH_TO/SAMTOOLS/1.5-DNANEXUS/INTEL/IMPI/bin/bgzip
export RSCRIPTBINDIR=/PATH_TO/R/3.5.1/INTEL/bin/

export R_LIBS_USER=/gpfs/home/bsc05/bsc05997/TOOLS/R_libs/
```

In this case, R_LIBS_USER specify the path into which are installed the R dependencies in case they are not installed in the default path. More precisely, the following packages should be available:

* data.table
* plyr
* dplyr
* reshape
* gap
* sfsmisc
* BiocManager
* IRanges

In addition, in the folder “RSCRIPTDIR” there shuould be all the scripts available [here](https://github.com/ramonamela/guidance/tree/master/src/main/R).

Finally, maven should be installed in the system in order to compile the GUIDANCE binary.

#### Installing

Once all the dependencies are in the system, you should checkout the GUIDANCE repository and generate the guidance binary. In order to do so, you should run the following command in the root of the repository:

```
mvn clean install
```
This will generate a binary guidance.jar.

</p>
</details> 

## Configuration file <a name="config"></a>

In order to run GUIDANCE the user has to edit a configuration file, where the basic input and output characteristics have to be specified. This file also allows the tuning of multiple parameters related to, among others, covariates, chunk size for genotype imputation, info scores, minor allele frequency and Hardy-Weinberg thresholds. The user can also decide to run several phenotypes and several combinations of phenotypes/covariates in a single run. 

A detailed description of each of the configurable parameters (configuration_file_parameters.pdf) is provided in this [file](https://github.com/ramonamela/guidance/tree/master/webpage/configurationfileparameters.pdf) or in the collapsable list placed at the end of this paragraph.

<details><summary>Show the configuration file parameters</summary>

<p>

    • wfDeep: Name that defines the number of stages to be executed. These stages are defined in Figures 1 and 2.
    • init_chromosome: First chromosome to analyse.
    • end_chromosome: Last chromosome to analyse.
    • maf_threshold: Minor allele frequency cut-off used to filter final results.
    • impute_threshold: IMPUTE2 info score cut-off used to filter final results.
    • minimac_threshold: MINIMAC Estimated imputation accuracy (R²) cut-off used to filter final results.
    • hwe_cohort_threshold: Hardy-Weinberg equilibrium p.value threshold for cohort.
    • hwe_cases_threshold: Hardy-Weinberg equilibrium p.value threshold for cases.
    • hwe_controls_threshold: Hardy-Weinberg equilibrium p.value threshold for controls.
    • exclude_cgat_snps: Logical. Whether or not G>C or A>T SNPs should be excluded. We strongly recommend activating this flag as to avoid strand orientation issues. Most of the genotyping arrays have a very small number of such SNPs, and their exclusion should not result in any noticeable loss of imputation performance.
    • imputation_tool: The name of the imputation tool to impute genotypes. To date, only “impute” to select IMPUTE2 and “minimac” to select MINIMAC4 are accepted. 
    • test_types: Names for the different analysis to be carried out by GUIDANCE, separated by commas. The association results for each “test_type” will be created in a directory with the same name inside the “associations” directory. Below this flag, different “test_types” have to be listed with the phenotype name and the covariates names to take into account in the association analysis (for instance, to analyse “test_types = DIA2,CARD” users should add: “DIA2 = DIA2:sex,BMI” and “CARD = CARD:sex,BMI” below, where sex and BMI are covariates).
    • chunk_size_analysis: Size of the chunks considered to partition the data.
    • file_name_for_list_of_stages: File into which all the commands launched in the workflow are stored.
    • input_format: (I think that now we only support BED input since we have not tried with the other formats since I am working in the project…).
    • mixed_cohort: Name of the cohort.
    • mixed_bed_file_dir: The path to the directory with genotype files.
    • mixed_bed/bim/fam/_file: Name of the file containing genotypes. 
    • mixed_sample_file_dir: Path to the directory where the sample file is located.
    • mixed_sample_file:. Name of the sample file.
    • genmap_file_dir: Path where genetic map files are located.
    • genmap_file_chr_n: Name of the genetic map file for each chromosome in every new line.
    • refpanel_number: Number of reference panels.
    • refpanel_combine: 'NO' if there is only one panel or imputed results from different reference panels should not be integrated; 'YES' when different reference panels are expected to be used in the analysis and also the integration of all the results is required. 
    • refpanel_type: Name of the reference panel.
    • refpanel_memory: Required amount of memory demanded by each particular panel. Currently, “HIGH”, “MEDIUM” and “LOW” are supported.
    • refpanel_file_dir: Path where the reference panel for each chromosome is located.
    • refpanel_hap_file_chr_n: Haplotypes files per chromosome of the reference panel provided in case IMPUTE2 is chosen as imputation tool and for the chrX in case Minimac4 is used.
    • refpanel_leg_file_chr_n: Legend files per chromosome of the reference panel provided in case IMPUTE2 is chosen as imputation tool and for the chrX in case Minimac4 is used.
    • refpanel_vcf_file_chr_n: VCF files per chromosome of the reference panel provided in case Minimac4 is used.
    • outputdir: The path of the directory where the results will be written.

</p>
</details> 


### Input description <a name="input"></a>
GUIDANCE accepts PLINK format files, e.g. gwas.bed, gwas.bim, gwas.fam), as input. It also accepts covariates, such as principal compontents, gender or any covariate defined in the sample file. As part of the input, GUIDANCE also requires the inclusion of one, or several reference panels for imputation accepting public (1KG phase 1 or 3, HapMap, DCEG, UK10K, GoNL) or private reference panels. Also a genetic map is needed.

#### Sample File
The user needs to also provide an adequate sample file. Please, check SNPTESTv2 webpage for information on how to prepare a suitable sample file.

#### Genetic Map File
All of IMPUTE2 reference panel download packages come with appropriate recombination map files. Check IMPUTE2 webpage for more information. Hence, this genetic map files can be used for phasing with SHAPEIT. A file per chromosome must be given.

On the other hand, when phasing with EAGLE, a single file must be given. Several compatible options are available on the [Broad Website](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/).

#### Reference panels

It must be noted that, when using IMPUTE2, the reference panels must be given with format .haps and .leg. On the other hand, Minimac4 only accepts M3VCF. [Regular VCF can be converted to M3VCF using Minimac3](https://genome.sph.umich.edu/wiki/Minimac4). 

In addition, we have encountered some problems when imputing the ChrX with Minimac4. This is why the imputation of this chromosome is always done with IMPUTE2, so even when using Minimac4 both .haps and .legend must be given for ChrX if this chromosome is included in the study.

### Configuration file examples <a name="examples_conf"></a>

We give some templates of configuration files for a study from Chr21 to ChrX, using HRC, 1kphase3, uk10k and gonl as reference panels and [Eagle-IMPUTE2](https://github.com/ramonamela/guidance/tree/master/utils/conf_examples/config_eagle_impute.file), [SHAPEIT-IMPUTE2](https://github.com/ramonamela/guidance/tree/master/utils/conf_examples/config_eagle_minimac.file), [Eagle-Minimac4](https://github.com/ramonamela/guidance/tree/master/utils/conf_examples/config_shapeit_impute.file) and [SHAPEIT-Minimac4](https://github.com/ramonamela/guidance/tree/master/utils/conf_examples/config_shapeit_minimac.file) as phasing and imputation tools.

### Environment file example <a name="environ_exam"></a>

<details><summary>Show instructions</summary>

<p>
   
```
#!/bin/bash

### PHASE 1 ###

export phasingMem="50.0"
export phasingCU="48"
export phasingBedMem="50.0"
export phasingBedCU="48"


### PHASE 2 ###

export qctoolMem="16.0"
export qctoolSMem="1.0"
export gtoolsMem="6.0"
export samtoolsBgzipMem="6.0"
export imputeWithImputeLowMem="8.0"
export imputeWithImputeMediumMem="12.0"
export imputeWithImputeHighMem="20.0"
export imputeWithMinimacLowMem="4.0"
export imputeWithMinimacMediumMem="8.0"
export imputeWithMinimacHighMem="32.0"
export filterByInfoImputeMem="12.0"
export filterByInfoMinimacMem="24.0"

### PHASE 3 ###

export createListOfExcludedSnpsMem="1.0"
export filterHaplotypesMem="1.0"
export filterByAllMem="1.0"
export jointFilteredByAllFilesMem="15.0"
export jointCondensedFilesMem="1.0"
export generateTopHitsAllMem="2.0"
export generateTopHitsMem="2.0"
export filterByMafMem="2.0"
export snptestMem="2.0"
export mergeTwoChunksMem="1.0"
export mergeTwoChunksInTheFirstMem="1.0"
export combinePanelsMem="1.0"
export combineCondensedFilesMem="1.0"
export combinePanelsComplex1Mem="1.0"

### PHASE 4 ###

export generateCondensedTopHitsCU="48"
export generateCondensedTopHitsMem="90.0"
export generateQQManhattanPlotsCU="24"
export generateQQManhattanPlotsMem="45.0"
export phenoMergeMem="80.0"
```
This constraints correspond to all the phases executed during an execution. The most part of them, should be leaved like here. Nevertheless, some of them should be tuned dependeing on the execution:
* `phasingMem`: when setting this parameter, it should be taken into account that only one task per chromosome will be created. Hence, it should be set in such a way that all chromosomes can start being phased from the beggining but, at the same time, holding as many resources as possible. 
* `imputeWithImputeX`: this is the amount of memory used by IMPUTE when imputing the different chunks. This parameter will depend on the size of the used panel as well as the size of the input. Indeed, the greater the cohor, the greater the amount of memory needed.
* `imputeWithMinimacX`: this is the amount of memory used by Minimac when imputing the different chunks. This parameter will depend on the size of the used panel as well as the size of the input. Indeed, the greater the cohor, the greater the amount of memory needed.
* `generateX`: this corresponds to the end-files generation. As in the first step, should be set as high as possible as long as all the possible executions can run at once.
   
</p>
</details>

## GUIDANCE execution <a name="execution"></a>

### Running on a Singularity image <a name="singularityexec"></a>
<details><summary>Show instructions</summary>

<p>
   
An example of a quite complete launch script with singularity is shown next:
   
```
#!/bin/bash -e

export COMPSS_PYTHON_VERSION="2"

module purge
module load intel/2018.1
module load singularity/2.4.2

base_dir=$(pwd)
work_dir=${base_dir}/logs/
worker_work_dir=${base_dir}/tmpForCOMPSs/
worker_work_dir=scratch

source $base_dir/set_environment.sh

exec_time=2880
num_nodes=50
tracing=true
graph=true
debug=off
cpus_per_node=48
worker_in_master_cpus=0
worker_in_master_memory=80000
qos=bsc_cs

mkdir -p ${base_dir}/outputs_directory

/path/to/COMPSs/Runtime/scripts/user/enqueue_compss \
  --qos=${qos} \
  --job_dependency=5804941 \
  --graph=${graph} \
  --tracing=${tracing} \
  --log_level=${debug} \
  --exec_time=${exec_time} \
  --num_nodes=${num_nodes} \
  --base_log_dir=${base_dir} \
  --worker_in_master_cpus=${worker_in_master_cpus} \
  --worker_in_master_memory=${worker_in_master_memory} \
  --cpus_per_node=${cpus_per_node} \
  --master_working_dir=${work_dir} \
  --worker_working_dir=${worker_work_dir} \
  --scheduler="es.bsc.compss.scheduler.fifodatanew.FIFODataScheduler" \
  --classpath=/path/to/guidance.jar \
  --jvm_workers_opts="-Dcompss.worker.removeWD=true" \
  --container_image=/path/to/guidance_singularity.img \
  guidance.Guidance -config_file ${base_dir}/config_GERA_5000_shapeit_impute_1_23_cloud.file
```
In the next list, the most important features are explained in the same order as they appear in the script:
* Modules needed to run the execution.
* `work_dir`: the folder where the log files are stored.
* `worker_work_dir`: the folder where all the temporary files will be stored. If `tmp`, each worker node will use its `/tmp`. Otherwise, a shared directory between all the nodes should be specified (in general, this means pointing to an `nfs`, `gpfs` or `lustre` directory.
* File with all the environment variables pointing out the memory constraints necessary to run the execution.
* General constraints for the queue system (COMPSs will correctly traduce them to whichever installation present in the cluster).
* Creating the folder where the output files will be placed. Should be equal to the one stated in the configuration file.

Afterwards, in the launch command, there are 3 important files:
* `guidance_25_09_03_20_0_1_1.jar`: GUIDANCE binary.
* `guidance_singularity.img`: generated singularity image.
* `config_GERA_5000_shapeit_impute_1_23_cloud.file`: configuration file.

It is important to keep in mind that the output directory created should be equal to the one specified in the configuration file.

</p>
</details> 

### Running on bare metal <a name="bareexec"></a>
<details><summary>Show instructions</summary>

<p>
   
```
#!/bin/bash -e

module load COMPSs
module load mkl
module load intel/2017.4
module load samtools/1.5
module load R/3.5.1
module load bcftools/1.8
module load gcc/5.4.0

base_dir=$(pwd)
work_dir=${base_dir}/logs/
worker_work_dir=${base_dir}/tmpForCOMPSs/

export BCFTOOLSBINARY=/path/to/BCFTOOLS/1.8/INTEL/bin/bcftools
export RSCRIPTBINDIR=/path/to/R/3.5.1/INTEL/bin/
export SAMTOOLSBINARY=/path/to/SAMTOOLS/1.5-DNANEXUS/INTEL/IMPI/bin

export PLINKBINARY=/path/to/TOOLS/apps_gwimp_compss/plink_1.9/plink
export EAGLEBINARY=/path/to/TOOLS/Eagle_v2.4.1/eagle
export RSCRIPTDIR=/path/to/R_SCRIPTS/
export QCTOOLBINARY=/path/to/TOOLS/qctool_v1.4-linux-x86_64/qctool
export SHAPEITBINARY=/path/to/TOOLS/shapeit.v2.r727.linux.x64
export IMPUTE2BINARY=/path/to/TOOLS/impute_v2.3.2_x86_64_static/impute2
export SNPTESTBINARY=/path/to/TOOLS/snptest_v2.5
export MINIMAC3BINARY=/path/to/TOOLS/Minimac3/bin/Minimac3
export MINIMAC4BINARY=/path/toTOOLS//Minimac4/release-build/minimac4
export TABIXBINARY=/path/to/SAMTOOLS/1.5-DNANEXUS/INTEL/IMPI/bin/tabix
export BGZIPBINARY=/path/to/SAMTOOLS/1.5-DNANEXUS/INTEL/IMPI/bin/bgzip

export R_LIBS_USER=/path/to/TOOLS/R_libs/

export LC_ALL="C"

source $base_dir/set_environment.sh

exec_time=700
num_nodes=25
tracing=true
graph=true
log_level=off
cpus_per_node=48
worker_in_master_cpus=0
worker_in_master_memory=80000
qos=bsc_cs

mkdir -p ${base_dir}/outputs_shapeit_impute_1909_erase_all

enqueue_compss \
  --qos=${qos} \
  --job_dependency=7403259 \
  --graph=${graph} \
  --tracing=${tracing} \
  --log_level=${log_level} \
  --exec_time=${exec_time} \
  --num_nodes=${num_nodes} \
  --base_log_dir=${base_dir} \
  --worker_in_master_cpus=${worker_in_master_cpus} \
  --worker_in_master_memory=${worker_in_master_memory} \
  --cpus_per_node=${cpus_per_node} \
  --master_working_dir=${work_dir} \
  --worker_working_dir=${worker_work_dir} \
  --scheduler="es.bsc.compss.scheduler.fifodatanew.FIFODataScheduler" \
  --classpath=${base_dir}/guidance_25_1909_erase.jar \
  --jvm_workers_opts="-Dcompss.worker.removeWD=true" \
  guidance.Guidance -config_file ${base_dir}/config_GERA_300_shapeit_impute_1909_erase_all.file
```
In the next list, the most important features are explained in the same order as they appear in the script:
* Modules needed to run COMPSs.
* `work_dir`: the folder where the log files are stored.
* `worker_work_dir`: the folder where all the temporary files will be stored. If `tmp`, each worker node will use its `/tmp`. Otherwise, a shared directory between all the nodes should be specified (in general, this means pointing to an `nfs`, `gpfs` or `lustre` directory.
* Environment variables pointing to where all the needed binaries are placed.
* File with all the environment variables pointing out the memory constraints necessary to run the execution.
* General constraints for the queue system (COMPSs will correctly traduce them to whichever installation present in the cluster).
* Creating the folder where the output files will be placed. Should be equal to the one stated in the configuration file.

Afterwards, in the launch command, there are 3 important files:
* `guidance_25_09_03_20_0_1_1.jar`: GUIDANCE binary.
* `config_GERA_5000_shapeit_impute_1_23_cloud.file`: configuration file.

It is important to keep in mind that the output directory created should be equal to the one specified in the configuration file.

</p>
</details> 

## Run guidance test example <a name="test"></a>
We have set up a test with a really little dataset to both verify that an installation is correct and to get more familiar with all the configuration files and the execution process in a local machine. At this point, this only works in a Debian based system. 

In order to get this repository and run the example, only `git` and `make` are needed. It is possible to install them with the following instruction:
`sudo apt-get install -y make git`

Afterwards, you can clone this repository with the following command:
`git clone -b add_docker_compose https://github.com/ramonamela/guidance.git`

Finally, get into the root of the repository and perform this three steps:
1. Run `make setup-test-ubuntu`
This will install all the dependencies
2. Get out and login again if you are accessing the machine by `ssh` or open a new terminal
This step is necessary to refresh the available ubuntu group in order to correctly run docker
3. RUn `make run-execution`
This will run the execution with the sample input dataset. All the generated files will be in the `outputs`folder.

## Authors <a name="authors"></a>

* **Marta Guindo Martinez** - *R scripts, Java code and Scientific Methodology* - [ORCID](https://orcid.org/0000-0002-8099-9170)
* **Ramon Amela Milian** - *R scripts and Java code* - [ORCID](https://orcid.org/0000-0001-5943-5824) [github](https://github.com/ramonamela)
* **Cristian Ramón-Cortés Vilarrodona** - *Java code* - [ORCID](https://orcid.org/0000-0003-4170-818X) [github](https://github.com/cristianrcv)
* **Montserrat Puiggròs** - *R scripts and Java code* - [ORCID](https://orcid.org/0000-0001-5034-7924)
* **Josep Maria Mercader** - *Scientific Methodology* - [ORCID](https://orcid.org/0000-0001-8494-3660)
* **David Torrents** - *Scientific Methodology* - [ORCID](https://orcid.org/0000-0002-6086-9037)

<!-- See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project. -->

## License <a name="license"></a>

This project is licensed under the BSC Dual License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments <a name="acknowledgments"></a>

This work has been sponsored by the grant SEV-2011-00067 of Severo Ochoa Program, awarded by the Spanish Government. This work was supported by an EFSD/Lilly research fellowship. Josep M. Mercader was supported by Sara Borrell Fellowship from the Instituto Carlos III. Sílvia Bonàs was FI-DGR Fellowship from FI-DGR 2013 from Agència de Gestió d’Ajuts Universitaris i de Recerca (AGAUR, Generalitat de Catalunya). This study also makes use of data generated by the UK10K Consortium, derived from samples from UK10K COHORT IMPUTATION (EGAS00001000713). A full list of the investigators who contributed to the generation of the data is available in [uk10k](www.UK10K.org). Funding for UK10K was provided by the Wellcome Trust under award WT091310. The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of the Office of the Director of the National Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. The data used for the analyses described in this manuscript were obtained from the GTEx Portal on 07/16/2019. We acknowledge PRACE for awarding us access to both MareNostrum supercomputer from the Barcelona Supercomputing Center, based in Spain at Barcelona, and the SuperMUC supercomputer of the Leibniz Supercomputing Centre (LRZ), based in Garching at Germany (proposals numbers 2016143358 and 2016163985). The technical support group from the Barcelona Supercomputing Center is gratefully acknowledged. We also thank the Scientific Computing Services team from the Broad Institute for their technical support in installing and testing GUIDANCE in their machines. Finally, we thank all the Computational Genomics group at the BSC for their helpful discussions and valuable comments on the manuscript. We also acknowledge Elias Rodriguez Fos for designing the GUIDANCE logo. 
