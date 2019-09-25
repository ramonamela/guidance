#!/bin/bash -e

mvn clean package
scp guidance.jar bsc05997@mn1.bsc.es:/gpfs/projects/bsc05/martagm/GWImp_COMPSs/GERA_sub200/
