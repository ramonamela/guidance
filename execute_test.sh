#!/bin/bash

FILE_DIR=/gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/
CLASSPATH=${FILE_DIR}/guidance.jar
$JAVA_HOME/bin/java -cp $CLASSPATH guidance.TestFunction "FILE1" "FILE2" "FILE3" "RESULTFILE" "0.5"
