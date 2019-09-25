#!/bin/bash

for elem in $(ls /home/ramela/git/guidance/singularity/ | grep -v docker | grep -v scripts)
do 
    echo $elem
    scp -r $elem pr1ees14@mn1.bsc.es:~/singularity
done

