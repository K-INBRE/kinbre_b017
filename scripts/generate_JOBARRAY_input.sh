#!/bin/bash

DIR=$(readlink -f ${2})

for file in ${DIR}/*_R1*.f*q*
do   
readlink -f  ${file} >> ${1}_sample_file_paths_ARRAY.txt 
done


