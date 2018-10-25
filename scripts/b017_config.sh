#!/usr/bin/env bash

# assign the dir names to variables
KINBRE_DIR=/panfs/pfs.local/scratch/sjmac/b961k922/kinbre

PROJECT="b017"
PROJECT_DIR=${KINBRE_DIR}/b017_macdonald && cd ${PROJECT_DIR}
REF_GENOME=${PROJECT_DIR}/refs/dmel_main_chr_r6.03_masked.fa

## organize the folders for the data
DATA_DIR=${PROJECT_DIR}/data

# parental lines files
PARENTAL_DATA_RAW=${DATA_DIR}/parental_lines/raw
PARENTAL_DATA_FILT=${DATA_DIR}/parental_lines/filtered

# data for the f18
F18_DATA_RAW=${DATA_DIR}/intercross_f18/raw
F18_DATA_FILT=${DATA_DIR}/intercross_f18/filtered

# data for the highly recombinant progeny
HR_DATA_RAW=${DATA_DIR}/intercross_hr/raw
HR_DATA_FILT=${DATA_DIR}/intercross_hr/filtered

# prepare the output directories
OUTPUTS_DIR=${PROJECT_DIR}/outputs
PARENTS_ALN_DIR=${PROJECT_DIR}/outputs/aln/parental_lines
F18_ALN_DIR=${PROJECT_DIR}/outputs/aln/intercross_f18
HR_ALN_DIR=${PROJECT_DIR}/outputs/aln/intercross_hr


SCRIPTS_DIR=${PROJECT_DIR}/scripts
INFO_DIR=${PROJECT_DIR}/info
REFS_DIR=${PROJECT_DIR}/refs
REPORTS_DIR=${PROJECT_DIR}/report_files

# create the dirs if they don't already exist
mkdir -p ${PROJECT_DIR} \
         ${DATA_DIR} \
         ${PARENTAL_DATA_RAW} \
         ${PARENTAL_DATA_FILT} \
         ${F18_DATA_RAW} \
         ${F18_DATA_FILT} \
         ${HR_DATA_RAW} \
         ${HR_DATA_FILT} \
         ${SCRIPTS_DIR} \
         ${OUTPUTS_DIR} \
         ${PARENTS_ALN_DIR} \
         ${F18_ALN_DIR} \
         ${HR_ALN_DIR} \
         ${INFO_DIR} \
         ${REFS_DIR} \
         ${REPORTS_DIR}
