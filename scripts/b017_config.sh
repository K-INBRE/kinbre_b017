#!/usr/bin/env bash

# assign the dir names to variables
KINBRE_DIR=/panfs/pfs.local/scratch/sjmac/b961k922/kinbre

PROJECT="b017"
PROJECT_DIR=${KINBRE_DIR}/b017_macdonald
cd ${PROJECT_DIR} || { echo "could not cd into ${PROJECT_DIR}"; exit 1; }
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

PARENTS_ALN_DIR=${OUTPUTS_DIR}/aln/parental_lines
PARENTS_SNPS_DIR=${OUTPUTS_DIR}/variants/parental_lines

F18_ALN_DIR=${OUTPUTS_DIR}/aln/intercross_f18
F18_SNPS_DIR=${OUTPUTS_DIR}/variants/intercross_f18

HR_ALN_DIR=${OUTPUTS_DIR}/aln/intercross_hr
HR_SNPS_DIR=${OUTPUTS_DIR}/variants/intercross_hr


REPORTS_DIR=${PROJECT_DIR}/report_files
ALN_REPORTS_DIR=${REPORTS_DIR}/aln
SNPS_REPORTS_DIR=${REPORTS_DIR}/var_calling

SCRIPTS_DIR=${PROJECT_DIR}/scripts
INFO_DIR=${PROJECT_DIR}/info
REFS_DIR=${PROJECT_DIR}/refs

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
         ${PARENTS_SNPS_DIR} \
         ${F18_ALN_DIR} \
         ${F18_SNPS_DIR} \
         ${HR_ALN_DIR} \
         ${HR_SNPS_DIR} \
         ${INFO_DIR} \
         ${REFS_DIR} \
         ${REPORTS_DIR} \
         ${ALN_REPORTS_DIR} \
         ${SNPS_REPORTS_DIR}
