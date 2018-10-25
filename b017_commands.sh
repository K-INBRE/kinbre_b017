#!/bin/bash

NTHREADS=8

PROJECTDIR=/panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald
RAWDATADIR=${PROJECTDIR}/data/parental_lines/raw
QCDIR=${PROJECTDIR}/data/parental_lines/filtered
REPORTSDIR=${PROJECTDIR}/report_files
REFGENOME=${PROJECTDIR}/refs/dmel_main_chr_r6.03_masked.fa
ALNDIR=${PROJECTDIR}/outputs/aln/parental_lines

mkdir -p QCDIR REPORTSDIR ALNDIR

parental_lines=( Db3852 Dsam )

activate_conda_environment () {
  ENVNAME=${1}_${2}
  # expected directory for the conda environment
  ENVDIR=/panfs/pfs.local/work/sjmac/software/conda_envs/${ENVNAME}
  # check if the directory exists; if it doesn't, create the env
  if [ ! -d "${ENVDIR}" ]; then
    # create a conda environment for the application
    # --yes = do not ask for confirmation during installation
    conda create --name ${ENVNAME} --yes
    # install the application in the newly created conda environment
    # --yes = do not ask for confirmation during installation
    conda install --name ${ENVNAME} ${2} --yes
  fi
  source activate ${ENVNAME}
}

#########################
# FASTP Quality Control #
#########################

cd ${PROJECTDIR}

activate_conda_environment b017 fastp

for parent in "${parental_lines[@]}"
do
  fastp -i ${RAWDATADIR}/${parent}_R1_sanger.fq.gz -I ${RAWDATADIR}/${parent}_R2_sanger.fq.gz \
        -o ${QCDIR}/${parent}_R1.filt.fq.gz -O ${QCDIR}/${parent}_R2.filt.fq.gz \
        -h ${REPORTSDIR}/${parent}.fastp.report.html \
        -j ${REPORTSDIR}/${parent}.fastp.report.json \
        --dont_overwrite \
        --cut_by_quality3 \
        --cut_window_size 5 \
        --cut_mean_quality 30 \
        --correction \
        --overrepresentation_analysis 2>&1 | tee ${REPORTSDIR}/${parent}.QC.fastp.log
done

source deactivate

##############################
# Index the reference genome #
##############################

activate_conda_environment b017 bwa

bwa index ${REFGENOME}

##########################
# Map the parental lines #
##########################

for parent in "${parental_lines[@]}"
do
  bwa sampe -t ${NTHREADS} -r "$(echo "@RG\tID:${parent}\tSM:${parent}\tLB:lib1")" \
            ${REFGENOME} \
            ${QCDIR}/${parent}_R1.filt.fq.gz ${QCDIR}/${parent}_R2.filt.fq.gz > ${ALNDIR}/${parent}.aln.sam
done

source deactivate

activate_conda_environment b017 samtools
for parent in "${parental_lines[@]}"
do
  samtools fixmate -O bam ${ALNDIR}/${parent}.aln.sam ${ALNDIR}/${parent}.aln.bam || exit 1
  samtools sort -t ${NTHREADS} -O bam -o ${ALNDIR}/${parent}.sorted.bam -T tmp_ ${ALNDIR}/${parent}.aln.bam || exit 1
  samtools index ${ALNDIR}/${parent}.sorted.bam
  rm ${ALNDIR}/${parent}.aln.sam
done

source deactivate
