#!/bin/bash
#MSUB -l nodes=1:ppn=8,mem=60g,walltime=02:00:00
#MSUB -o /home/b961k922/scratch/kinbre/b017_macdonald/report_files/aln/b017_aln_$(JOBID)-${MOAB_JOBARRAYINDEX}.out
#MSUB -l walltime=6:00:00
#MSUB -q sixhour
#MSUB -j oe
#MSUB -N b017_aln_parents
#MSUB -t [1-2]

echo "Sourcing the config file"
source /panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald/scripts/b017_config.sh

ARRAY_INFILE=/panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald/b017_sample_file_paths_ARRAY.txt

SEQ_FILE_PATH=$(sed -n "$MOAB_JOBARRAYINDEX"p ${ARRAY_INFILE})
echo "Current file: ${SEQ_FILE_PATH}"

parent=$(basename ${SEQ_FILE_PATH} _R1_sanger.fq.gz)
echo "File prefix: ${parent}"

N_THREADS=8

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


echo "project ID: ${PROJECT}"

cd ${PROJECT_DIR}



#########################
# FASTP Quality Control #
#########################

# check if the filtered FASTQ file exists already
# if yes, skip the alignment
if [[ ! -f "${PARENTAL_DATA_FILT}/${parent}_R1.filt.fq.gz" ]]; then
	echo "Starting QC for sample: ${parent}"

  activate_conda_environment ${PROJECT} fastp

  fastp -i ${PARENTAL_DATA_RAW}/${parent}_R1_sanger.fq.gz -I ${PARENTAL_DATA_RAW}/${parent}_R2_sanger.fq.gz \
        -o ${PARENTAL_DATA_FILT}/${parent}_R1.filt.fq.gz -O ${PARENTAL_DATA_FILT}/${parent}_R2.filt.fq.gz \
        -h ${REPORTS_DIR}/${parent}.fastp.report.html \
        -j ${REPORTS_DIR}/${parent}.fastp.report.json \
        -w ${N_THREADS} \
        --dont_overwrite \
        --cut_by_quality3 \
        --cut_window_size 5 \
        --cut_mean_quality 30 \
        --correction \
        --overrepresentation_analysis || { echo "fastp failed" ; exit 1; }
  source deactivate
fi



##########################
# Map the parental lines #
##########################

# check if the alignment file (SAM) already exists
# if yes, skip the alignment
if [[ ! -f "${PARENTS_ALN_DIR}/${parent}.aln.sam" ]]; then
  echo "Starting alignment"
  activate_conda_environment ${PROJECT} bwa

  bwa mem -t ${N_THREADS} ${REF_GENOME} \
          -R "@RG\tID:${parent}\tSM:${parent}\tLB:lib1" \
          -o ${PARENTS_ALN_DIR}/${parent}.aln.sam \
          ${PARENTAL_DATA_FILT}/${parent}_R1.filt.fq.gz \
          ${PARENTAL_DATA_FILT}/${parent}_R2.filt.fq.gz
  source deactivate
fi

# check if the sorted BAM file already exists
# if yes, skip the conversion
if [[ ! -f "${PARENTS_ALN_DIR}/${parent}.sorted.bam" ]]; then
  echo "Starting SAM conversation"
  activate_conda_environment ${PROJECT} samtools
  samtools fixmate -O bam ${PARENTS_ALN_DIR}/${parent}.aln.sam \
                   -@ $(( ${N_THREADS} - 1 )) \
                   ${PARENTS_ALN_DIR}/${parent}.aln.bam || { echo "samtools fixmate" failed ; exit 1; }

  samtools sort -t ${N_THREADS} \
                -O bam -o ${PARENTS_ALN_DIR}/${parent}.sorted.bam \
                -T ${MOAB_JOBARRAYINDEX}_tmp_ ${PARENTS_ALN_DIR}/${parent}.aln.bam || { echo "samtools sort failed" ; exit 1; }

  samtools index ${PARENTS_ALN_DIR}/${parent}.sorted.bam ||  { echo "samtools index failed" ; exit 1; }
  source deactivate
fi
