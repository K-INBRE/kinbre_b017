#!/bin/bash
#MSUB -l nodes=1:ppn=8,mem=60g,walltime=02:00:00
#MSUB -o /home/b961k922/scratch/kinbre/b017_macdonald/report_files/b017_process_radtags_f18_$(JOBID)-${MOAB_JOBARRAYINDEX}.out
#MSUB -l walltime=6:00:00
#MSUB -q sixhour
#MSUB -j oe
#MSUB -N b017_process_radtags_f18
#MSUB -t [1-10]

echo "Sourcing the config file"
source /panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald/scripts/b017_config.sh

ARRAY_INFILE=${SCRIPTS_DIR}/b017_f18_plates_ARRAY_input.txt

plate=$(sed -n "$MOAB_JOBARRAYINDEX"p ${ARRAY_INFILE})
echo "Current file: ${plate}"

F18_BARCODES=${INFO_DIR}/f18_barcodes.txt

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

activate_conda_environment b017 stacks

echo "*** Separate the relevant barcodes ***"
grep ${plate} ${F18_BARCODES} | cut -f 2,3 > ${INFO_DIR}/${PROJECT}_${plate}_barcodes.txt

echo "*** Processing radtags for plate ${plate}... "
process_radtags -p ${F18_DATA_RAW}/Sample_Sam3852_F18_${plate}/ \
                -o ${F18_DATA_RAW}/Sample_Sam3852_F18_${plate} \
                -b ${INFO_DIR}/${PROJECT}_${plate}_barcodes.txt \
                -i 'gzfastq' -e mseI -r -c -q || { echo "process_radtags failed" ; exit 1; }
