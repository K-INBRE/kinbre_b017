#!/bin/bash
#MSUB -l nodes=1:ppn=8,mem=60g,walltime=02:00:00
#MSUB -o /home/b961k922/scratch/kinbre/b017_macdonald/report_files/b017_bcf_parents_$(JOBID).out
#MSUB -l walltime=6:00:00
#MSUB -q sixhour
#MSUB -j oe
#MSUB -N b017_bcf_parents

echo "Sourcing the config file"
source /panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald/scripts/b017_config.sh

while [[ ! -f "${PROJECT_DIR}/status/parents.align.done" ]]; do
  sleep 1m
done


parental_lines=( Db3852 Dsam )
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


activate_conda_environment ${PROJECT} bcftools

OUTPUT_PREFIX=${parental_lines[0]}_${parental_lines[1]}

bcftools mpileup -Ob -o ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.bcf \
                 --threads $(( ${N_THREADS} - 1 )) \
                 -f ${REF_GENOME} \
                 ${PARENTS_ALN_DIR}/${parental_lines[0]}.sorted.bam \
                 ${PARENTS_ALN_DIR}/${parental_lines[1]}.sorted.bam || { echo "bcftools mpileup failed" ; exit 1; }

bcftools call -vmO z --threads $(( ${N_THREADS} - 1 )) \
              -o ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.vcf.gz \
              ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.bcf || { echo "bcftools call failed" ; exit 1; }

tabix -p vcf ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.vcf.gz || { echo "tabix failed" ; exit 1; }

bcftools stats -F ${REF_GENOME} \
               -s - ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.vcf.gz > ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.vcf.gz.stats || { echo "bcftools stats failed" ; exit 1; }

mkdir ${PARENTS_SNPS_DIR}/plots || { echo "failed to created dir: ${PARENTS_SNPS_DIR}/plots" ; exit 1; }

plot-vcfstats -p ${PARENTS_SNPS_DIR}/plots ${PARENTS_SNPS_DIR}/${OUTPUT_PREFIX}.vcf.gz.stats || { echo "plot-vcfstats failed" ; exit 1; }
