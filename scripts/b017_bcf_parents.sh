#!/bin/bash
#MSUB -l nodes=1:ppn=8,mem=60g,walltime=02:00:00
#MSUB -o /home/b961k922/scratch/kinbre/b017_macdonald/report_files/b017_bcf_parents_$(JOBID).out
#MSUB -l walltime=6:00:00
#MSUB -q sixhour
#MSUB -j oe
#MSUB -N b017_bcf_parents

echo "Sourcing the config file"
source /panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald/scripts/b017_config.sh

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


activate_conda_environment ${PROJECT} bcftools

bcftools mpileup -Ou -f ${REFGENOME} \
                 ${ALNDIR}/${parental_lines[0]}.sorted.bam \
                 ${ALNDIR}/${parental_lines[1]}.sorted.bam | \
                 bcftools call -vmO z -o ${ALNDIR}/${parent}.vcf.gz || { echo "bcftools failed" ; exit 1; }
source deactivate
