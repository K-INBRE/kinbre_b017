# b017: STITCH Pipeline

## 1. Create project directories

```
source /panfs/pfs.local/scratch/sjmac/b961k922/kinbre/b017_macdonald/scripts/b017_config.sh
```

The following data were moved from the _franklin_ RAID drive:

-   `${PARENTAL_DATA_RAW}` contains `xiaofei/KINBRE_Projects/SJM_data/Dmel_sam_b3852/WGS/fq/WGS_cat_sanger/*`  
-   `${F18_DATA_RAW}` contains `xiaofei/KINBRE_Projects/SJM_data/Dmel_sam_b3852/data2/fq/Project_F18_sam3852/Sample_*`  
-   `${HR_DATA_RAW}` contains `xiaofei/KINBRE_Projects/SJM_data/Dmel_sam_b3852/highly_Rec/Sample_*`  
-   `${INFO_DIR}` contains the following:
    -   `xiaofei/KINBRE_Projects/SJM_data/Dmel_sam_b3852/F18_samx3852_96_well_plates.xls`  
    -   `xiaofei/KINBRE_Projects/SJM_data/Dmel_sam_b3852/Sam3852_DNA_plate_layout.xls`

-   `${REFS_DIR}` contains `dmel_main_chr_r6.03_masked.fa`


Here is a helper function to run in terminal to simplify creating and activating conda environments.

```
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

```


## 2. Index the reference genome

```
activate_conda_environment b017 bwa

bwa index ${REFS_DIR}/dmel_main_chr_r6.03_masked.fa

source deactivate
```

## 3. QC and align Parental lines

```
# generate the input file for the ARRAY job
sh ${SCRIPTS_DIR}/generate_JOBARRAY_input.sh b017 ${PARENTAL_DATA_RAW}

# submit ARRAY job to align all parental samples
msub ${SCRIPTS_DIR}/b017_align_progeny_ARRAY.sh

```

## 4. Determine SNPs between the two parental lines

```
activate_conda_environment b017 bcftools

parental_lines=( Db3852 Dsam )

bcftools mpileup -Ou -f ${REF_GENOME} ${PARENTS_ALN_DIR}/${parental_lines[0]}.sorted.bam ${PARENTS_ALN_DIR}/${parental_lines[1]}.sorted.bam | bcftools call -vmO z -o ${PARENTS_ALN_DIR}/${parental_lines[0]}_${parental_lines[1]}.vcf.gz

bcftools mpileup -Ou -f ${REF_GENOME} ${PARENTS_ALN_DIR}/${parental_lines[0]}.sorted.bam ${PARENTS_ALN_DIR}/${parental_lines[1]}.sorted.bam | bcftools call -vmO z -o ${PARENTS_ALN_DIR}/${parental_lines[0]}_${parental_lines[1]}.vcf.gz
# tee ${REPORTS_DIR}/${parental_lines[0]}_${parental_lines[1]}.bcftools.log
source deactivate

```



```

cut -f 1 ${INFO_DIR}/f18_barcodes.txt | sort | uniq > ${SCRIPTS_DIR}/b017_f18_plates_ARRAY_input.txt
msub ${SCRIPTS_DIR}/b017_process_radtags_f18_ARRAY.sh

sed -i 's/p/Plate/g' ${INFO_DIR}/hr_barcodes.txt
cut -f 1 ${INFO_DIR}/hr_barcodes.txt | sort | uniq > ${SCRIPTS_DIR}/b017_hr_plates_ARRAY_input.txt

```
