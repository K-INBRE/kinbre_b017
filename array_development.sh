for N in {1..6}
do
  LINE=$(sed -n "$N"p test.txt)
  echo ${LINE}

  myArray=(${LINE})
  echo "${myArray[@]}"
done



activate_conda_environment ${PROJECT} bcftools
bcftools mpileup -Ou -f ${REFGENOME} \
                     ${ALNDIR}/${parental_lines[0]}.sorted.bam \
                     ${ALNDIR}/${parental_lines[1]}.sorted.bam | \
                     bcftools call -vmO z \
                     -o ${ALNDIR}/${parent}.vcf.gz 2>&1 | tee ${REPORTSDIR}/parents.bcftools.log
source deactivate
