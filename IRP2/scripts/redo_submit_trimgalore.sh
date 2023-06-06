for d in redo_*/;
do
  cd "$d"
  for R1 in *R1_001.fastq.gz
  do
  echo $R1
  R2=`echo ${R1} |sed 's/_R1_/_R2_/'`
  echo $R2
  if [ -f "${R1}" ] && [ -f "${R2}" ]; then
  echo "
  sbatch ../trimgalor.sh  ${R1} ${R2}"
  sbatch ../trimgalor.sh  ${R1} ${R2}
 else
  echo "ERROR. DATA FILES NOT FOUND."
  exit 4
 fi
  done
cd ..
done
