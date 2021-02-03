#!/bin/bash

for i in "$@"
do
case $i in
    -n=*|--run-name=*)
    RUN_NAME="${i#*=}"
    shift # past argument=value
    ;;
    -i=*|--input-dir=*)
    ROOT_INPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--output-dir=*)
    ROOT_OUTPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;
    --cpus=*)
    CPUS="${i#*=}"
    shift # past argument=value
    ;;
    --mem=*)
    MEM="${i#*=}"
    shift # past argument=value
    ;;
    --queue=*)
    QUEUE="${i#*=}"
    shift # past argument=value
    ;;
    --time=*)
    TIME="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

INPUT_BAM=${ROOT_INPUT_DIR}/${RUN_NAME}.bam

if [[ ! -f ${INPUT_BAM} ]];then
  echo "BAM FILE ${INPUT_BAM} NOT FOUND!"
  exit 1
fi

if [[ ! -d ${ROOT_OUTPUT_DIR}/jobs ]];then
  mkdir -p ${ROOT_OUTPUT_DIR}/jobs
fi

TMP_FOLDER=${ROOT_OUTPUT_DIR}/${RUN_NAME}_split
if [[ -d ${TMP_FOLDER} ]];then
  rm -r ${TMP_FOLDER}
fi

mkdir -p ${TMP_FOLDER}

JOB_NAME=scifi_pipeline.${RUN_NAME}.split_bam
JOB=${ROOT_OUTPUT_DIR}/jobs/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/jobs/${JOB_NAME}.%a.log

SAMPLES_DIR_SUFFIX="_samples"
if [[ $RUN_NAME =~ .*10xGenomics.* ]];then
  SAMPLES_DIR_SUFFIX=${SAMPLES_DIR_SUFFIX}"_10x"
fi

echo '#!/bin/env bash' > $JOB

echo "date" >> $JOB
echo '' >> $JOB

echo "#RUN_NAME               = ${RUN_NAME}" >> $JOB
echo "#INPUT_BAM              = ${INPUT_BAM}" >> $JOB
echo "#ROOT OUTPUT DIRECTORY  = ${ROOT_OUTPUT_DIR}" >> $JOB
echo "#TMP OUTPUT DIRECTORY   = ${TMP_FOLDER}" >> $JOB
echo "#SLURM PARAMETERS       = $CPUS, $MEM, $QUEUE, $TIME" >> $JOB
echo '' >> $JOB

echo "samtools split -@${CPUS} -f ${TMP_FOLDER}'/%!.%.' ${INPUT_BAM}" >> $JOB
echo '' >> $JOB

cat <<EOT >> $JOB
ls ${TMP_FOLDER}/*.bam|while read file;do
  filename=\$(basename \${file})
  lane=\$(echo \$filename|cut -d# -f1)
  run=\$(echo \$lane|cut -d_ -f1,2,3)
  samplesdir=${ROOT_OUTPUT_DIR}/\${run}/\${lane}${SAMPLES_DIR_SUFFIX}

  if [[ ! -d \${samplesdir} ]];then
    mkdir -p \${samplesdir}
  fi

  if [[ -f \${samplesdir}/\${filename} ]];then
    rm \${samplesdir}/\${filename}
  fi

  mv \${file} \${samplesdir}/
done

rm -r ${TMP_FOLDER}
EOT

echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
