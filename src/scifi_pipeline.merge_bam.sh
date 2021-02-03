#! /bin/bash

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
    -f=*|--flowcells=*)
    FLOWCELLS="${i#*=}"
    shift # past argument=value
    ;;
    -l=*|--n-barcodes=*)
    N_BARCODES="${i#*=}"
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
    -a=*|--annotation=*)
    BARCODE_ANNOTATION="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

IFS=';' read -ra flowcells_arr <<< "$FLOWCELLS"

if [ $BARCODE_ANNOTATION = "nan" ]; then
  SAMPLES=$RUN_NAME
else
  SAMPLES=`tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 1`
fi

SAMPLES_DIR_SUFFIX="_samples"
if [[ $RUN_NAME =~ .*10xGenomics.* ]];then
  SAMPLES_DIR_SUFFIX=${SAMPLES_DIR_SUFFIX}"_10x"
fi

BAMS=()
for SAMPLE_NAME in $SAMPLES; do
  for FLOWCELL in "${flowcells_arr[@]}"
  do
    DIR=${ROOT_INPUT_DIR}/$(echo $FLOWCELL|cut -d_ -f1,2,3)

    if [[ $N_BARCODES -gt 1 ]]; then
        BAMS+=(`eval echo $DIR/${FLOWCELL}${SAMPLES_DIR_SUFFIX}/${FLOWCELL}#${SAMPLE_NAME}_{01..$N_BARCODES}.bam`)
    else
        SAMPLE_FILE=`echo $DIR/${FLOWCELL}${SAMPLES_DIR_SUFFIX}/${FLOWCELL}#${SAMPLE_NAME}.bam`
        if [ ! -f $SAMPLE_FILE ]; then
            SAMPLE_FILE=`echo $DIR/${FLOWCELL}${SAMPLES_DIR_SUFFIX}/${FLOWCELL}#${SAMPLE_NAME}_01.bam`
        fi

        BAMS+=(${SAMPLE_FILE})
    fi
  done
done

if [[ ! -d ${ROOT_OUTPUT_DIR}/jobs ]]; then
  mkdir -p ${ROOT_OUTPUT_DIR}/jobs
fi

JOB_NAME=scifi_pipeline.${RUN_NAME}.merge_bam
JOB=${ROOT_OUTPUT_DIR}/jobs/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/jobs/${JOB_NAME}.%a.log

echo '#!/bin/env bash' > $JOB

echo "date" >> $JOB
echo '' >> $JOB

echo "#RUN_NAME         = ${RUN_NAME}" >> $JOB
echo "#FLOWCELLS        = ${FLOWCELLS}" >> $JOB
echo "#BARCODE NUMBER   = ${N_BARCODES}" >> $JOB
echo "#ROOT INPUT DIRECTORY   = ${ROOT_INPUT_DIR}" >> $JOB
echo "#ROOT INPUT DIRECTORY   = ${ROOT_OUTPUT_DIR}" >> $JOB
echo "#SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME" >> $JOB
echo '' >> $JOB

echo "samtools merge -nrpf -@${CPUS} ${ROOT_OUTPUT_DIR}/${RUN_NAME}.bam ${BAMS[@]}" >> $JOB

echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
