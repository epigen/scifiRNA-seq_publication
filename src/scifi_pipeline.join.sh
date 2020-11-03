#!/bin/bash

for i in "$@"
do
case $i in
    -n=*|--run-name=*)
    RUN_NAME="${i#*=}"
    shift # past argument=value
    ;;
    -v=*|--variables=*)
    VARIABLES="${i#*=}"
    shift # past argument=value
    ;;
    --species-mixture=*)
    SPECIES_MIXTURE="${i#*=}"
    shift # past argument=value
    ;;
    --output-dir=*)
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
    --correct-barcodes=*)
    CORRECT_BARCODES="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

echo "RUN_NAME         = ${RUN_NAME}"
echo "VARIABLES        = ${VARIABLES}"
echo "SPECIES_MIXTURE  = ${SPECIES_MIXTURE}"
echo "ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}"
echo "SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME"

# Start
mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR


if [[ $SPECIES_MIXTURE = "1" ]]; then
    ADDITIONAL_FIELDS="human,mouse,total,max,ratio,sp_ratio,doublet,unique_fraction,human_norm,total_norm,max_norm,ratio_norm,sp_ratio_norm,doublet_norm"
    ADDITIONAL_FIELDS="human,mouse,total,max,ratio,sp_ratio,doublet,read_human,read_mouse,read_total,read_max,read_ratio,read_sp_ratio,read_doublet,unique_fraction,human_norm,total_norm,max_norm,ratio_norm,sp_ratio_norm,doublet_norm,read_human_norm,read_total_norm,read_max_norm,read_ratio_norm,read_sp_ratio_norm,read_doublet_norm"
else
    ADDITIONAL_FIELDS="unique_fraction"
fi


OUTPUT_SUFFIX=""
if [[ ! -z "$CORRECT_BARCODES" ]]; then
    OUTPUT_SUFFIX="_corrected"
fi

ADDITIONAL_FIELDS="${ADDITIONAL_FIELDS},$VARIABLES"

# Simply concatenate all files
JOB_NAME=scifi_pipeline.${RUN_NAME}.concatenate_metrics
JOB=${ROOT_OUTPUT_DIR}/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/${JOB_NAME}.log
echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "find .  -mindepth 2 -name '*.metrics${OUTPUT_SUFFIX}.csv.gz' ! -name '*exon*' -exec cat {} \; > ${ROOT_OUTPUT_DIR}/${RUN_NAME}.metrics${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "echo 'r2,read,unique_umi,umi,gene,$ADDITIONAL_FIELDS' > ${RUN_NAME}_header1" >> $JOB
echo "gzip ${RUN_NAME}_header1" >> $JOB
echo "cat ${RUN_NAME}_header1.gz ${ROOT_OUTPUT_DIR}/${RUN_NAME}.metrics${OUTPUT_SUFFIX}.csv.gz > ${RUN_NAME}_tmp1" >> $JOB
echo "mv ${RUN_NAME}_tmp1 ${ROOT_OUTPUT_DIR}/${RUN_NAME}.metrics${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "rm ${RUN_NAME}_header1.gz" >> $JOB
echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB
sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB

JOB_NAME=scifi_pipeline.${RUN_NAME}.concatenate_expression
JOB=${ROOT_OUTPUT_DIR}/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/${JOB_NAME}.log
echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "find .  -mindepth 2 -name '*.expression${OUTPUT_SUFFIX}.csv.gz' ! -name '*exon*' -exec cat {} \; > ${ROOT_OUTPUT_DIR}/${RUN_NAME}.expression${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "echo 'r2,gene,umi,$VARIABLES' > ${RUN_NAME}_header2" >> $JOB
echo "gzip ${RUN_NAME}_header2" >> $JOB
echo "cat ${RUN_NAME}_header2.gz ${ROOT_OUTPUT_DIR}/${RUN_NAME}.expression${OUTPUT_SUFFIX}.csv.gz > ${RUN_NAME}_tmp2" >> $JOB
echo "mv ${RUN_NAME}_tmp2 ${ROOT_OUTPUT_DIR}/${RUN_NAME}.expression${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "rm ${RUN_NAME}_header2.gz" >> $JOB
echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB
sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB

JOB_NAME=scifi_pipeline.${RUN_NAME}.concatenate_metrics.exon
JOB=${ROOT_OUTPUT_DIR}/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/${JOB_NAME}.log
echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "find .  -mindepth 2 -name '*.exon.metrics${OUTPUT_SUFFIX}.csv.gz' -exec cat {} \; > ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.metrics${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "echo 'r2,read,unique_umi,umi,gene,$ADDITIONAL_FIELDS' > ${RUN_NAME}_header3" >> $JOB
echo "gzip ${RUN_NAME}_header3" >> $JOB
echo "cat ${RUN_NAME}_header3.gz ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.metrics${OUTPUT_SUFFIX}.csv.gz > ${RUN_NAME}_tmp3" >> $JOB
echo "mv ${RUN_NAME}_tmp3 ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.metrics${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "rm ${RUN_NAME}_header3.gz" >> $JOB
echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB
sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB

JOB_NAME=scifi_pipeline.${RUN_NAME}.concatenate_expression.exon
JOB=${ROOT_OUTPUT_DIR}/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/${JOB_NAME}.log
echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "find .  -mindepth 2 -name '*.exon.expression${OUTPUT_SUFFIX}.csv.gz' -exec cat {} \; > ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.expression${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "echo 'r2,gene,umi,$VARIABLES' > ${RUN_NAME}_header4" >> $JOB
echo "gzip ${RUN_NAME}_header4" >> $JOB
echo "cat ${RUN_NAME}_header4.gz ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.expression${OUTPUT_SUFFIX}.csv.gz > ${RUN_NAME}_tmp4" >> $JOB
echo "mv ${RUN_NAME}_tmp4 ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.expression${OUTPUT_SUFFIX}.csv.gz" >> $JOB
echo "rm ${RUN_NAME}_header4.gz" >> $JOB
echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB
sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
