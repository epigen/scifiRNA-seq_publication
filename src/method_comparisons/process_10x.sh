#!/usr/bin/env bash


function join_by { local IFS="$1"; shift; echo "$*"; }


#RUNS=(PD2XX1_10xscRNA_Human_Tcells_2S3Qmixed_Stimulated PD2XX1_10xscRNA_Human_Tcells_2S3Qmixed_Unstimulated PD2XX1_10xscRNA_Human_PBMCs_2S3Qmixed)
RUNS=(T46_10xGenomics_1_4lines_7650_nuclei T46_10xGenomics_2_4lines_7650_MeOH-cells T46_10xGenomics_3_4lines_7650_intact-cells)

for RUN_NAME in ${RUNS[@]}; do
FLOWCELL=BSF_0774_HNNGMDMXX
N_LANES=2
VARIABLES=plate_well
CELL_BARCODES=r1  # only round1 cell barcode has been tagged in the BAM file as "r1"
ROOT_OUTPUT_DIR=`pwd`/data/pipeline_out
STAR_EXE=/home/dbarreca/bin/STAR

# Human
STAR_DIR=/data/groups/lab_bock/shared/resources/genomes/hg38/indexed_STAR-2.7.0e/
GTF_FILE=/data/groups/lab_bock/shared/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf

# Mouse
# STAR_DIR=/home/dbarreca/resources/mm10_e99/indexed_STAR-2.7.0e/
# GTF_FILE=/home/dbarreca/resources/mm10_e99/mm10_e99.gtf

# Mix
# STAR_DIR=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
# GTF_FILE=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
# ADDITIONAL_ARGS="--species-mixture"

SUMMARIZER=`pwd`/src/scifi_pipeline.summarizer.py
CPUS=4
MEM=200G
QUEUE=mediumq
TIME=2-00:00:00

RAW_BASEDIR=`pwd`/data/raw/run

BAMS=(
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_1_samples_10x/${FLOWCELL}_1#${RUN_NAME}_01.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_1_samples_10x/${FLOWCELL}_1#${RUN_NAME}_02.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_1_samples_10x/${FLOWCELL}_1#${RUN_NAME}_03.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_1_samples_10x/${FLOWCELL}_1#${RUN_NAME}_04.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_2_samples_10x/${FLOWCELL}_2#${RUN_NAME}_01.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_2_samples_10x/${FLOWCELL}_2#${RUN_NAME}_02.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_2_samples_10x/${FLOWCELL}_2#${RUN_NAME}_03.bam
${RAW_BASEDIR}/${FLOWCELL}/${FLOWCELL}_2_samples_10x/${FLOWCELL}_2#${RUN_NAME}_04.bam
)
INPUT_BAM=`join_by , ${BAMS[@]}`
SAMPLE_NAME=${RUN_NAME}
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
PREFIX=${SAMPLE_DIR}/${SAMPLE_NAME}
mkdir -p $SAMPLE_DIR

JOB_DESCR=full
JOB_NAME=scifi_pipeline.${RUN_NAME}.${JOB_DESCR}
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

echo '#!/bin/env bash' > $JOB

echo "date" >> $JOB
echo '' >> $JOB

echo "#RUN_NAME         = ${RUN_NAME}" >> $JOB
echo "#FLOWCELL         = ${FLOWCELL}" >> $JOB
echo "#NUMBER OF LANES  = ${N_LANES}" >> $JOB
echo "#ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}" >> $JOB
echo "#STAR EXECUTABLE  = ${STAR_EXE}" >> $JOB
echo "#STAR DIRECTORY   = ${STAR_DIR}" >> $JOB
echo "#GTF FILE         = ${GTF_FILE}" >> $JOB
echo "#SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME" >> $JOB
echo '' >> $JOB

# align with STAR >=2.7.0e
echo '' >> $JOB
echo "$STAR_EXE \\
--runThreadN $CPUS \\
--genomeDir $STAR_DIR \\
--clip3pAdapterSeq AAAAAA \\
--outSAMprimaryFlag AllBestScore \\
--outSAMattributes All \\
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \\
--outSAMunmapped Within \\
--outSAMtype BAM Unsorted \\
--readFilesType SAM SE \\
--readFilesCommand samtools view -h \\
--outFileNamePrefix ${PREFIX}.STAR. \\
--readFilesIn $INPUT_BAM" >> $JOB

echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

# count all reads overlapping a gene
echo '' >> $JOB
echo "featureCounts \\
-T $CPUS \\
-F GTF \\
-t gene \\
-g gene_id \\
--extraAttributes gene_name \\
-Q 30 \\
-s 0 \\
-R BAM \\
-a $GTF_FILE \\
-o ${PREFIX}.STAR.featureCounts.quant_gene.tsv \\
${PREFIX}.STAR.Aligned.out.bam" >> $JOB

echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

# Same as above but just for exons
echo '' >> $JOB
echo "cd ${SAMPLE_DIR}" >> $JOB
echo "ln -s ${SAMPLE_NAME}.STAR.Aligned.out.bam \\
${SAMPLE_NAME}.STAR.Aligned.out.exon.bam" >> $JOB
echo 'cd -' >> $JOB
echo '' >> $JOB

echo "featureCounts \\
-T $CPUS \\
-F GTF \\
-t exon \\
-g gene_id \\
--extraAttributes gene_name \\
-Q 30 \\
-s 0 \\
-R BAM \\
-a $GTF_FILE \\
-o ${PREFIX}.STAR.featureCounts.quant_gene.exon.tsv \\
${PREFIX}.STAR.Aligned.out.exon.bam" >> $JOB


echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

# Filtering
echo "python3 -u $SUMMARIZER \\
--r1-barcode-as-r1-tag \\
--r1-attributes $VARIABLES \\
--cell-barcodes $CELL_BARCODES \\
--only-summary \\
--no-save-intermediate \\
--min-umi-output 3 \\
--expected-cell-number 10000 \\
--save-gene-expression \\
$ADDITIONAL_ARGS \\
--sample-name $SAMPLE_NAME \\
${SAMPLE_DIR}/${SAMPLE_NAME}.STAR.Aligned.out.bam.featureCounts.bam \\
${SAMPLE_DIR}/${SAMPLE_NAME}" >> $JOB

echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

echo "python3 -u $SUMMARIZER \\
--r1-barcode-as-r1-tag \\
--r1-attributes $VARIABLES \\
--cell-barcodes $CELL_BARCODES \\
--only-summary \\
--no-save-intermediate \\
--min-umi-output 3 \\
--expected-cell-number 10000 \\
--save-gene-expression \\
$ADDITIONAL_ARGS \\
--sample-name $SAMPLE_NAME \\
${SAMPLE_DIR}/${SAMPLE_NAME}.STAR.Aligned.out.exon.bam.featureCounts.bam \\
${SAMPLE_DIR}/${SAMPLE_NAME}.exon" >> $JOB

# Footer
echo '' >> $JOB
echo "date" >> $JOB
echo '' >> $JOB

# Submit
sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
done
