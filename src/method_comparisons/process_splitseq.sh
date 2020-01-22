#!/usr/bin/env bash


function join_by { local IFS="$1"; shift; echo "$*"; }


RUNS=(
splitseq_300
splitseq_3000
)

for RUN_NAME in ${RUNS[@]}; do
FLOWCELL=GSE110823
N_LANES=1
VARIABLES=plate_well
ADDITIONAL_ARGS="--species-mixture"
CELL_BARCODES="r1 r2 r3"
ROOT_OUTPUT_DIR=data
STAR_EXE=/home/arendeiro/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR
STAR_DIR=/home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
GTF_FILE=/home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
SUMMARIZER=`pwd`/src/scifi_pipeline.summarizer.py
CPUS=4
MEM=60000
QUEUE=shortq
TIME=08:00:00

BAMS=(
/scratch/lab_bock/shared/projects/sci-rna/data/external/${FLOWCELL}/${RUN_NAME}.annotated.bam
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