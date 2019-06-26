.DEFAULT_GOAL := all

# Change these parameters
RUN_NAME       := SCI023-TCells-S2
# RUN_NAME can be any name you want, but the comment cannot be at any location you desire
FLOWCELL       := BSF_0624_HJY5CDMXX
N_LANES        := 2
ANNOTATION     := /home/pdatlinger/projects/scifiRNA-seq/scifi-RNA-seq_analysis1_25-06-2019/sciRNA-seq.SCI024.oligos_2019-05-17_four-entries.csv
ROOT_OUTPUT_DIR:= /home/pdatlinger/projects/scifiRNA-seq/scifi-RNA-seq_analysis1_25-06-2019/$(RUN_NAME)
STAR_EXE       := /home/arendeiro/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR
STAR_DIR       := /home/arendeiro/resources/genomes/hg38/indexed_STAR-2.7.0e/
# for mixing experiments use: /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
GTF_FILE       := /home/arendeiro/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
# for mixing experiments use: /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
PYTHON_SUMMARIZER_SCRIPT     := /home/pdatlinger/projects/scifiRNA-seq/src/scifi_pipeline.summarizer.py
# this script is used in the 'filter' step
PYTHON_REPORT_SCRIPT         := /home/pdatlinger/projects/scifiRNA-seq/src/scifi_pipeline.report.py
# this script is used in the 'plot' step
ROUND2_WHITELIST   := /home/pdatlinger/projects/scifiRNA-seq/scifi-RNA-seq_analysis1_25-06-2019/737K-cratac-v1.reverse_complement.csv
# whitelist for 10x Genomics Single Cell ATAC barcodes, obtained from 10x Genomics website

# Don't change below unless you need different resources
map:
	echo "scifi_pipeline: map"
	sh src/scifi_pipeline.map.sh \
	--run-name=$(RUN_NAME) \
	--flowcell=$(FLOWCELL) \
	--n-lanes=$(N_LANES) \
	--annotation=$(ANNOTATION) \
	--cpus=4 \
	--mem=50000 \
	--queue=shortq \
	--time=08:00:00 \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--star-exe=$(STAR_EXE) \
	--star-dir=$(STAR_DIR) \
	--gtf=$(GTF_FILE)

filter:
	echo "scifi_pipeline: filter"
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--barcode_annotation=$(ANNOTATION) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=08:00:00 \
	--python-summarizer-script=$(PYTHON_SUMMARIZER_SCRIPT) \
	--round2_whitelist=$(ROUND2_WHITELIST)

join:
	echo "scifi_pipeline: join"
	sh src/scifi_pipeline.join.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--cpus=1 \
	--mem=12000 \
	--queue=shortq \
	--time=08:00:00

plot:
	echo "scifi_pipeline: plot"
	python3 -u $(PYTHON_REPORT_SCRIPT) \
	--metric-file ${ROOT_OUTPUT_DIR}/${RUN_NAME}.metrics.csv.gz \
	--output-prefix ${ROOT_OUTPUT_DIR}/plots/${RUN_NAME}. \
	--plotting-attributes donor_id \
	--expected-cell-number 250000

	python3 -u $(PYTHON_REPORT_SCRIPT) \
	--metric-file ${ROOT_OUTPUT_DIR}/${RUN_NAME}.exon.metrics.csv.gz \
	--output-prefix ${ROOT_OUTPUT_DIR}/plots/${RUN_NAME}.exon. \
	--plotting-attributes donor_id \
	--expected-cell-number 250000

all: map filter join plot

clean:
	find . -name "*bam" -delete

.PHONY: map filter join plot clean
.SILENT: all
