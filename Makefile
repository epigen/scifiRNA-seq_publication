.DEFAULT_GOAL := all

# STAR CONFIGURATION
STAR_EXE := /home/dbarreca/bin/STAR

HUMAN_STAR_DIR :=  /data/groups/lab_bock/shared/resources/genomes/hg38/indexed_STAR-2.7.0e/
HUMAN_GTF_FILE := /data/groups/lab_bock/shared/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf

MOUSE_STAR_DIR :=  /home/dbarreca/resources/mm10_e99/indexed_STAR-2.7.0e/
MOUSE_GTF_FILE := /home/dbarreca/resources/mm10_e99/mm10_e99.gtf

HUMAN_MOUSE_STAR_DIR := /data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
HUMAN_MOUSE_GTF_FILE := /data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf

# OTHER RESOURCES
R2_BARCODES := $(shell pwd)/metadata/737K-cratac-v1.reverse_complement.csv

# INPUT/OUTPUT CONFIGURATION

## FOLDER CONTAINING DEMULTIPLEXED DATA (ONE BAM FILE PER WELL)
RAW_DATA_DIR := $(shell pwd)/data/raw/run

## FOLDER CONTAINING MERGED BAM FILES (i.e. DOWNLOADED FROM GEO)
ROOT_MERGE_BAM_DIR ?= $(shell pwd)/data/merge_bam

## FOLDER CONTAINING SPLIT BAM FILES
### This is to be used if using the demultiplexed files
#### ROOT_INPUT_DIR := $(RAW_DATA_DIR)
### This is to be used if using the downloaded and split files
ROOT_INPUT_DIR := $(shell pwd)/data/split_bam

ROOT_OUTPUT_DIR ?= $(shell pwd)/data/pipeline_out/$(RUN_NAME)
ROOT_REPORTS_DIR ?= $(shell pwd)/data/pipeline_out/reports

parse:
	@[ "${RUN_NAME}" ] || ( echo "'RUN_NAME' is not set"; exit 1 )
	@[ "${N_BARCODES}" ] || ( echo "'N_BARCODES' is not set"; exit 1 )
	@[ "${FLOWCELLS}" ] || ( echo "'FLOWCELLS' is not set"; exit 1 )
	@[ "${ANNOTATION}" ] || ( echo "'ANNOTATION' is not set"; exit 1 )


EXPECTED_CELL_NUMBER ?= 200000
MIN_UMI_OUTPUT ?= 3
VARIABLES ?= "plate_well"
ARRAY_SIZE ?= 24

STAR_DIR := $(HUMAN_STAR_DIR)
GTF_FILE := $(HUMAN_GTF_FILE)

IS_MOUSE ?= 0
ifeq ($(IS_MOUSE), 1)
  STAR_DIR := $(MOUSE_STAR_DIR)
  GTF_FILE := $(MOUSE_GTF_FILE)
endif

SPECIES_MIXING ?= 1
SPECIES_MIX_FLAG :=
ifeq ($(SPECIES_MIXING), 1)
	SPECIES_MIX_FLAG := --species-mixture
	STAR_DIR := $(HUMAN_MOUSE_STAR_DIR)
  GTF_FILE := $(HUMAN_MOUSE_GTF_FILE)
endif

merge_bam: parse
	$(info "scifi_pipeline: merge_bam")
	sh src/scifi_pipeline.merge_bam.sh \
	--run-name=$(RUN_NAME) \
	--flowcells="$(FLOWCELLS)" \
	--n-barcodes=$(N_BARCODES) \
	--annotation=$(ANNOTATION) \
	--cpus=8 \
	--mem=8000 \
	--queue=shortq \
	--time=08:00:00 \
	--output-dir=$(ROOT_MERGE_BAM_DIR) \
	--input-dir=$(RAW_DATA_DIR)
	$(info "scifi_pipeline: done")

split_bam: parse
	$(info "scifi_pipeline: split_bam")
	sh src/scifi_pipeline.split_bam.sh \
	--run-name=$(RUN_NAME) \
	--cpus=8 \
	--mem=8000 \
	--queue=mediumq \
	--time=2-00:00:00 \
	--output-dir=$(ROOT_INPUT_DIR) \
	--input-dir=$(ROOT_MERGE_BAM_DIR)
	$(info "scifi_pipeline: done")

map: parse
	$(info "scifi_pipeline: map")
	sh src/scifi_pipeline.map.sh \
	--run-name=$(RUN_NAME) \
	--flowcells="$(FLOWCELLS)" \
	--n-barcodes=$(N_BARCODES) \
	--annotation=$(ANNOTATION) \
	--cpus=4 \
	--mem=60000 \
	--queue=shortq \
	--time=08:00:00 \
	--array-size=$(ARRAY_SIZE) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--input-dir=$(ROOT_INPUT_DIR) \
	--star-exe=$(STAR_EXE) \
	--star-dir=$(STAR_DIR) \
	--gtf=$(GTF_FILE)
	$(info "scifi_pipeline: done")

filter: parse
	$(info "scifi_pipeline: filter")
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--annotation=$(ANNOTATION) \
	--expected-cell-number=$(EXPECTED_CELL_NUMBER) \
	--min-umi-output=${MIN_UMI_OUTPUT} \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=01:00:00 \
	--array-size=$(ARRAY_SIZE) \
	--r2-barcodes=${R2_BARCODES}
	$(info "scifi_pipeline: done")

join: parse
	$(info "scifi_pipeline: join")
	sh src/scifi_pipeline.join.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=12000 \
	--queue=shortq \
	--time=08:00:00
	$(info "scifi_pipeline: done")

report: parse
	$(info "scifi_pipeline: report")

	mkdir -p $(ROOT_REPORTS_DIR)/$(RUN_NAME)

	sbatch -J scifi_pipeline.report.$(RUN_NAME) \
	-o $(ROOT_REPORTS_DIR)/$(RUN_NAME)/scifi_pipeline.report.log \
	-p shortq --mem 250000 --cpus 4 --time 0-08:00:00 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics.csv.gz \
	$(ROOT_REPORTS_DIR)/$(RUN_NAME)/$(RUN_NAME). \
	--plotting-attributes $(VARIABLES) $(SPECIES_MIX_FLAG)"

	sbatch -J scifi_pipeline.report-exon.$(RUN_NAME) \
	-o $(ROOT_REPORTS_DIR)/$(RUN_NAME)/scifi_pipeline.report-exon.log \
	-p shortq --mem 250000 --cpus 4 --time 0-08:00:00 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).exon.metrics.csv.gz \
	$(ROOT_REPORTS_DIR)/$(RUN_NAME)/$(RUN_NAME).exon. \
	--plotting-attributes $(VARIABLES) $(SPECIES_MIX_FLAG)"
	$(info "scifi_pipeline: done")


all: map filter join  report

clean:
	find . -name "*bam" ! -name external -delete

.PHONY: map filter join  report clean
.SILENT: all
