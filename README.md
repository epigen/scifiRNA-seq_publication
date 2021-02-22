scifi-RNA-seq
===================

Pipeline and analysis scripts for the scifi-RNA-seq publication.

> :warning: This repository includes source code to reproduce the analysis in the manuscript and won't be maintained to support the analysis of other datasets. For that, [see the general-purpose data processing pipeline of scifi-RNA-seq data](https://github.com/epigen/scifiRNA-seq).

The pipeline works by launching a single executable, the `scifi` python script.
This entrypoint script parses an annotation file passed as input and then
calls the [Makefile](Makefile) for each sample, which in turns calls the relevant shell script.

The annotations used for the publication are in the [metadata folder](metadata/). The main annotation file is [annotation.csv](metadata/annotation.csv)

All the environment variables (paths, reference, software, etc... ) are set at the top of the [Makefile](Makefile). Please adjust to your environment.

Please note that the various shell scripts in [src/scifi_pipeline.{step}.sh](src/) use SLURM to submit the step. For the map and filter steps an array job will be submitted. Please
adjust the queue names and slurm settings in the Makefile for each step.

The steps to reproduce the analysis must be ran one by one in the following order:

 - **split_bam** splits the merged bams downloaded from GEO into well-specific bam files and reconstructs the input folder structure from demultiplexing
 - **map** performs the mapping of the well-specific bam files to the reference transcriptome
 - **filter** selects the mapped reads in the mapped bam files and generates count matrices and metrics
 - **join** merges together the well-specific count matrices and metrics into sample-specific files
 - **report** generates various QC plots as well as scanpy's anndata objects serialized in h5ad format

To run the pipeline we suggest to use the following command: `./scifi {step} -t ./metadata/annotation.csv`. This will run the step specified only for those samples which have the toggle field set to 1 in the [annotation.csv](metadata/annotation.csv) file. We also suggest to setup a dedicate virtualenv to run the scripts, using the dependencies specified in the [requirements.txt](requirements.txt) file.

Other **requiremets** are:
  - STAR (and indexed genomes)
  - [featureCounts](http://subread.sourceforge.net/) binary installed on your PATH
  - [samtools](http://www.htslib.org/) binary installed on your PATH


Please note that 10x samples cannot be processed by the pipeline (except for the 'split' step). To process them in a consisten way, please use the [process_10x](src/method_comparisons/process_10x.sh) script to run those samples after running the 'split' step.

Check the [python executable](scifi), the [Makefile](Makefile) and [source files in src](src/) for more details.

Scripts used in the downstream analysis of the data are:
 - [monte_carlo_simulations.py](src/monte_carlo_simulations.py): for the theoretical "best-case scenario" simulation experiments;
 - [droplet_modeling.py](src/droplet_modeling.py): for the modeling of the Chromium device and prediction of collision rates;
 - [analysis.4lines_CROP-seq.py](src/analysis.4lines_CROP-seq.py): for the cell line mixture and CROP-seq experiments;
 - [analysis.PBMC_Tcell.py](src/analysis.PBMC_Tcell.py): for the experiments with primary human data;
 - [method_comparison.ipynb](src/method_comparison.ipynb): for the comparison across various methods (10X Chromium, sci-rna, sciPlex, SPLiT-seq, and scifi-RNA-seq).
