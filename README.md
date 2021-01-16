scifi-RNA-seq publication
===================

The repository with source code used in the development of [scifi-RNA-seq (Datlinger et al.)](https://www.biorxiv.org/content/10.1101/2019.12.17.879304v1).

> :warning: This repository includes source code to reproduce the analysis in the manuscript and won't be maintained to support the analysis of other datasets. For that, [see the general-purpose data processing pipeline of scifi-RNA-seq data](https://github.com/epigen/scifiRNA-seq).

This repository contains [scripts used in the processing of data and its analysis](src/). For processing the data, the [Makefile](Makefile) runs discrete steps, but a full run can be done using the [submission script](scifi).

Metadata registering the experiments and their barcode annotation is [also avaialable](metadata/), and software required is listed in the [requirements file](requirements.txt).


Scripts used in the downstream analysis of the data are:
 - [monte_carlo_simulations.py](src/monte_carlo_simulations.py): for the theoretical "best-case scenario" simulation experiments;
 - [droplet_modeling.py](src/droplet_modeling.py): for the modeling of the Chromium device and prediction of collision rates;
 - [analysis.4lines_CROP-seq.py](src/analysis.4lines_CROP-seq.py): for the cell line mixture and CROP-seq experiments;
 - [analysis.PBMC_Tcell.py](src/analysis.PBMC_Tcell.py): for the experiments with primary human data;
 - [method_comparison.ipynb](src/method_comparison.ipynb): for the comparison across various methods (10X Chromium, sci-rna, sciPlex, SPLiT-seq, and scifi-RNA-seq).
