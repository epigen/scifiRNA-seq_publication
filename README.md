scifi-RNA-seq publication
===================

The repository with source code used in the development of [scifi-RNA-seq (Datlinger et al.)](https://www.biorxiv.org/content/10.1101/2019.12.17.879304v1).

> :warning: This repository includes source code to reproduce the analysis in the manuscript and won't be maintained to support the analysis of other datasets. For that, [see the general-purpose data processing pipeline of scifi-RNA-seq data](https://github.com/epigen/scifiRNA-seq).

This repository contains [scripts used in the processing of data and its analysis](src/). For processing the data, the [Makefile](Makefile) runs discrete steps, but a full run can be done using the [submission script](scifi).

Metadata registering the experiments and their barcode annotation is [also avaialable](metadata/), and software required is listed in the [requirements file](requirements.txt).
