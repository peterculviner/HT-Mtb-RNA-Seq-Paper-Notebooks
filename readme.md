# Overview
This repository includes Jupyter notebooks used to produce analysis outputs, figures and tables for the forthcoming paper "Evolution of Mycobacterium tuberculosis transcription regulation is associated with increased transmission and drug resistance".

Analysis folders and notebooks are ordered by number. The first folder, `01_isolate_genomics_data`, includes multiple notebooks used to prepare WGS and RNA-Seq data for downstream analyses.

To run analyses, additional datasets are requred. These are uploaded to Mendeley data:
> Culviner, Peter (2025), “HT Mtb RNA-Seq Paper Data Files”, Mendeley Data, V1, doi: 10.17632/h9dd6hsvc2.1

The Mendeley dataset includes all trees used in this paper (`datasets/trees`) and all variant calling and ancestral reconstruction files (`datasets/variants`).

In addition to libraries availible via `conda`, the analyses in this paper require other packages written for this paper:
- [mtbvartools](https://github.com/peterculviner/mtbvartools) - tools for storage and access of large bacterial variant call datasets produced in this study.
- [phyoverlap2](https://github.com/peterculviner/phyoverlap2) - tools re-implementating the [phyoverlap](https://github.com/nathan-d-hicks/phyOverlap) algorithm for very large trees.
- [mtbrnaseq](https://github.com/peterculviner/mtbrnaseq) - tools for demultiplexing RNA-Seq datasets using our 

# Accessing Variant Call Data
Of particular interest to the Tuberculosis research community is the summary data from the parsimony ancestral reconstruction. These are at `datasets/variants/global/241104_event_calls.vcb/event_calls.annotated.csv`. Please note that these `csv` files are very long and may not be openable in excel or desktop software. For examples of accessing these datasets with python see notebooks such as `08_whiB6_analyses.ipynb > Cumulative WhiB6 variants` or `09_espACD_analyses.ipynb > EspACD upstream region map`.
