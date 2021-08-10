# LoCOP: The molecular basis, genetic control and pleiotropic effects of local gene co-expression
This repository contains scripts for data processing, analysis and figure generation data for our paper:
Ribeiro D. M. <em>et al.</em> (2021) Nature Communications https://www.nature.com/articles/s41467-021-25129-x

## Analysis scripts
- **CODer.py** : Script to identify local co-expressed gene pairs (COPs) given a gene expression matrix. Example usage: python3 CODer.py expression_matrix.bed output_folder 1000 --fdrCutoff 0.01
- **FilterExpressionMatrix.py** : Script to preprocess input gene expression matrices. Example usage: python3 expression_matrix.bed expression_matrix_filtered.bed ---minimumSamples 0.5 --minimumValue 0.0
- **manuscript_figures** : Folder with all scripts used to produce figures for the paper. This includes the code to perform logistic regression analysis.

## Data availability
Data on co-expressed genes and shared eQTLs discovered here are available for consultation and download through the [LoCOP DB](http://glcoex.unil.ch) database developed here. 

## Geuvadis
EBI ArrayExpress (accession code E-GEUV-1) for RNA-seq data
[1000 Genomes](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) for the genotype data. 
## GTEx 
RNA-seq and genotype data are available from dbGaP (accession: [phs000424.v8.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2)). 

## License
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
LoCOP is available under a MIT license. For more information please see the [LICENSE](LICENSE).
