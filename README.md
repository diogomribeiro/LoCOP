# LoCOP: The molecular basis, genetic control and pleiotropic effects of local gene co-expression
This repository contains scripts for data processing, analysis and figure generation data for our paper:
**TODO: ADD LINK TO PAPER**

## Analysis scripts
- **expression_quantification**/FilterExpressionMatrix.py : used for preprocessing input gene expression matrices
- **cod_identification/CODer.py** : script to identify local co-expressed gene pairs (COPs) given a gene expression matrix
- **cod_analysis** : folder with multiple scripts used for analysis of COPs (including regression, parsing molecular features, variant pleiotropy, etc)
- **paper_figures** : contains all R scripts used to plot all figures present in the manuscript.

## Data availability
Data on co-expressed genes and shared eQTLs discovered here are available for consultation and download through the LoCOP DB **(LINK)** database developed here. 

## Geuvadis
EBI ArrayExpress (accession code E-GEUV-1) for RNA-seq data
1000 Genomes (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) for the genotype data. 
## GTEx 
RNA-seq and genotype data are available from dbGaP (accession: phs000424.v8.p2). 

[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## License
LoCOP is available under a MIT license. For more information please see the [LICENSE](LICENSE).
 

