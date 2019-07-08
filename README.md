Chemoproteomics Data Analysis 
========================================

Authors: Yuting Xu, Andy Liaw, Huijun Wang

Copyright© 2019 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.

Last updated: 07/03/2019


Requirements
----------------
* R (>= 3.4.0) 
* R packages: drc, truncnorm, magrittr, dplyr, knitr, ggplot2, gridExtra, reshape2, tidyr, limma


Install
----------------
- Option 1) Install from github
``` r
if('rCPDMS' %in% installed.packages()){remove.packages('rCPDMS')}
# install.packages("devtools")
devtools::install_github("Merck/rCPDMS")
```
- Option 2) Install from R-package source file rCPDMS_0.1.3.tar.gz in this repository
``` r
if('rCPDMS' %in% installed.packages()){remove.packages('rCPDMS')}
install.packages("rCPDMS_0.1.3.tar.gz", repos = NULL, type = "source")
```


Usage
----------------
Follow the scripts in demo/demo_CPDMS.R to set parameters for model fitting.
Example datasets are included in folder data_raw.



Workflows
----------------
* Pull down 
* Cellular Thermal Shift Assay (CETSA) 
* Iso-Thermal Dose Response (2dCETSA) 


Notes
----------------
The example code in miscellaneous/processMaxquantProtein.R demonstrates how to take a maxquant protein output (protein.txt) and the mapping file (MappingTemplaete.csv) to generate input data file for rCPDMS package.
