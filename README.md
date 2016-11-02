# DGSA

This package implements distance based generalized sensitivity analysis method (DGSA). For more information about the method please read:  

Fenwick, Darryl, Céline Scheidt, and Jef Caers. "Quantifying asymmetric parameter interactions in sensitivity analysis: 
Application to reservoir modeling." Mathematical Geosciences 46.4 (2014): 493-511.

To install the package you can follow two paths:

1. Clone to your own computer, open Rstudio or R and type: install.packages("package_Location", repos = NULL, type=source). You also have 
   an option to load your clone as a project in Rstudio and then proceed with project build Ctrl+Shift+B (cmd+shift+B). This is particularly 
   useful if you plan on adding/modifying features.

2. Install directly from github. 
   First you will need to make sure that you have devtools installed (library(devtools)), then you can install with the following command:
   devtools::install_github("ogru/DGSA"). After that you can proceed to load the routines in the usual way with library(DGSA).

If you prefer to work with MATLAB, you can find the original MATLAB implementation here: https://github.com/scheidtc/dGSA. However, without 
the matrix plot for visualization of all sensitivities/importances. 

