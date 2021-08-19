# Regenie WDL

This is a workflow for regenie. 

- The workflow regenie.wdl can be found on Dockstore: https://dockstore.org/workflows/github.com/briansha/Regenie_WDL/regenie:master?tab=info
- The workflow regenie_alternate can be found on Dockstore:  https://dockstore.org/workflows/github.com/briansha/Regenie_WDL/regenie_alternate:master?tab=info

regenie: https://rgcgithub.github.io/regenie/

Step1 and Step2 use docker images with regenie compiled with Boost Iostream installed: https://github.com/rgcgithub/regenie/wiki/Using-docker

The docker image for the Plot task uses the provided Dockerfile.

## regenie.wdl
  - Step1 - Made with BGEN files in mind to use as the input genetic data file. - (BGEN version 1.2, 8-bit probabilities).
  - Step2 - Use separate .bed, .bim, .fam files for each chromosome. (If testing chr 1-22, there should be 22 separate files).
  - PLINK can be used to convert bed, bim, and fam files to BGEN files outside of this workflow.
  - PLINK can also be used to convert pgen, pvar, and psam files to BGEN files outside of this workflow.
   
## regenie_alternate.wdl
  - Step1 - Made with BGEN files in mind to use as the input genetic data file. - (BGEN version 1.2, 8-bit probabilities)
  - Step2 - Use a BGEN file containing all chromosomal data (chr 1-22, etc.)
  - PLINK can be used to convert bed, bim, and fam files to BGEN files outside of this workflow.
  - PLINK can also be used to convert pgen, pvar, and psam files to BGEN files outside of this workflow.
