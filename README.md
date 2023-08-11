This repository contains code to reproduce the results in

> Bass AJ, Bian S, Wingo AP, Wingo TS, Culter DJ, Epstein MP. Identifying latent genetic interactions in genome-wide association studies using multiple traits. *Submitted*; 2023.

Overview of directories:

-   `simulations/`:
    -   `00-functions.R`: contains functions to simulate genotypes and traits.
    -   `01-lit-null.R`: type I error rate simulations.
    -   `02-lit-alt.R`: power simulations.
    -   `03-lit-time.R`: total computational time.
    -   `04-PC.R`: testing the principal components of the SQ/CP matrix.
    -   `05-figures.R`: code to create figures.
-   `analysis/`:
    -   `01-preprocessing.R`: preprocessing of the UK Biobank data.
    -   `02-run-plink.sh`: run PLINK to filter SNPs and individuals from the UK Biobank data. Note that we split the PLINK files into chunks using code from [here](https://bhoom.wordpress.com/2011/04/22/split-whole-genome-plink-binary-files-to-small-chunks/) to distribute across cores.
    -   `03-run-lit.R`: run LIT on the PLINK files generated from `02-run-plink.sh`.
    -   `04-run-marginal.R`: run Marginal (SQ/CP) and Marginal (SQ) on the PLINK files generated from `02-run-plink.sh`.
    -   `05-run-gamut.R`: run GAMuT on the PLINK files generated from `02-run-plink.sh`.
    -   `UKB_analysis.R`: code to analyze results.
