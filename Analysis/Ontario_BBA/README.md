# Ontario Breeding Bird Analysis 3

 This folder contains scripts that are in development, aimed at analysis of data from the [3rd Ontario Breeding Bird Analysis](https://www.birdsontario.org/), which began in 2021 and will conclude data collection in 2025.
 
# Methods
 
 *Note: Methods have not been formally peer-reviewed*
 
## Overview
  The analysis aims to describe spatial patterns of relative abundance and species occurrence across the province of Ontario over the time period from 2021-2025.
  
  The analysis will also evaluate change in abundance between the second Ontario Breeding Bird Atlas (2001-2005) and the third (2021-2025).
  
### Bird Data
Data were downloaded from two sources: 1) [WildTrax](https://wildtrax.ca/), and 2) [NatureCounts](https://naturecounts.ca/).

### Covariates
Covariate layers were obtained from the following sources:

- [2020 Land Cover of Canada](https://open.canada.ca/data/en/dataset/ee1580ab-a23d-4f86-a09b-79763677eb47)

### Model Structure

Models were fit using the integrated nested Laplace approximation (INLA), using R packages `INLA` and `inlabru`.  


