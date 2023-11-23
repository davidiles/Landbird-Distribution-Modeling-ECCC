# Saskatchewan Breeding Bird Analysis

 This folder contains the current analysis of relative abundance for the [Saskatchewan Breeding Bird Analysis](https://sk.birdatlas.ca/), which conducted data collection from 2017-2021.
 
# Methods
 
 *Note: Methods have not been formally peer-reviewed, but have been heavily tested and partially reviewed by David Iles and colleagues in ECCC and Birds Canada.*
 
## Overview
  The analysis aims to describe spatial patterns of relative abundance and species occurrence across the province of Saskatchewan over the time period from 2017-2021.
  
### Bird Data
Data were downloaded from two sources: 1) [WildTrax](https://wildtrax.ca/), and 2) [NatureCounts](https://naturecounts.ca/).

### Covariates
Covariate layers were obtained from the following sources:

- [2020 Land Cover of Canada](https://open.canada.ca/data/en/dataset/ee1580ab-a23d-4f86-a09b-79763677eb47)
- Digital elevation map (provided by Birds Canada)
- Stand Canopy Closure (extracted from [Beaudoin et al. 2014](https://cdnsciencepub.com/doi/10.1139/cjfr-2013-0401))


### Model Structure

Models were fit using the integrated nested Laplace approximation (INLA), using R packages `INLA` and `inlabru`.  


