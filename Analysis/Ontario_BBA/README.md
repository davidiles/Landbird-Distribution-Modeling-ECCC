# Analysis of 3rd Ontario Breeding Bird Atlas 

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




### Issues

Nov 26, 2024

to do:

- Create a "data visualization" script, prior to analysis to examine things like spatial/temporal coverage of surveys, changes in survey effort, etc
- store goodness-of-fit statistics produced by DHARMa package
   - this will also help to evaluate which error distributions are best for each species
- evaluate predictive accuracy using cross-validation
- store pixel-level predictions so that maps can be re-generated quickly
- incorporate time-varying covariates to represent habitat change explicitly
- incorporate checklists
- incorporate ARU effect, and ensure that ARU-based surveys are correctly labeled in NatureCounts dataframe
- confirm time zones / Hours-since-Sunrise are accurately calculated (are time zones correct?)
- simulations to confirm change analysis is working properly
