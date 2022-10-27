# HKU_Covid_Serosurvey_2022
This repository contains code and publicly available data for "Population-based sero-epidemiological estimates of real-world vaccine effectiveness against Omicron infection in an infection-naive population, Hong Kong, January to July 2022" 

Required software: MATLAB 2022a with Parallel Computing and Econometrics Toolboxes and R 4.2.1 with the tidyverse, pROC and cowplot packages

## Code files for analysis
Code files:
- main_estimate.m: Primary analysis for VE and waning
- logLikelihood.m: Log-likelihood function for MCMC inference
- loge.m: Log function floored at 1e-100
- logPriorDist.m: Flat prior for MCMC inference
- MCMC_Gibbs.m: Gibbs sampler

## Code files in "/IAR_calculation/" folders
(Note: separate folders have been created for different assumptions on delay to VE taking effect: 0, 7 and 14 days)
- main_p_estimate_vF_age.m: Primary analysis for IAR and population immunity
- infection_calc_vF.m: Function calculating probability of infection and immunity for each individual in Hong Kong

## Code files for Figures and Extended Data Figures
- Figure_1.R: Creates Figure 1
- /Figure_2/Figure_2_age.R: Creates Figure 2, Extended data Figure 5 and Extended data Figure 6 (#)
- Extended_data_figure_1.m: Creates Extended data Figure 1 (#)
- /controls/Extended_data_figure_3.R: Creates Extended Data Figure 3 (#)
- Extended_data_figure_4.R: Creates Extended Data Figure 4

## Notes on data availability
- Note 1: All source data have been included for files marked above with (#). 
- Note 2: Source data for files not marked above with (#) cannot be shared due to confidentiality undertakings to the Department of Health and the Environmental Protection Department, both of the Government of the HKSAR. Interested parties may contact these departments to request for access.
