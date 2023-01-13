# HK Covid Serosurvey 2022
This repository contains code and available source data for "Real-world COVID-19 vaccine effectiveness against the Omicron BA.2 variant in a SARS-CoV-2 infection-naive population" 

Last revision: Jan 12, 2023

Required software: MATLAB 2022a with Parallel Computing and Econometrics Toolboxes and R 4.2.1 with the tidyverse, pROC, cowplot and janitor packages

## Code files for analysis
Code files:
- main_estimate.m: Primary analysis for VE and waning
- logLikelihood.m: Log-likelihood function for MCMC inference
- loge.m: Log function floored at 1e-100
- logPriorDist.m: Flat prior for MCMC inference
- MCMC_Gibbs.m: Gibbs sampler
- /Extended_Data_Fig_8/logLikelihood.m: Log-likelihood function for MCMC inference under seroconversion rate sensitvity scenarios as presented on Extended Data Figure 8 (all other code are the same as those used in the primary analysis)

## Code files in "/IAR_calculation/" folders
(Note: separate folders have been created for different assumptions on delay to VE taking effect: 7 (base case), 14 and 21 days)
- main_p_estimate_vF_age_[x]delay_weekly.m: Primary analysis for IAR and population immunity pursuant to assumptions on waning of infection-acquired immunity, where [x] denotes the delay to VE taking effect
- infection_calc_vF_weekly.m: Function calculating probability of infection and immunity for each individual in Hong Kong, incorporating different assumptions on waning of infection-acquired immunity

## Code files for Figures and Extended Data Figures
- Fig_1.R: Creates Figure 1
- Fig_2.R: Creates Figure 2
- /Figure_3_4_5_6_Extended_Fig_5_6_7/Figure_3_4_5_6_Extended_Figs_5_6_7.R: Creates Figures 3 through 6, and Extended Data Figures 5 through 7 (#)
- Extended_Data_Fig_1.m: Creates Extended Data Figure 1 (#)
- /controls/Extended_Data_Fig_3.R: Creates Extended Data Figure 3
- /Extended_Data_Fig_4/Extended_Data_Fig_4.R: Creates Extended Data Figure 4 (#)
- /Extended_Data_Fig_8/Extended_Data_Fig_8.R: Creates Extended Data Figure 8 (#)

## Notes on data availability
- Note 1: All source data have been included for files marked above with (#). 
- Note 2: Source data for files not marked above with (#) cannot be shared due to confidentiality undertakings to the Department of Health, the Office of the Government Chief Information officer, the Centre for Health Protection and the Environmental Protection Department, all of the Government of the HKSAR. Interested parties may contact these departments to request for access.
