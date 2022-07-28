Claudia Del Vecchio, Bethan Cracknell Daniels, Giuseppina Brancaccio et al. Impact of antigen test target failure and testing strategies on the transmission of SARS-CoV-2 variants, 28 March 2022, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-1481444/v1]

# About this repository

This repository contains all the code needed to reproduce the results presented in Impact of antigen test target failure and testing strategies on the transmission of SARS-CoV-2 variants, in full. All data used is made available in the repository. 

All required packages and their versions are available in the DESCRIPTION file. Instructions to install download Rstan can be found [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). 



# Repository structure

This repository is divided into 3 sections: 

* *hospital_based_surveillance* - code to analyse antigen assay performance observed during hospital surveillance in Veneto.

* *model_fitting* - code to fit a multivariant model to reconstruct the transmission dynamics of the dominant variants carrying the A220V and M234I-A376T mutations, all other co-circulating variants, and the alpha variant (B.1.1.7) for Veneto and the rest of Italy. **Note that the code is currently set up to run each model both locally and on a HPC for 200 iterations with the first 100 discarded as burn in. This is to demo the models quickly (~4 hours on a normal computer). In order to obtain the results provided in the manuscript, 4 chains were run for 2000 iterations, discarding the first 1000 iterations as burnin. This takes ~6  hours per model to run on the HPC.** 

* *counterfactuals*- code to run models of the baseline and counterfactual testing scenarios using posterior samples obtained from *model_fitting* as input. Also includes code to run the sensitivity analyses, code to calculate test performance metrics and estimate variant detection by genomic surveillance, and code to produce Figures 3,5,6 and 7 from the manuscript. 


## hospital_based_surveillance

**This folder contains**

* *S2.xlsx* -  data on frequency of variants carying concordant and discordant mutations during the study period in Veneto and the rest of Italy 
* *S2.Rmd* - Rmarkdown file to run analysis antigen assay performance observed during hospital surveillance in Veneto. Also includes code to produce figures S3, S4, S5, S9 and S10, and tables 1, S1 and S2 from the manuscript.

## *model_fitting*

**This folder contains**

### Data

* *Dataset_Italy_A_v5.csv*, *Dataset_Italy_Alpha_v1.csv*, *Dataset_Italy_M_v5.csv*, *Dataset_Italy_O_v1.csv*: 
For Italy excluding Veneto, data on the number of sequences tested, the number of sequences testing positive by variant, obtained from the GISAID databank (1), data on the total monthly reported incidence and mean daily reported incidence in month M, obtained from the Civil Protection (2).  

* *Dataset_Veneto_A_v5.csv*, *Dataset_Veneto_Alpha_v1.csv*, *Dataset_Veneto_M_v5.csv*, *Dataset_Veneto_O_v1.csv*:
For Veneto, data on the number of sequences tested, the number of sequences testing positive by variant, obtained from the GISAID databank (1), data on the total monthly reported incidence and mean daily reported incidence in month M, obtained from the Civil Protection (2).  

* *Italy_daily_test_data.csv*, *Italy_monthly_test_data.csv*: 
For Italy excluding Veneto, data on the daily reported number of antigen and molecular tests and mean daily reported number of antigen and molecular tests in month M, obtained from the Civil Protection (2). 

* *Veneto_daily_test_data.csv*, *Veneto_monthly_test_data.csv*: 
For Veneto, data on the daily reported number of antigen and molecular tests and mean daily reported number of antigen and molecular tests in month M, obtained from the Civil Protection (2). 

* *dailyReportedIncidence_italy.csv*:
For Italy excluding Veneto, data on the daily reported incidence, obtained from the Civil Protection (2).  

* *dailyReportedIncidence_veneto.csv*:
For Veneto, data on the daily reported incidence, obtained from the Civil Protection (2).  

* *Vac_Italy_For_Month.csv* *data_vac_italy_day*
For Italy excluding Veneto, data on the mean daily per-capita number of second doses administered in month M / day i, obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

* *Vac_Veneto_For_Month.csv* *data_vac_veneto_day*
For Veneto, data on the mean daily per-capita number of second doses administered in month M / day i, obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

### Scripts 
* *run_model_fitting.R* is the main script to fit a multivariant model to epidemiological and genomic data from Veneto / the rest of Italy. We fit parameters using a Markov chain Monte Carlo (MCMC) framework using the No-U-Turn sampler via Rstan. Stan models can be run on HPC or locally. 

### Model 
* Rstan models to fit a multivariant model to epidemiological and genomic data from Veneto and the rest of Italy, assuming a negative binomial likelihood:

 - *est_test_symp.stan* Model variant which assumes that only symptomatic individuals undergo diagnostic testing and that they isolate if they receive a positive result (main model). 
 
 - *est_test_asymp.stan* Model variant which assumes that all symptomatic individuals isolate and that asymptomatic individuals test with a probability ρ and isolate if they receive a positive result.
 
 - *est_test_asymp2.stan* Model variant which assumes that symptomatic individuals test and isolate if they receive a positive result, and that asymptomatic individuals test with a probability ρ and isolate if they receive a positive result. 
 
 - *est_test_asymp_and_symp.stan* Model variant which assumes that the probability of taking a diagnostic test is independent of symptom occurrence, and that asymptomatic and symptomatic individuals test with the same probability ρ and isolate if they receive a positive result. 
 
  - *est_test_symp_alpha.stan* Model variation of  *est_test_symp.stan* which assumes that that symptomatic infectious individuals who did not test or who receive a false negative test result limit their contacts by a factor (1-α). 
  

### R
The R folder contains the functions used by the main script

 * *format_stan_data.R* takes as input raw data from *data* folder and parameter values and ouputs a list of data to input to Stan model for fitting. 
 * *sample_stan_model.R* functions to draw initial values for MCMC chains and run sampler. 
 * *diagnose_stan_fit.R* function to run diagnostics on Stan models and save likelihood. 
 * *plot_model_fit.R* function to plot Stan model fit against observed data. 
 * *save_post_chains.R* function to save posterior chains for all parameters. 
 * *run_model_fitting.R* function which sources all other functions and runs them. 


## *counterfactuals*

**This folder contains**

### Data 
The data folder contains all observed data - see above for description (model_fitting/data)

### MF_fitting results
 
This folder contains all results obtained by running *run_model_fitting.R* for each model variant (see folders: *symp_test*, *asymp_test*, *asymp2_test*, *asymp_symp_test*). 


### Scripts 

The order of scripts presented here indicates the order in which the analysis was conducted and is presented in the paper. Some scripts output results used in later scripts, however data and results are pre-provided so each script can also be run as a standalone. 

* *script_plot_modelling_data.R* Script to plot data on number of antigen and molecular tests, the reported incidece, the reported GISAID variant prevalence and the reconstructed variant reported incidence in Veneto and the rest of Italy (Figure 3). 

* *script_compare_main_models.R* Script to calcualte the DIC for each of the 4 model variants and summarise their posterior distributions. 

* *script_compare_main_models.R* Script to calcualte the DIC for each of the 4 model variants and summarise their posterior distributions (Table 2, Supplementary Table 6).
 
* *script_run_main_analysis.R* Script to run the baseline model of testing in Veneto and the rest of Italy (molecular tests follows positive antigen test) to reproduce observed data presented in the main analysis and the probability of variant specific case detection (i.e., outputs the data presented Figures 5a-c)

* *script_run_testing_cf.R* Script to run the counterfactual testing scenarios for Veneto (assuming proportion of antigen and molecular tests conducted in Veneto was the same as in the rest of Italy over the study period; assuming a molecular test followed a negative antigen test; assuming an antigen only testing strategy and an antigen test sensitivity of 68.9%; assuming an antigen only testing strategy and an antigen test sensitivity of 87.5%. Outputs the data presented in Figure 5d. 

* *script_run_transmission_cf.R* Script to run the counterfactual transmission scenarios for Veneto (i.e., assuming different scaling factors of R0M). Outputs the data presented in Figure 5e. 

* *script_plot_model_fit_and_cf.R* Script to plot the transmission dynamics of all variants under the baseline testing scenario, theprobability of variant specific case detection and cumulative incidence obtained under counterfactual scenarios. (Figure 5). 

* *script_run_and_plot_test_performance.R* Script to estimate and plot the performance of alternative antigen and molecular based testing strategies in diagnosing concordant and discordant variants (Figure 6).

* *script_run_and_plot_genom_surv.R* Script to estimate and plot the proportion of genomic detection of concordant and discordant variants under alternative antigen and molecular based testing strategies (Figure 7). 

* *script_posterior_sensitivity_analyses.R* Script to summarise the posterior distributions of each sensitivity analysis (Supplementary Table 7)


### Models 


Fixed-parameter Rstan models which take random samples of the posterior distribution obtained by running the multivariant Rstan model in *model_fitting* as input.

 - *est_test_symp_cf.stan* Fixed model variant which assumes that only symptomatic individuals undergo diagnostic testing and that they isolate if they receive a positive result (main model). 
 
 - *est_test_asymp_cf.stan* Fixed model variant which assumes that all symptomatic individuals isolate and that asymptomatic individuals test with a probability ρ and isolate if they receive a positive result.
 
 - *est_test_asymp2_cf.stan* Fixed model variant which assumes that symptomatic individuals test and isolate if they receive a positive result, and that asymptomatic individuals test with a probability ρ and isolate if they receive a positive result. 
 
 - *est_test_asymp_and_symp_cf.stan* Fixed model variant which assumes that the probability of taking a diagnostic test is independent of symptom occurrence, and that asymptomatic and symptomatic individuals test with the same probability ρ and isolate if they receive a positive result.
 
- *est_test_symp_cf_PCR_after_ag.stan* - counterfactual model of testing in Veneto where a molecular tests follows a negative antigen test. 
  
- *est_test_symp_cf_pPCR_0.stan* - counterfactual model of testing in Veneto where only antigen tests are used (i.e., no molecular testing). 


### R 

The R folder contains the function used by the scripts: 

 * *format_stan_data_cf.R* takes as input raw data from *data* folder and parameter values and ouputs a list of data to input for fixed Stan models. 

* *sample_posterior_chains.R* Function to sample from posterior chains in order to reconstruct variant specific incidences and run counterfactual scenarions; Function to sample from posterior chains in order to calculate the DIC of the main model variants; Function to sample from the posterior chains, calculate the R0 and format results to present in manuscript. Requires function *calculate_R0.R*. 

* *calculate_R0.R* Functions to calculate the R0 for all model varaints. 

* *calc_dic.R* Function to calculate the DIC. 

* *run_cf.R* Function to run fixed Stan models, taking as input parameter draws from the posterior distribution, obtained during model fitting; functions to extract and format the output of the running fixed Stan models; function to calculate the probability of detecting a case. 

* *plot_model_fit.R* function to plot the model transmission dynamics against the observed data presented in the main analysis (I.e., Figure 5a-b). 

* *plot_PPV_NPV_p_t+.R*  Function to calculate PPV/NPV, p(T+), dependent on test sensitivity and specificity, variant prevalence, and alternative antigen and molecular based testing strategies.  

* *plot_variant_detection_prob.R*  Function to calculate proportion of variant detected by genomic surveillance, dependent on % samples sequenced, test sensitivity, variant prevalence, and alternative antigen and molecular based testing strategies. 


### Results
 
This folder contains all results pertaining to the baseline model, model variants and counterfactual scenarios:

- *summarise_models* folder contains a summary of the posterior distributions for each of the 4 model variants and each of the 8 sensitivity analysis. 
- *baseline* folder contains the incidence and probability of detection for each variant, in each location for the main model, baseline testing scenario.
- *baseline* folder contains the incidence and probability of detection for each variant, in each location for the main model, counterfactual testing / transmission scenarios. 


### Figures 

This folder contains the main Figures (Figure 3, 5,6,7) generated by running the scripts contained within the counterfactuals folder. 



## References 
(1) Elbe, S. & Buckland‐Merrett, G. Data, disease, and diplomacy: GISAID’s innovative contribution to global health. Global challenges 1, 33-46 (2017).

(2) Civile, P.d.C.d.M.-D.d.P. Dati COVID-19 Italia. 2021. [cited 3 March 2022] Available from:https://github.com/pcm-dpc/COVID-19 

(3) Commissario straordinario per l'emergenza Covid-19 - Presidenza del Consiglio dei Ministri. Open Data su consegna e somministrazione dei vaccini anti COVID-19 in Italia. 2021 16th November 2021 [cited 3 March 2022] Available from: https://github.com/italia/covid19
