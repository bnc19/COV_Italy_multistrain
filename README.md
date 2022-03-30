Claudia Del Vecchio, Bethan Cracknell Daniels, Giuseppina Brancaccio et al. Impact of antigen test target failure and testing strategies on the transmission of SARS-CoV-2 variants, 28 March 2022, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-1481444/v1]

# Repository structure

This repository is divided into 3 sections: 

* *hospital_based_surveillance* - code to analyse antigen assay performance observed during hospital surveillance in Veneto. Also includes code to produce figures S3, S4, S5, S9 and S10, and tables 1, S1 and S2 from the manuscript. 
* *model_fitting* - code to fit a multivariant model to reconstruct the transmission dynamics of the dominant variants carrying the A220V and M234I-A376T mutations, all other co-circulating variants, and the alpha variant (B.1.1.7) for Veneto (between July 2020-May 2021) and the rest of Italy (between May 2020-2021), separately. 
* *counterfactuals*- code to run models of the baseline and counterfactual testing scenarios using posterior samples obtained from *model_fitting* as input. Also includes code to run the sensitivity analyses, code to calculate test performance metrics and estimate variant detection by genomic surveillance, and code to produce figures 1,5,6,7S7 and S8 from the manuscript. 

* N.B. the scripts *fit_SEIR_model_on_cluster.R*, *script_run_main_analysis.R* and *script_run_and_plot_sensitivity_analysis.R* require Rstan installation to run. See instructions: https://mc-stan.org/users/interfaces/rstan

## hospital_based_surveillance

**This folder contains**

* *S2.xlsx* -  data on frequency of variants carying concordant and discordant mutations during the study period in Veneto and the rest of Italy 
* *S2.Rmd* - Rmarkdown file to run anal

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

* *Vac_Italy_For_Month.csv*
For Italy excluding Veneto, data on the mean daily per-capita number of second doses administered in month M, obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

* *Vac_Veneto_For_Month.csv*
For Veneto, data on the mean daily per-capita number of second doses administered in month M, obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

### Scripts 
* *fit_SEIR_model_on_cluster.R* is the main script to fit a multivariant model to epidemiological and genomic data from Veneto / the rest of Italy. We fit parameters using a Markov chain Monte Carlo (MCMC) framework using the No-U-Turn sampler via Rstan. 

### Model 
* *neg_bin_sens.stan* Rstan model to fit a multivariant model to epidemiological and genomic data from Veneto / the rest of Italy, assuming a negative binomial likelihood. 

### R
The R folder contains the function used by the main script
 * *run_testing_rstan_sens.R* 



## *counterfactuals*

**This folder contains**

### Italy 
The Italy folder contains data pertaining to Italy exluding the Veneto region. 

* *posteriorChains.csv*, *posteriorChains_sens.csv*:
Posterior distributions from fitting the multivariant model to epidemiological and genomic data from Italy excluding Veneto in Rstan, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis). 

* *100_posterior_samples_italy.csv*, *100_posterior_samples_italy_sens.csv* 100 random samples of the posterior distributions, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis), obtained using the script *script_sample_posterior.R*. 

* *Italy_variant_data.csv* Mean and 95% CI daily variant specific incidence in month M, assumed to be the 1st of each month, derived from the available epidemiological and genomic data.

* *average_daily_vaccination_i_italy.csv* Mean daily vaccination in month M during Italy study period. 

* *daily_Ag_i_italy.csv* Daily antigen tests administered during Italy study period. 

* *daily_PCR_i_italy.csv* Daily molecular tests administered during Italy study period. 

* *x_i_data_italy.csv* integer data provided to Rstan. Includes: Number of months to run the model, vector position of the first of each month, monthly average number of antigen and molecular tests administered, population assumed to be susceptible at the start of the study period (S comp), population assumed to have immunity at the start of the study period (R comp), number of days to run the model, time point interventions were implemented, time point to seed M234I-A376T and alpha variants.  

### Veneto 
The Veneto folder contains data pertaining to the Italian region of Veneto. 

* *posteriorChains.csv*, *posteriorChains_sens.csv*:
Posterior distributions from fitting the multivariant model to epidemiological and genomic data from Veneto in Rstan, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis). 

* *100_posterior_samples_veneto.csv*, *100_posterior_samples_veneto_sens.csv* 100 random samples of the posterior distributions, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis), obtained using the script *script_sample_posterior.R*. 

* *Veneto_variant_data.csv* Mean and 95% CI daily variant specific incidence in month M, assumed to be the 1st of each month, derived from the available epidemiological and genomic data.

* *average_daily_vaccination_i.csv* Mean daily vaccination in month M during Veneto study period. 

* *daily_Ag_i.csv* Daily antigen tests administered during Veneto study period. 

* *daily_PCR_i.csv* Daily molecular tests administered during Veneto study period. 

* *x_i_data.csv* integer data provided to Rstan. Includes: Number of months to run the model, vector position of the first of each month, monthly average number of antigen and molecular tests administered in Veneto, population assumed to be susceptible at the start of the study period (S comp), population assumed to have immunity at the start of the study period (R comp), number of days to run the model, time point interventions were implemented, time point to seed M234I-A376T and alpha variants.  

Data for counterfactual analysis where we assume that the proportion of ANCOV and DNCOV tests conducted in Veneto was the same as in the rest of Italy over the study period: 

* *x_i_data_Veneto_italy_test.csv* integer data provided to Rstan. Includes: Number of months to run the model, vector position of the first of each month, monthly average number of antigen and molecular tests administered in Italy, population assumed to be susceptible at the start of the study period (S comp), population assumed to have immunity at the start of the study period (R comp), number of days to run the model, time point interventions were implemented, time point to seed M234I-A376T and alpha variants.  

* *daily_Ag_i_Itest.csv* Daily antigen tests administered in Italy during Veneto study period. 

* *daily_PCR_i_Itest.csv* Daily molecular tests administered in Italy during Veneto study period. 


### Scripts 

The order of scripts presented here indicates the order in which the analysis was conducted and is presented in the paper. Some scripts output results used in later scripts, however data and results are pre-provided so each script can also be run as a standalone. 

* *script_plot_proportion_antigen_test.R* Script to plot data on number of antigen and molecular tests in Veneto and Italy (Figure 1)

* *script_format_cf_data.R* Script to format data to run counterfactuals. Takes external genomic and epidemiological data and user defined dates (e.g., start and end of the study period, variant seeding, and interventions) and returns data in the required format to run the Rstan models.  
 
* *script_sample_posterior.R* Script to sample from posterior distribution obtained by running the multivariant Rstan model in *model_fitting*. Also produces mean and 95% CrI from the 2.5 and 97.5 percentiles of the posterior distribution sample.
 
* *script_run_main_analysis.R* Script to run the baseline model of testing in Veneto and the rest of Italy (molecular tests follows positive antigen test) to reproduce observed data presented in the main analysis. Also runs the counterfactual testing scenarios for Veneto (assuming proportion of antigen and molecular tests conducted in Veneto was the same as in the rest of Italy over the study period; assuming a molecular test followed a negative antigen test; assuming an antigen-based mass testing strategy and an antigen test sensitivity of 64.3%; assuming an antigen-based mass testing strategy and an antigen test sensitivity of 87.5%; and assuming different scaling factors of R0M). 

* *script_plot_model_fit_and_cf.R* Script to plot the transmission dynamics of all variants under the baseline testing scenario, the % of incidence that is reported, cumulative incidence obtained under counterfactual scenarios (Figure 5). 

* *script_run_and_plot_test_performance.R* Script to estimate and plot the performance of alternative antigen and molecular based testing strategies in diagnosing concordant and discordant variants (Figure 6).

* *script_run_and_plot_genom_surv.R* Script to estimate and plot the proportion of genomic detection of concordant and discordant variants under alternative antigen and molecular based testing strategies (Figure 7). 

* *script_plot_R0.R* Script to plot R0 as a function of proportion of molecular tests conducted (Figure S7)

* *script_run_and_plot_sensitivity_analysis.R* Script to run and plot the results of the baseline model of testing in Veneto and the rest of Italy (molecular tests follows positive antigen test) to reproduce observed data for a sensitivity analysis where antigen test sensitivity is assumed to be 87.5% (Figure S8). 


### Models 

Fixed-parameter Rstan models which take random samples of the posterior distribution obtained by running the multivariant Rstan model in *model_fitting* as input.

* *stan_model.stan* - baseline model of testing in Veneto and the rest of Italy during the study period (molecular test follows a positive antigen test).

* *stan_model_PCR_after_ag.stan* - counterfactual model of testing in Veneto where a molecular tests follows a negative antigen test. 

* *stan_model_pPCR_0.stan* - counterfactual model of testing in Veneto where only antigen tests are used (i.e., no molecular testing). 

### R 
The R folder contains the function used by the scripts

* *Format_data_for_SEIQR_model.R* Function to format posterior chains and genomic and epidemiological data to run multivariant models assuming fixed parameters. Required to run *script_format_cf_data.R*. 

* *Format_variant_data.R* Function to format variant-specific genomic data to run multivariant models assuming fixed parameters. Required to run *script_format_cf_data.R*. 

* *replicate_rstan_fit.R* Functions to run multivariant models assuming fixed parameters. Required to run *script_run_main_analysis.R*.

* *plot_PPV_NPV_p_t+.R*  Function to calculate PPV/NPV, p(T+), dependent on test sensitivity and specificity, variant prevalence, and alternative antigen and molecular based testing strategies. Required to run *script_run_and_plot_test_performance.R*. 

* *plot_variant_detection_prob.R*  Function to calculate proportion of variant detected by genomic surveillance, dependent on % samples sequenced, test sensitivity, variant prevalence, and alternative antigen and molecular based testing strategies. Required to run *script_run_and_plot_genom_surv.R*. 

### Results 
Results from running *script_run_main_analysis.R*  and *script_sample_posterior.R* 


## References 
(1) Elbe, S. & Buckland‐Merrett, G. Data, disease, and diplomacy: GISAID’s innovative contribution to global health. Global challenges 1, 33-46 (2017).

(2) Civile, P.d.C.d.M.-D.d.P. Dati COVID-19 Italia. 2021. [cited 3 March 2022] Available from:https://github.com/pcm-dpc/COVID-19 

(3) Commissario straordinario per l'emergenza Covid-19 - Presidenza del Consiglio dei Ministri. Open Data su consegna e somministrazione dei vaccini anti COVID-19 in Italia. 2021 16th November 2021 [cited 3 March 2022] Available from: https://github.com/italia/covid19
