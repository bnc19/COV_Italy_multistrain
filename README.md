# COV_Italy_multistrain

## Repository structure

This repository is divided into 2 sections: 

* *model_fitting* - code to fit a multivariant model to reconstruct the transmission dynamics of the dominant variants carrying the A220V and M234I-A376T mutations, all other co-circulating variants, and the alpha variant (B.1.1.7) for Veneto (between July 2020-May 2021) and the rest of Italy (between May 2020-2021), separately. 
* *counterfactuals*- code to run models of the baseline and counterfactual testing scenarios using posterior samples obtained from *model_fitting* as input. Also includes code to run the sensitivitiy analyses, code to calculate test performance metrics and estimate varaint detection by genomic surveillance, and code to produce figures from the manuscript. 

### model_fitting

**This folder contains**

#### Data

* *Dataset_Italy_A_v5.csv*, *Dataset_Italy_Alpha_v1.csv*, *Dataset_Italy_M_v5.csv*, *Dataset_Italy_O_v1.csv* : 
For Italy exluding Veneto, data on the number of sequences tested, the number of sequences testing positive by variant, obtained from the GISAID databank (1), data only the total monthly reported incidence and mean daily reported incidence in month M, obtained obtained from the Civil Protection (2).  

* *Dataset_Veneto_A_v5.csv*, *Dataset_Veneto_Alpha_v1.csv*, *Dataset_Veneto_M_v5.csv*, *Dataset_Veneto_O_v1.csv*:
For Veneto, data on the number of sequences tested, the number of sequences testing positive by variant, obtained from the GISAID databank (1), data only the total monthly reported incidence and mean daily reported incidence in month M, obtained obtained from the Civil Protection (2).  

* *Italy_daily_test_data.csv*, *Italy_monthly_test_data.csv* : 
For Italy exluding Veneto,  data on the daily reported number of antigen and molecular tests and mean daily reported number of antigen and molecular tests in month M, obtained from the Civil Protection (2). 

* *Veneto_daily_test_data.csv*, *Veneto_monthly_test_data.csv* : 
For Veneto,  data on the daily reported number of antigen and molecular tests and mean daily reported number of antigen and molecular tests in month M, obtained from the Civil Protection (2). 

* *dailyReportedIncidence_italy.csv* :
For Italy exluding Veneto, data on the daily reported incidence, obtained obtained from the Civil Protection (2).  

* *dailyReportedIncidence_veneto.csv* :
For Veneto, data on the daily reported incidence, obtained obtained from the Civil Protection (2).  

* *Vac_Italy_For_Month.csv*
For Italy exluding Veneto, data on the mean daily per-capita number of second doses administered in month M, obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

* *Vac_Veneto_For_Month.csv*
For Veneto, data on the mean daily per-capita number of second doses administered in month M, obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

#### Scripts 
* fit_SEIR_model_on_cluster.R is the main script to to fit a multivariant model to epidemiological and genomic data from Veneto / the rest of Italy. We fit parameters using a Markov chain Monte Carlo (MCMC) framework using the No-U-Turn sampler via Rstan. 

#### Model 
* *neg_bin_sens.stan* R stan model to fit a multivariant model to epidemiological and genomic data from Veneto / the rest of Italy, assuming a negative binomial likelihood. 

#### R
The R folder contains the function used by the main script
 * *run_testing_rstan_sens.R* 



### counterfactuals 

**This folder contains**

#### Italy 
The Italy folder contains data pertaining to Italy exluding the Veneto region. 

* *posteriorChains.csv*, *posteriorChains_sens.csv* :
Posterior distributions from fitting the multivariant model to epidemiological and genomic data from Italy excluding Veneto in Rstan, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis). 

* *100_posterior_samples_italy.csv*, *100_posterior_samples_italy_sens.csv* 100 random samples of the posterior distributions, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis), obtained using the script *script_sample_posterior.R*. 

* *Italy_variant_data.csv*  Mean and 95% CI daily variant specific incidence in month M, assumed to be the 1st of each month, derived from the available epidemiological and genomic data.

* *average_daily_vaccination_i_italy.csv* Mean daily vaccination in month M during Italy study period. 

* *daily_Ag_i_italy.csv* Daily antigen tests administered during Italy study period. 

* *daily_PCR_i_italy.csv* Daily molecular tests administered during Italy study period. 

* *x_i_data_italy.csv* integer data provided to rstan. Includes: Number of months to run the model, vector position of the first of each month, monthly average number of antigen and molecular tests administered, population assumed to be susceptible at the start of the study period (S comp), population assumed to have immunity at the start of the study period (R comp), number of days to run the model, time point interventions were implemented, time point to seed M234I-A376T and alpha variants.  

#### Veneto 
The Veneto folder contains data pertaining to the Italian region of Veneto. 

* *posteriorChains.csv*, *posteriorChains_sens.csv* :
Posterior distributions from fitting the multivariant model to epidemiological and genomic data from Veneto in Rstan, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis). 

* *100_posterior_samples_veneto.csv*, *100_posterior_samples_veneto_sens.csv* 100 random samples of the posterior distributions, assuming an antigen test sensitivity of 64.3% (main analysis) or 87.5% (sensitivity analysis), obtained using the script *script_sample_posterior.R*. 

* *Veneto_variant_data.csv*  Mean and 95% CI daily variant specific incidence in month M, assumed to be the 1st of each month, derived from the available epidemiological and genomic data.

* *average_daily_vaccination_i.csv* Mean daily vaccination in month M during Veneto study period. 

* *daily_Ag_i.csv* Daily antigen tests administered during Veneto study period. 

* *daily_PCR_i.csv* Daily molecular tests administered during Veneto study period. 

* *x_i_data.csv* integer data provided to rstan. Includes: Number of months to run the model, vector position of the first of each month, monthly average number of antigen and molecular tests administered in Veneto, population assumed to be susceptible at the start of the study period (S comp), population assumed to have immunity at the start of the study period (R comp), number of days to run the model, time point interventions were implemented, time point to seed M234I-A376T and alpha variants.  


Data for counterfactual analysis where we assume that the proportion of ANCOV and DNCOV tests conducted in Veneto was the same as in the rest of Italy over the study period: 

* *x_i_data_Veneto_italy_test.csv* integer data provided to rstan. Includes: Number of months to run the model, vector position of the first of each month, monthly average number of antigen and molecular tests administered in Italy, population assumed to be susceptible at the start of the study period (S comp), population assumed to have immunity at the start of the study period (R comp), number of days to run the model, time point interventions were implemented, time point to seed M234I-A376T and alpha variants.  

* * *daily_Ag_i_Itest.csv* Daily antigen tests administered in Italy during Veneto study period. 

* *daily_PCR_i_Itest.csv* Daily molecular tests administered in Italy during Veneto study period. 




### References 
(1) Elbe, S. & Buckland‐Merrett, G. Data, disease and diplomacy: GISAID’s innovative contribution to global health. Global challenges 1, 33-46 (2017).
(2) Civile, P.d.C.d.M.-D.d.P. Dati COVID-19 Italia. 2021.
(3) Commissario straordinario per l'emergenza Covid-19 - Presidenza del Consiglio dei Ministri. Open Data su consegna e somministrazione dei vaccini anti COVID-19 in Italia.  2021 16th November 2021 [cited  3 March 2022]Available from: https://github.com/italia/covid19-opendata-vaccini
