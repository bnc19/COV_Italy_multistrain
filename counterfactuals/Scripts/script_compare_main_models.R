# Script to calculate DIC and return parameter estimates for all models --------


# Set up -----------------------------------------------------------------------

# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/calc_dic.R")
source("R/sample_posterior_chains.R")


# Incidence models -------------------------------------------------------------

# symp inc ---------------------------------------------------------------------

# inc DIC / WAIC
calc_dic(
  posterior_chains = read.csv("MF_results/symp_test/1/posterior_chains.csv"),
  log_lik = read.csv("MF_results/symp_test/1/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_symp_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689)


# M1 - (baseline) -   1125.51


# symp summary 
sum_symp_inc = summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/1/posterior_chains.csv"))

write.csv(sum_symp_inc, "MF_results/comp_main_model/1/sum_symp_inc.csv")

# asymp inc --------------------------------------------------------------------

# DIC  
calc_dic(
  posterior_chains = read.csv("MF_results/asymp_test/1/posterior_chains.csv"),
  log_lik = read.csv("MF_results/asymp_test/1/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
)


# M1   1174.05


# asymp summary 
sum_asymp_inc = summarise_posterior_chains_not_symp(
  posterior_chains = read.csv("MF_results/asymp_test/1/posterior_chains.csv"),
  asymp = T)

write.csv(sum_asymp_inc, "MF_results/comp_main_model/1/sum_asymp_inc.csv")


# asymp 2 inc -------------------------------------------------------------------

# DIC  
 calc_dic(
  posterior_chains = read.csv("MF_results/asymp_test_2/1/posterior_chains.csv"),
  log_lik = read.csv("MF_results/asymp_test_2/1/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp2_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
)



# M1 - 1128.76

# asymp 2 summary 
sum_asymp_2_inc = summarise_posterior_chains_not_symp(
  posterior_chains = read.csv("MF_results/asymp_test_2/1/posterior_chains.csv"),
  asymp2= T)

write.csv(sum_asymp_2_inc, "MF_results/comp_main_model/1/sum_asymp_inc_2.csv")

# asymp symp inc --------------------------------------------------------------------

# DIC 
calc_dic(
  posterior_chains = read.csv("MF_results/asymp_and_symp_test/1/posterior_chains.csv"),
  log_lik = read.csv("MF_results/asymp_and_symp_test/1/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_and_symp_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
)



# M1 -  1140.87 


# asymp symp summary 
sum_asymp_symp_inc = summarise_posterior_chains_not_symp(
  posterior_chains =  read.csv("MF_results/asymp_and_symp_test/1/posterior_chains.csv"),
  asymp_symp = T)

write.csv(sum_asymp_symp_inc, "MF_results/comp_main_model/1/sum_asymp_symp_inc.csv")

