# Set up -----------------------------------------------------------------------


# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/calc_dic.R")



# symp prev DIC 
DIC = calc_dic(
posterior_chains = read.csv("Results/prev_scale_5_non_log_1000_it_no_div/symp_test/posterior_chains.csv"),
log_lik = read.csv("Results/prev_scale_5_non_log_1000_it_no_div/symp_test/log_lik.csv")[,-1],
fixed_model_path = "models/est_test_symp_prev_cf.stan"
)


