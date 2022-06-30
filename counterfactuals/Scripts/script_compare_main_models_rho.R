# Script to calculate DIC and return parameter estimates for all models --------


# Set up -----------------------------------------------------------------------


# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

library(tidyverse)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/calc_dic.R")
source("R/sample_posterior_chains.R")


# Incidence models -------------------------------------------------------------

# symp inc ---------------------------------------------------------------------

# inc DIC / WAIC
DIC_symp_inc = calc_dic(
  posterior_chains = read.csv("MF_results/inc/symp_test/same_omega_bound_seed/bound_rho/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/symp_test/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_symp_cf.stan",
  scale_time_step = 5
)
DIC_symp_inc # 729.3521

waic(as.matrix(read.csv("MF_results/inc/symp_test/same_omega_bound_seed/log_lik.csv")[,-1]))
# 738.6

# symp summary 
sum_symp_inc = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/inc/symp_test/same_omega_bound_seed/bound_rho/posterior_chains.csv"))

write.csv(sum_symp_inc, "CF_results/comp_main_model/sum_symp_inc_bound_rho.csv")


# asymp inc --------------------------------------------------------------------

# DIC  / WAIC
DIC_asymp_inc = calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_test/same_omega_bound_seed/bound_rho/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_test/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_cf.stan",
  scale_time_step = 5
)

DIC_asymp_inc #  1064.624


waic(as.matrix(read.csv("MF_results/inc/asymp_test/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1]))
# 1512.3 


# asymp summary 
sum_asymp_inc = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/inc/asymp_test/same_omega_bound_seed/bound_rho/posterior_chains.csv"))

write.csv(sum_asymp_inc, "CF_results/comp_main_model/sum_asymp_inc_bound_rho.csv")


# asymp 2 inc -------------------------------------------------------------------

# DIC  / WAIC
DIC_asymp_inc_2 = calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp2_cf.stan",
  scale_time_step = 5
)

DIC_asymp_inc_2 #  1063.687


waic(as.matrix(read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1]))

#   1407.2 

# asymp 2 summary 
sum_asymp_2_inc = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/posterior_chains.csv"))

write.csv(sum_asymp_2_inc, "CF_results/comp_main_model/sum_asymp_inc_2_bound_rho.csv")

# asymp symp inc --------------------------------------------------------------------

# DIC / WAIC
DIC_asymp_symp_inc = calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_and_symp_test/same_omega_bound_seed/bound_rho/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_and_symp_test/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_and_symp_cf.stan",
  scale_time_step = 5
)


DIC_asymp_symp_inc # 1067.064

waic(as.matrix(read.csv("MF_results/inc/asymp_and_symp_test/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1]))
# 1372.9 

# asymp symp summary 
sum_asymp_symp_inc = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/inc/asymp_and_symp_test/same_omega_bound_seed/bound_rho/posterior_chains.csv"))

write.csv(sum_asymp_symp_inc, "CF_results/comp_main_model/sum_asymp_symp_inc_bound_rho.csv")


# Prevalence models ------------------------------------------------------------

# symp prev --------------------------------------------------------------------

#DIC / WAIC

DIC_symp_prev = calc_dic(
posterior_chains = read.csv("MF_results/prev/symp_test/same_omega_bound_seed/posterior_chains.csv"),
log_lik = read.csv("MF_results/prev/symp_test/same_omega_bound_seed/log_lik.csv")[,-1],
fixed_model_path = "models/est_test_symp_prev_cf.stan",
scale_time_step = 5
)

DIC_symp_prev # 724.28


waic(as.matrix(read.csv("MF_results/prev/symp_test/same_omega_bound_seed/log_lik.csv")[,-1]))

# symp prev summary 
sum_symp_prev = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/prev/symp_test/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_symp_prev, "CF_results/comp_main_model/sum_symp_prev.csv")



# asymp prev DIC / WAIC
DIC_asymp_prev = calc_dic(
  posterior_chains = read.csv("MF_results/prev/asymp_test/same_omega_bound_seed/posterior_chains.csv"),
  log_lik = read.csv("MF_results/prev/asymp_test/same_omega_bound_seed/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_cf.stan",
  scale_time_step = 2
)

DIC_asymp_prev #  989.9785



waic(as.matrix(read.csv("MF_results/prev/asymp_test/same_omega_bound_seed/log_lik.csv")[,-1]))

# asymp prev summary 
sum_asymp_prev = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/prev/asymp_test/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_asymp_prev, "CF_results/comp_main_model/sum_asymp_prev.csv")


# asymp 2 prev------------------------------------------------------------------

# DIC  / WAIC
DIC_asymp_prev_2 = calc_dic(
  posterior_chains = read.csv("MF_results/prev/asymp_test_2/same_omega_bound_seed/posterior_chains.csv"),
  log_lik = read.csv("MF_results/prev/asymp_test_2/same_omega_bound_seed/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_prev2_cf.stan",
  scale_time_step = 4
)

DIC_asymp_prev_2 # 1001.813

waic(as.matrix(read.csv("MF_results/prev/asymp_test_2/same_omega_bound_seed/log_lik.csv")[,-1]))

#   

# asymp 2 summary 
sum_asymp_2_prev = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/prev/asymp_test_2/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_asymp_inc, "CF_results/comp_main_model/sum_asymp_prev_2.csv")


# asymp symp prev --------------------------------------------------------------

# DIC / WAIC
DIC_asymp_symp_prev = calc_dic(
  posterior_chains = read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/posterior_chains.csv"),
  log_lik = read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_and_symp_prev_cf.stan",
  scale_time_step = 2
)

DIC_asymp_symp_prev # 890.5929

waic(as.matrix(read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/log_lik.csv")[,-1]))

# asymp symp summary 
sum_asymp_symp_prev = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_asymp_symp_prev, "CF_results/comp_main_model/sum_asymp_symp_prev.csv")





# test 
DIC_test= calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp2_cf.stan",
  scale_time_step = 2
)

DIC_test


# asymp symp summary 
sum_test = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/bound_rho/posterior_chains.csv"))


