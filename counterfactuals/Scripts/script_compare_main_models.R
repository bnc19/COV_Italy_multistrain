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
calc_dic(
  posterior_chains = read.csv("MF_results/inc/symp_test/17/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/symp_test/17/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_symp_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
  # ,
  # sigma = 1 / (5.6 - 1.31),
  # gamma = 1 / (7.2 - 5.6)
)

# 1 1196.377

# 6 1101.816 

# M17 - (baseline) - 1121.50   NO DIV, "Rhat good" :) 
# 16 - 1024.95


# 10 68.9% sens 1146.156  
# 11  1017.97
# 12 1091.8
# 13 1573.40
# 14 1108.9744 

waic(as.matrix(read.csv("MF_results/inc/symp_test/6/log_lik.csv")[,-1]))


# symp summary 
sum_symp_inc = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/inc/symp_test/17/posterior_chains.csv"))

write.csv(sum_symp_inc, "MF_results/comp_main_model/17/sum_symp_inc.csv")



# asymp inc --------------------------------------------------------------------

# DIC  / WAIC
DIC_asymp_inc = calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_test/6/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_test/6/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
)

DIC_asymp_inc 
# 1  1196.533
# 2  1163.329 
# 6   1174.052  - -564.195

waic(as.matrix(read.csv("MF_results/inc/asymp_test/6/log_lik.csv")[,-1]))



# asymp summary 
sum_asymp_inc = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/inc/asymp_test/6/posterior_chains.csv"))

write.csv(sum_asymp_inc, "MF_results/comp_main_model/6/sum_asymp_inc.csv")


# asymp 2 inc -------------------------------------------------------------------

# DIC  / WAIC
DIC_asymp_inc_2 = calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_test_2/9/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_test_2/9/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp2_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
)


DIC_asymp_inc_2 

#  1 - 1195.755
#  2 - 1158.526
#  6 - 1128.766   - -547.1309
#  9 - 1139.3618 
waic(as.matrix(read.csv("MF_results/inc/asymp_test_2/6/log_lik.csv")[,-1]))

#   

# asymp 2 summary 
sum_asymp_2_inc = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/inc/asymp_test_2/6/posterior_chains.csv"))

write.csv(sum_asymp_2_inc, "MF_results/comp_main_model/6/sum_asymp_inc_2.csv")

# asymp symp inc --------------------------------------------------------------------

# DIC / WAIC
DIC_asymp_symp_inc = calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_and_symp_test/6/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_and_symp_test/6/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_and_symp_cf.stan",
  scale_time_step = 2,
  start_date_veneto =  "01-05-2020",
  time_seed_M_veneto = "01-08-2020",
  phi_Ag =  0.689
)


DIC_asymp_symp_inc 

# 1 - 1225.759
# 2 - 1160.763
# 6 - 1138.211 
waic(as.matrix(read.csv("MF_results/inc/asymp_and_symp_test/6/log_lik.csv")[,-1]))


# asymp symp summary 
sum_asymp_symp_inc = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/inc/asymp_and_symp_test/6/posterior_chains.csv"))

write.csv(sum_asymp_symp_inc, "MF_results/comp_main_model/6/sum_asymp_symp_inc.csv")


# Prevalence models ------------------------------------------------------------

# symp prev --------------------------------------------------------------------

#DIC / WAIC

DIC_symp_prev = calc_dic(
posterior_chains = read.csv("MF_results/prev/symp_test/same_omega_bound_seed/posterior_chains.csv"),
log_lik = read.csv("MF_results/prev/symp_test/same_omega_bound_seed/log_lik.csv")[,-1],
fixed_model_path = "models/est_test_symp_prev_cf.stan",
scale_time_step = 5
)

DIC_symp_prev # 723.8621


waic(as.matrix(read.csv("MF_results/prev/symp_test/same_omega_bound_seed/log_lik.csv")[,-1]))

# symp prev summary 
sum_symp_prev = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/prev/symp_test/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_symp_prev, "CF_results/comp_main_model/sum_symp_prev.csv")



# asymp prev ------------------------------------------------------------------------------

#DIC / WAIC
DIC_asymp_prev = calc_dic(
  posterior_chains = read.csv("MF_results/prev/asymp_test/same_omega_bound_seed/posterior_chains.csv"),
  log_lik = read.csv("MF_results/prev/asymp_test/same_omega_bound_seed/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_cf.stan",
  scale_time_step = 5
)

DIC_asymp_prev #   990.0746



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
  scale_time_step = 5
)

DIC_asymp_prev_2 # 1001.813

waic(as.matrix(read.csv("MF_results/prev/asymp_test_2/same_omega_bound_seed/log_lik.csv")[,-1]))

#   

# asymp 2 summary 
sum_asymp_2_prev = summarise_posterior_chains(
  posterior_chains = read.csv("MF_results/prev/asymp_test_2/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_asymp_2_prev, "CF_results/comp_main_model/sum_asymp_prev_2.csv")


# asymp symp prev --------------------------------------------------------------

# DIC / WAIC
DIC_asymp_symp_prev = calc_dic(
  posterior_chains = read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/posterior_chains.csv"),
  log_lik = read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp_and_symp_prev_cf.stan",
  scale_time_step = 5
)

DIC_asymp_symp_prev # 949.34

waic(as.matrix(read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/log_lik.csv")[,-1]))

# asymp symp summary 
sum_asymp_symp_prev = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/prev/asymp_and_symp_test/same_omega_bound_seed/posterior_chains.csv"))

write.csv(sum_asymp_symp_prev, "CF_results/comp_main_model/sum_asymp_symp_prev.csv")










# test 
DIC_test= calc_dic(
  posterior_chains = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/same_rho/posterior_chains.csv"),
  log_lik = read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/same_rho/log_lik.csv")[,-1],
  fixed_model_path = "models/est_test_asymp2_cf.stan",
  scale_time_step = 2
)

DIC_test

waic(as.matrix(read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/same_rho/log_lik.csv")[,-1]))


# asymp symp summary 
sum_test = summarise_posterior_chains(
  posterior_chains =  read.csv("MF_results/inc/asymp_test_2/same_omega_bound_seed/same_rho/posterior_chains.csv"))


