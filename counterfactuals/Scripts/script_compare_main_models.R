################################################################################
# Script to calculate the log-lik, posterior means and 95% CrI for for all model   #
# variants (Table 2, Table S6)                                                 #
################################################################################


# Set up -----------------------------------------------------------------------

# rm(list = ls())
# setwd("Q:/COV_Italy_multistrain/counterfactuals")

# packages
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# functions 
source("R/calc_LL.R")
source("R/sample_posterior_chains.R")

# output 
file_path = "Results/summarise_models/1"
dir.create("Results")
dir.create("Results/summarise_models")
dir.create(paste0(file_path))

# symp inc ---------------------------------------------------------------------

# LL
calculate_log_like(read.csv("MF_results/symp_test/1/log_lik.csv")[, -1])

# symp summary
sum_symp_inc = summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/1/posterior_chains.csv"))

write.csv(sum_symp_inc, paste0(file_path, "/sum_symp_inc.csv"))

# asymp inc --------------------------------------------------------------------

# LL
calculate_log_like(read.csv("MF_results/asymp_test/1/log_lik.csv")[, -1])


# asymp summary
sum_asymp_inc = summarise_posterior_chains_not_symp(
  posterior_chains = read.csv("MF_results/asymp_test/1/posterior_chains.csv"),
  asymp = T
)

write.csv(sum_asymp_inc, paste0(file_path, "/sum_asymp_inc.csv"))


# asymp 2 inc -------------------------------------------------------------------

# LL
calculate_log_like(read.csv("MF_results/asymp_test_2/1/log_lik.csv")[, -1])


# asymp 2 summary
sum_asymp_2_inc = summarise_posterior_chains_not_symp(
  posterior_chains = read.csv("MF_results/asymp_test_2/1/posterior_chains.csv"),
  asymp2 = T
)

write.csv(sum_asymp_2_inc, paste0(file_path, "/sum_asymp2_inc.csv"))

# asymp symp inc --------------------------------------------------------------------

# LL
calculate_log_like(read.csv("MF_results/asymp_and_symp_test/1/log_lik.csv")[, -1])

# asymp symp summary
sum_asymp_symp_inc = summarise_posterior_chains_not_symp(
  posterior_chains =  read.csv("MF_results/asymp_and_symp_test/1/posterior_chains.csv"),
  asymp_symp = T
)

write.csv(sum_asymp_symp_inc,
          paste0(file_path, "/sum_asymp_symp_inc.csv"))
