# Script to calculate return parameter estimates for all sensitiity analyses ---


# Set up -----------------------------------------------------------------------

# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

library(tidyverse)

source("R/sample_posterior_chains.R")
source("R/calculate_R0.R")

# 64.3% sens (2) --------------------------------------------------------------------

sum_symp_inc_2 = summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/2/posterior_chains.csv"))


write.csv(sum_symp_inc_2, "MF_results/comp_main_model/sum_symp_inc_2.csv")

# 87.5% sens (3) --------------------------------------------------------------------

# asymp summary 
sum_symp_inc_3 = summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/3/posterior_chains.csv"))

write.csv(sum_symp_inc_3, "MF_results/comp_main_model/sum_symp_inc_3.csv")


# 5.6 days exposure to symp (4) -----------------------------------------------

# asymp summary 
sum_symp_inc_4 = summarise_posterior_chains_symp(
  sigma = 1 / (5.6 - 1.31),
  gamma = 1 / 1.6 ,
  posterior_chains = read.csv("MF_results/symp_test/4/posterior_chains.csv"))

write.csv(sum_symp_inc_4, "MF_results/comp_main_model/sum_symp_inc_4.csv")




# 2.15 day inc  (5) --------------------------------------------------------------------

# asymp summary 
sum_symp_inc_5= summarise_posterior_chains_symp(
  sigma = 1 / (5.1 - 2.15),
  gamma = 1 / 2.1 ,
  posterior_chains = read.csv("MF_results/symp_test/5/posterior_chains.csv"))

write.csv(sum_symp_inc_5, "MF_results/comp_main_model/sum_symp_inc_5.csv")


# 70% VE  (6) --------------------------------------------------------------------

# asymp summary 
sum_symp_inc_6 = summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/6/posterior_chains.csv"))

write.csv(sum_symp_inc_6, "MF_results/comp_main_model/sum_symp_inc_6.csv")

# 80% VE  (7) --------------------------------------------------------------------

# asymp summary 
sum_symp_inc_7= summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/7/posterior_chains.csv"))

write.csv(sum_symp_inc_7, "MF_results/comp_main_model/sum_symp_inc_7.csv")


# alpha (8) --------------------------------------------------------------------

# asymp summary 
sum_symp_inc_8 = summarise_posterior_chains_symp(
  posterior_chains = read.csv("MF_results/symp_test/8/posterior_chains.csv"),
  alpha = TRUE)

write.csv(sum_symp_inc_8, "MF_results/comp_main_model/sum_symp_inc_8.csv")

