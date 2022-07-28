# Script to run counterfactual testing scenarios in Veneto ---------------------

# Set up -----------------------------------------------------------------------


# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/sample_posterior_chains.R")
source("R/run_cf.R")
source("R/plot_model_fit.R")


posterior_chains = read.csv("MF_results/symp_test/1/posterior_chains.csv")


file_path = "Results/counterfactuals"
dir.create("Results")
dir.create(paste0(file_path))

baseline_model = stan_model("models/est_test_symp_cf.stan")
antigen_only_model = stan_model("models/cf/est_test_symp_cf_pPCR_0.stan")
molecular_after_antigen_neg_model =stan_model("models/cf/est_test_symp_cf_PCR_after_Ag.stan")


# Sample from posterior distribution -------------------------------------------

posterior_samples = sample_posterior_chains(
  posterior_chains = posterior_chains,
  number_of_samples = 100
)

# posterior_samples_cf2 = posterior_samples
# 
# posterior_samples_cf2$rho_v = 1

# CF 1: antigen testing only, 0.689 sens, rho = est ----------------------------

model_posts_cf1 = sapply(1:nrow(posterior_samples), function(i) {
  replicate_rstan_fixed(model = antigen_only_model,
                        posterior_sample_row = posterior_samples[i, ],
                        scale_time_step = 2,
                        start_date_veneto =  "01-05-2020",
                        time_seed_M_veneto = "01-08-2020")
})



Veneto_posts_cf1 = sapply(1:nrow(posterior_samples), function(i) {
  extract_fit_results(posts = model_posts_cf1[, i],
                      location = "Veneto")
}, simplify = "array")

# summarise results 

Veneto_ratio_cf1 = calculate_ratio_reported(Veneto_posts_cf1, 
                                            S0 = 4753625)

write.csv(Veneto_ratio_cf1, paste0(file_path, "/Veneto_ratio_cf1.csv"))

# # CF 2: antigen testing only, 0.689 sens, rho = 1 ------------------------------
# 
# model_posts_cf2 = sapply(1:nrow(posterior_samples_cf2), function(i) {
#   replicate_rstan_fixed(model = antigen_only_model,
#                         posterior_sample_row = posterior_samples_cf2[i, ],
#                         scale_time_step = 2,
#                         start_date_veneto =  "01-05-2020",
#                         time_seed_M_veneto = "01-08-2020")
# })
# 
# 
# 
# Veneto_posts_cf2 = sapply(1:nrow(posterior_samples_cf2), function(i) {
#   extract_fit_results(posts = model_posts_cf2[, i],
#                       location = "Veneto")
# }, simplify = "array")
# 
# # summarise results 
# 
# 
# Veneto_ratio_cf2 = calculate_ratio_reported(Veneto_posts_cf2, 
#                                             S0 = 4753625)
# 
# write.csv(Veneto_ratio_cf2, paste0(file_path, "/Veneto_ratio_cf2.csv"))
# 

# CF 3: antigen testing only, 0.875 sens, rho = est ------------------------------

model_posts_cf3 = sapply(1:nrow(posterior_samples), function(i) {
  replicate_rstan_fixed(model = antigen_only_model,
                        posterior_sample_row = posterior_samples[i, ],
                        phi_Ag = 0.875,
                        scale_time_step = 2,
                        start_date_veneto =  "01-05-2020",
                        time_seed_M_veneto = "01-08-2020")
})

Veneto_posts_cf3 = sapply(1:nrow(posterior_samples), function(i) {
  extract_fit_results(posts = model_posts_cf3[, i],
                      location = "Veneto")
}, simplify = "array")

# summarise results 


Veneto_ratio_cf3 = calculate_ratio_reported(Veneto_posts_cf3, 
                                            S0 = 4753625)

write.csv(Veneto_ratio_cf3, paste0(file_path, "/Veneto_ratio_cf3.csv"))


# CF 4: molecular follows negative antigen test, 0.689 sens, rho = est ---------

model_posts_cf4 = sapply(1:nrow(posterior_samples), function(i) {
  replicate_rstan_fixed(model = molecular_after_antigen_neg_model,
                        posterior_sample_row = posterior_samples[i, ],
                        scale_time_step = 2,
                        start_date_veneto =  "01-05-2020",
                        time_seed_M_veneto = "01-08-2020")
})


Veneto_posts_cf4 = sapply(1:nrow(posterior_samples), function(i) {
  extract_fit_results(posts = model_posts_cf4[, i],
                      location = "Veneto")
}, simplify = "array")

# summarise results 

Veneto_ratio_cf4 = calculate_ratio_reported(Veneto_posts_cf4, 
                                            S0 = 4753625)

write.csv(Veneto_ratio_cf4, paste0(file_path, "/Veneto_ratio_cf4.csv"))

# CF 5: Italy testing in Veneto, 0.689 sens, rho = est -------------------------


model_posts_cf5 = sapply(1:nrow(posterior_samples), function(i) {
  replicate_rstan_fixed(model = baseline_model,
                        posterior_sample_row = posterior_samples[i, ],
                        italy_testing_in_veneto = T,
                        scale_time_step = 2,
                        start_date_veneto =  "01-05-2020",
                        time_seed_M_veneto = "01-08-2020")
})


Veneto_posts_cf5 = sapply(1:nrow(posterior_samples), function(i) {
  extract_fit_results(posts = model_posts_cf5[, i],
                      location = "Veneto")
}, simplify = "array")

# summarise results 


Veneto_ratio_cf5 = calculate_ratio_reported(Veneto_posts_cf5, 
                                            S0 = 4753625)

write.csv(Veneto_ratio_cf5, paste0(file_path, "/Veneto_ratio_cf5.csv"))

