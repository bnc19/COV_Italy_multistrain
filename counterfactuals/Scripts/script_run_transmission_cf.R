# Script to run counterfactual M variant transmission scenarios in Veneto ------


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



R0_scales = c(".8", "1", "1.2", "1.4", "1.6")

# Sample from posterior distribution -------------------------------------------

posterior_samples = sample_posterior_chains(
  posterior_chains = posterior_chains,
  number_of_samples = 100
)


posterior_samples.8   = posterior_samples
posterior_samples.1.2 = posterior_samples
posterior_samples.1.4 = posterior_samples
posterior_samples.1.6 = posterior_samples

posterior_samples.8$beta_M   = (posterior_samples$beta_M * .8)
posterior_samples.1.2$beta_M = (posterior_samples$beta_M * 1.2)
posterior_samples.1.4$beta_M = (posterior_samples$beta_M * 1.4)
posterior_samples.1.6$beta_M = (posterior_samples$beta_M * 1.6)

beta_scale_posterior_samples_list = list(
  posterior_samples.8,
  posterior_samples,
  posterior_samples.1.2,
  posterior_samples.1.4,
  posterior_samples.1.6
)


# CF 6:baseline testing, 0.643 sens, rho = est, scale B_M -----------------

model_posts_cf6 = Veneto_posts_cf6  = Veneto_ratio_cf6_list  = list()

for(x in 1:length(R0_scales)){
  
  model_posts_cf6[[x]] =  sapply(1:nrow(posterior_samples), function(i) {
    replicate_rstan_fixed(model = baseline_model,
                          posterior_sample_row = 
                            beta_scale_posterior_samples_list[[x]][i, ],
                          start_date_veneto =  "01-05-2020",
                          time_seed_M_veneto = "01-08-2020")
  })
  
  Veneto_posts_cf6[[x]] = sapply(1:nrow(posterior_samples), function(i) {
      extract_fit_results(posts = model_posts_cf6[[x]][, i],
                          location = "Veneto")
    }, simplify = "array")

  Veneto_ratio_cf6_list[[x]] = calculate_ratio_reported(Veneto_posts_cf6[[x]],
                                                         S0 = 4753625)
}


Veneto_ratio_cf6 = bind_rows(Veneto_ratio_cf6_list) 
Veneto_ratio_cf6$R0_scale = rep(R0_scales, each = nrow(Veneto_ratio_cf6_list[[1]]))

write.csv(Veneto_ratio_cf6, paste0(file_path, "/Veneto_ratio_cf6.csv"))

# CF 7: antigen testing only, 0.643 sens, rho = est, scale B_M -----------------


model_posts_cf7 = Veneto_posts_cf7  = Veneto_ratio_cf7_list  = list()

for(x in 1:length(R0_scales)){
  
  model_posts_cf7[[x]] =  sapply(1:nrow(posterior_samples), function(i) {
    replicate_rstan_fixed(model = antigen_only_model,
                          posterior_sample_row = 
                          beta_scale_posterior_samples_list[[x]][i, ],
                          start_date_veneto =  "01-05-2020",
                          time_seed_M_veneto = "01-08-2020")
  })
  
  Veneto_posts_cf7[[x]] = sapply(1:nrow(posterior_samples), function(i) {
    extract_fit_results(posts = model_posts_cf7[[x]][, i],
                        location = "Veneto")
  }, simplify = "array")
  
  
  Veneto_ratio_cf7_list[[x]] = calculate_ratio_reported(Veneto_posts_cf7[[x]],
                                                        S0 = 4753625)
}


Veneto_ratio_cf7 = bind_rows(Veneto_ratio_cf7_list) 
Veneto_ratio_cf7$R0_scale = rep(R0_scales, each = nrow(Veneto_ratio_cf7_list[[1]]))

write.csv(Veneto_ratio_cf7, paste0(file_path, "/Veneto_ratio_cf7.csv"))


# CF 8: Italy testing, 0.643 sens, rho = est, scale B_M -----------------

model_posts_cf8 = Veneto_posts_cf8  = Veneto_ratio_cf8_list  = list()

for(x in 1:length(R0_scales)){
  
  model_posts_cf8[[x]] =  sapply(1:nrow(posterior_samples), function(i) {
    replicate_rstan_fixed(model = baseline_model,
                          posterior_sample_row = 
                          beta_scale_posterior_samples_list[[x]][i, ],
                          italy_testing_in_veneto = TRUE,
                          start_date_veneto =  "01-05-2020",
                          time_seed_M_veneto = "01-08-2020")
  })
  
  Veneto_posts_cf8[[x]] = sapply(1:nrow(posterior_samples), function(i) {
    extract_fit_results(posts = model_posts_cf8[[x]][, i],
                        location = "Veneto")
  }, simplify = "array")
  
  Veneto_ratio_cf8_list[[x]] = calculate_ratio_reported(Veneto_posts_cf8[[x]],
                                                        S0 = 4753625)
}


Veneto_ratio_cf8 = bind_rows(Veneto_ratio_cf6_list) 
Veneto_ratio_cf8$R0_scale = rep(R0_scales, each = nrow(Veneto_ratio_cf8_list[[1]]))

write.csv(Veneto_ratio_cf8, paste0(file_path, "/Veneto_ratio_cf8.csv"))


