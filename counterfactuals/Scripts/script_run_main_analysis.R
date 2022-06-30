# Script to run the main analysis (baseline testing) in Veneto and rest of Italy 


# Set up -----------------------------------------------------------------------


# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")
library(scales)  
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/sample_posterior_chains.R")
source("R/run_cf.R")
source("R/plot_model_fit.R")


posterior_chains = read.csv("MF_results/inc/symp_test/6/posterior_chains.csv")

file_path = "MF_results/baseline"
dir.create(paste0(file_path))

baseline_model = stan_model("models/est_test_symp_cf.stan")
# Sample from posterior distribution -------------------------------------------

posterior_samples = sample_posterior_chains(
  posterior_chains = posterior_chains,
  number_of_samples = 100
)

parameter_summary = summarise_posterior_chains(posterior_samples)


# Run model across posterior samples -------------------------------------------


model_posts = sapply(1:nrow(posterior_samples), function(i) {
  replicate_rstan_fixed(model = baseline_model,
                        posterior_sample_row = posterior_samples[i, ],
                        scale_time_step = 2,
                        start_date_veneto =  "01-05-2020",
                        time_seed_M_veneto = "01-08-2020")
})

Italy_posts = sapply(1:nrow(posterior_samples), function(i) {
  extract_fit_results(posts = model_posts[, i],
                      location = "Italy",
                      baseline = T)
}, simplify = "array")


Veneto_posts = sapply(1:nrow(posterior_samples), function(i) {
  extract_fit_results(posts = model_posts[, i],
                      location = "Veneto",
                      baseline = T)
}, simplify = "array")


# summarise results and plot ---------------------------------------------------

Veneto_post_df = summarise_results(Veneto_posts,
                                   start_date = "01-05-2020",
                                   end_date = "31-05-2021", 
                                   S0 = 4753625)

Italy_post_df = summarise_results(Italy_posts,
                                  start_date = "01-05-2020",
                                  end_date = "31-05-2021",
                                  S0 = 53021564)



Veneto_ratio = calculate_ratio_reported(Veneto_posts,  S0 = 4753625)
Veneto_ratio$Location = "Veneto"
Italy_ratio = calculate_ratio_reported(Italy_posts,S0 = 53021564)
Italy_ratio$Location = "Rest of Italy"

# Save files  ---------------------------------------------------------
write.csv(Veneto_post_df, paste0(file_path, "/Veneto_post_df.csv"))
write.csv(Italy_post_df, paste0(file_path, "/Italy_post_df.csv"))
write.csv(Veneto_ratio, paste0(file_path, "/Veneto_ratio.csv"))
write.csv(Italy_ratio, paste0(file_path, "/Italy_ratio.csv"))

