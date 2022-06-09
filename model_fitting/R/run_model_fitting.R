

run_model_fitting = function(file_path,
                             model_path,
                             # initial conditions
                             n_pop_it = 59257566 - 4847026,
                             n_recov_it = 1482377 - 93401,
                             
                             # index data Italy
                             index_M_it = 6:14,
                             index_A_it = 3:14,
                             index_O_it =  3:9,
                             index_Al_it = 9:14,
                             
                             # initial conditions
                             n_pop_veneto = 4847026,
                             n_recov_veneto = 93401,
                             
                             # index data Veneto
                             index_M_veneto = 8:14,
                             index_A_veneto = c(5, 7:14),
                             index_O_veneto = c(5, 7:9),
                             index_Al_veneto = 9:14,
                             
                             #5.1 days = incubation period (no symptoms)
                             epsilon = 1 / 2.8,
                             ## mean of values given in Vo table S5 assuming R0 = 2.7
                             
                             # 1/ sigma = 5.1 - 2.8
                             sigma = 1 / (5.1 - 2.8),
                             gamma = 1 / 2.1 ,
                             mu = 0.59 ,
                             # probability symptomatic
                             phi_PCR = 0.920,
                             phi_Ag =  0.643,
                             
                             # Dates Italy
                             start_date_it =  "01-05-2020",
                             end_date_it = "31-05-2021",
                             time_intervention_it = c("15-11-2020" , "15-03-2021") ,
                             
                             time_seed_alpha_it = "01-11-2020",
                             time_seed_M_it = "01-08-2020" ,
                             time_vac_it = "27-12-2020",
                             # date first vaccine (first datapoint)
                             
                             # Dates Veneto
                             start_date_veneto =  "01-07-2020",
                             end_date_veneto = "31-05-2021",
                             time_intervention_veneto = c("15-11-2020" , "15-03-2021"),
                             
                             time_seed_alpha_veneto = "01-11-2020",
                             time_seed_M_veneto = "01-10-2020",
                             time_vac_veneto = "27-12-2020",
                             # date first vaccine (first datapoint),
                             n_chains = 4,
                             n_warmups = 2000,
                             n_iter = 4000,
                             n_thin = 1,
                             pars = c("lp__",
                                      "beta",
                                      "rho_it" ,
                                      "rho_ven" ,
                                      "omega",
                                      "I0_it",
                                      "I0_ven",
                                      "kappa"),
                             seed_values = c(1, 4, 5, 97)
                                   ){
  
# Create file to save folder 
  
  folder = file_path
  
  if (file.exists(folder)) {
    
    cat("The folder already exists")
    
  } else {
    
    dir.create(folder)
    
  }

# Load packages 
  library(bayesplot)
  library(tidyverse)
  library(Hmisc)
  library(cowplot)
  
# Source function 
  source("R/sample_stan_model.R")
  source("R/diagnose_stan_fit.R")
  source("R/plot_model_fit.R")
  
  
# Run model 
  fit = sample_stan_model(modelPath=model_path,
                       n_pop_it = n_pop_it,
                       n_recov_it=n_recov_it,
                       index_M_it = index_M_it,
                       index_A_it = index_A_it,
                       index_O_it =  index_O_it,
                       index_Al_it = index_Al_it,
                       n_pop_veneto= n_pop_veneto,
                       n_recov_veneto= n_recov_veneto,
                       index_M_veneto = index_M_veneto,
                       index_A_veneto = index_A_veneto,
                       index_O_veneto = index_O_veneto,
                       index_Al_veneto = index_Al_veneto,
                       epsilon = epsilon, 
                       sigma = sigma,
                       gamma = gamma ,
                       mu = mu ,
                       phi_PCR = phi_PCR,
                       phi_Ag =  phi_Ag,
                       start_date_it = start_date_it, 
                       end_date_it = end_date_it,
                       time_intervention_it= time_intervention_it,
                       time_seed_alpha_it = time_seed_alpha_it,
                       time_seed_M_it= time_seed_M_it ,
                       time_vac_it = time_vac_it,
                       start_date_veneto=  start_date_veneto,
                       end_date_veneto= end_date_veneto,
                       time_intervention_veneto= time_intervention_veneto, 
                       time_seed_alpha_veneto= time_seed_alpha_veneto,
                       time_seed_M_veneto= time_seed_M_veneto,
                       time_vac_veneto= time_vac_veneto,
                       n_chains=n_chains ,
                       n_warmups=n_warmups,
                       n_iter=n_iter ,
                       n_thin=n_thin ,
                       pars=pars ,
                       seed_values=seed_values)
  
# Diagnostics 
  diagnostics = diagnose_stan_fit(fit=fit,file_path = file_path)
  
# Plot 
  plot = plot_stan_fit(fit=fit,file_path =file_path)
  
  return(list(fit,diagnostics,plot))
}
