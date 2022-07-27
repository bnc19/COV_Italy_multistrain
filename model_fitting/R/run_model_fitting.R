

run_model_fitting = function(file_path,
                             model_path,
                             scale_time_step = 1,
                             n_pop_it = 59257566 - 4847026,
                             n_recov_it = 1482377 - 93401,
                             index_M_it = 6:14,
                             index_A_it = 3:14,
                             index_O_it =  3:9,
                             index_Al_it = 9:14,
                             n_pop_veneto = 4847026,
                             n_recov_veneto = 93401,
                             index_M_veneto = 8:14,
                             index_A_veneto = c(5,7:14),
                             index_O_veneto = c(5,7:9),
                             index_Al_veneto = 9:14,
                             epsilon = 1 / 1.31,
                             sigma = 1 / (5.1 - 1.31),
                             gamma = 1 / 2.1 ,
                             mu = 0.59 ,
                             phi_PCR = 0.920,
                             phi_Ag =  0.689,
                             VE = 1, 
                             start_date_it =  "01-05-2020",
                             end_date_it = "31-05-2021",
                             time_intervention_it = c("15-11-2020" , "15-03-2021") ,
                             time_seed_alpha_it = "01-11-2020",
                             time_seed_M_it = "01-08-2020" ,
                             time_vac_it = "27-12-2020",
                             start_date_veneto =  "01-07-2020",
                             end_date_veneto = "31-05-2021",
                             time_intervention_veneto = c("15-11-2020" , "15-03-2021"),
                             time_seed_alpha_veneto = "01-11-2020",
                             time_seed_M_veneto = "01-10-2020",
                             time_vac_veneto = "27-12-2020",
                             n_chains = 4,
                             n_warmups = 500,
                             n_iter = 1000,
                             n_thin = 1,
                             pars = c("lp__", 
                                      "beta[1]","beta[2]","beta[3]","beta[4]",
                                      "rho_it" , "rho_ven" ,
                                      "omega[1]", "omega[2]","omega[3]","omega[4]",
                                      "I0_it[1]", "I0_it[2]","I0_it[3]","I0_it[4]",
                                      "I0_ven[1]", "I0_ven[2]","I0_ven[3]","I0_ven[4]",
                                      "k"),
                             seed_values = c(1, 4, 5, 97),
                             prev ,
                             adapt_delta = 0.8
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
  library(loo)
  
# Source function 
  source("R/sample_stan_model.R")
  source("R/diagnose_stan_fit.R")
  source("R/plot_model_fit.R")
  source("R/save_post_chains.R")
  
  
# Run model 
  fit = sample_stan_model(
                       modelPath=model_path,
                       scale_time_step = scale_time_step,
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
                       VE = VE, 
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
                       seed_values=seed_values,
                       prev = prev,
                       adapt_delta=adapt_delta)
  
  
  # Plot 
  
  plot = plot_stan_fit(fit=fit,file_path =file_path,
                       n_pop_it=n_pop_it,
                       n_recov_it=n_recov_it,
                       n_pop_veneto=n_pop_veneto,
                       n_recov_veneto=n_recov_veneto,
                       index_M_it=index_M_it,
                       index_A_it=index_A_it,
                       index_O_it =index_O_it,
                       index_Al_it=index_Al_it,
                       index_M_veneto=index_M_veneto,
                       index_A_veneto =index_A_veneto,
                       index_O_veneto=index_O_veneto,
                       index_Al_veneto=index_Al_veneto,  
                       start_date_it=start_date_it, 
                       end_date_it=end_date_it,
                       start_date_veneto=start_date_veneto,
                       end_date_veneto=end_date_veneto)
  
  # Save posterior 
  
  save_post_chains(fit,file_path = file_path)
  
  

  # Diagnostics 
  
  
  check_divergences(fit)
  
  
  
  log_lik= extract_log_lik(fit)
  write.csv(log_lik,  paste0(file_path, "/log_lik.csv"))
  

  stan_fit_post= as.array(fit)
  
  markov_trace = mcmc_trace(stan_fit_post, pars=c("lp__", pars))

 
  ggsave(
    markov_trace,
    file =  paste0(file_path, "/trace_plot.png"),
    height = 20,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
 # diagnostics = diagnose_stan_fit(fit=fit,pars=pars,file_path = file_path)  
  
 return(fit)
  
}
