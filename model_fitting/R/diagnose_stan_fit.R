# Function to run diagnostics on a stan fit -----------------------------------# 


# Input:

# - stan_fit: an object of S4 class which contains the model fit results.
# - pars: names of parameters we want to check 
# - file_path: file path to save output 

# Output: 

# - Diagnostic plots
# - data.frame of summary statistics 

diagnose_stan_fit = function(
  fit,
  file_path ,
  pars
){
  
  # required package 
  library(bayesplot)
  library(loo)
  library(rstan)
  
  

  rstan::check_divergences(fit)

  # markov chain trace plots   
  
  stan_fit_post= as.array(fit)
  
  markov_trace = mcmc_trace(stan_fit_post, pars=c("lp__", pars))
  # 
  # bivar = pairs(fit, pars = c("tau", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "rho_it", "rho_ven"))
  # bivar2 = pairs(fit, pars = c("tau", "I0_it[1]", "I0_it[2]", "I0_ven[1]", "I0_ven[2]", "omega[1]", "omega[2]"))

  np = nuts_params(fit)
  mcmc_parcoord=  mcmc_parcoord(stan_fit_post, np = np, pars = pars)
  
  # summary of parameter values, effective sample size and Rhat 
  param_sum = summary(fit, pars = pars)$summary
  
  # save files 
  write.csv(param_sum, file=paste0(file_path,"/param_sum.csv"))
  
  ggsave(
    markov_trace,
    file =  paste0(file_path, "/trace_plot.png"),
    height = 20,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
  ggsave(
    mcmc_parcoord,
    file =  paste0(file_path, "/mcmc_parcoord.png"),
    height = 20,
    width = 50,
    unit = "cm",
    dpi = 720
  )

  # 
  # ggsave(
  #   bivar2,
  #   file =  paste0(file_path, "/bivar2.png"),
  #   height = 20,
  #   width = 50,
  #   unit = "cm",
  #   dpi = 720
  # )
  # WAIC and LOO
  
  log_lik= extract_log_lik(fit)
  waic = waic(log_lik)
  write.csv(waic$estimates,  paste0(file_path, "/WAIC.csv"))
  write.csv(log_lik,  paste0(file_path, "/log_lik.csv"))
  
}
