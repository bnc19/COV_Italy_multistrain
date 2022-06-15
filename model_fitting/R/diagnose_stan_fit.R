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
  file_path
){
  
  # required package 
  library(bayesplot)
  library(loo)
  library(rstan)
  
  

  rstan::check_divergences(fit)

  # markov chain trace plots   
  
  stan_fit_post= as.array(fit)
  
  markov_trace = mcmc_trace(stan_fit_post, pars=c("lp__", pars))
  
  
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
  
  # WAIC and LOO
  
  log_lik= extract_log_lik(fit)
  waic = waic(log_lik)
  write.csv(waic$estimates,  paste0(file_path, "/WAIC.csv"))

  return(list(markov_trace,  param_sum))
  
}
