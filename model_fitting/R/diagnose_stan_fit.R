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
  pars = c("lp__", 
           "beta[1]","beta[2]","beta[3]","beta[4]",
           "rho_it" , "rho_ven" ,
           "omega[1]", "omega[2]","omega[3]","omega[4]",
           "I0_it[1]", "I0_it[2]","I0_it[3]","I0_it[4]",
           "I0_ven[1]", "I0_ven[2]","I0_ven[3]","I0_ven[4]",
           "k"),
  file_path
){
  
  # required package 
  library(bayesplot)
  library(loo)
  
  
  stan_fit_post= as.array(fit)
  
  # markov chain trace plots   
  markov_trace = mcmc_trace(stan_fit_post, pars=c("lp__", pars))
  
  # bivariate marginal posterior distributions
  bivar = pairs(fit, pars = pars[1:round(length(pars)/2)], cex.labels=1.5, font.labels=9, condition = "accept_stat__")  
  bivar2 = pairs(fit, pars = pars[(round(length(pars)/2)+1):length(pars)], cex.labels=1.5, font.labels=9, condition = "accept_stat__")  
  
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
    bivar,
    file =  paste0(file_path,"/bivar_plot.png"),
    height = 20,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
  
  ggsave(
    bivar2,
    file =  paste0(file_path,"/bivar2_plot.png"),
    height = 20,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
  
  # WAIC and LOO
  
  log_lik= extract_log_lik(fit)
  waic = waic(log_lik)
  print(waic)
  # loo
  loo =loo(fit, save_psis = TRUE)
  
  print(loo)
  

  return(list(markov_trace, bivar, param_sum))
  
}
