

draw_init_values = function(seed=1){
  
  set.seed(seed)
  list(  beta = replicate(4,runif(1,0,3)),
         I0_it = replicate(4, runif(1, 1,20)),
         I0_ven = replicate(4, runif(1, 1,20)),
         omega = replicate(4,runif(1,0.2,0.8)),
         rho_it = runif(1,0.2,0.8),
         rho_ven = runif(1,0.2,0.8),
         k = runif(1,0.01,2)
  )}


run_stan_model = function(
  modelPath,
  n_chains =4,
  n_warmups =2000,
  n_iter = 4000,
  n_thin = 1,
  pars = c("lp__", "beta", "rho_it" , "rho_ven" ,
           "omega", "I0_it","I0_ven", "kappa"), 
  seed_values = c(1,4,5,97)
){
  
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  source("R/format_stan_data.R")
  
  list_of_inits = list()
  
  for(i in 1:n_chains)  {
    list_of_inits[[i]] =  draw_init_values(seed = seed_values[i]
    )}
  

  model = stan_model(paste(modelPath))
  stan_data = format_stan_data()
  
  time.start = Sys.time()
  
  fit = sampling(
    model,
    data = stan_data,
    init = list_of_inits,
    control = list(
      adapt_delta = 0.85, 
      max_treedepth = 12
    )
  ) 
  
  time.end = Sys.time()
  
  print(time.end - time.start)

  return(fit)

  }
