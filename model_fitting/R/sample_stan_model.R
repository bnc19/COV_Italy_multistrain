

draw_init_values = function(seed=1){
  
  set.seed(seed)
  list(  beta = replicate(4,runif(1,0,3)),
         I0_it = replicate(4, runif(1, 1,100)),
         I0_ven = replicate(4, runif(1, 1,100)),
         #omega = replicate(4,runif(1,0.2,0.8)),
         rho_it = runif(1,0.2,0.8),
         #rho_ven = runif(1,0.2,0.8),
         k = runif(1,0.01,1)
  )}

draw_init_values_prev = function(seed=1){
  
  set.seed(seed)
  list(  beta = replicate(4,runif(1,0,1)),
         I0_it = replicate(4, runif(1, 1,100)),
         I0_ven = replicate(4, runif(1, 1,100)),
         omega = replicate(4,runif(1,0.2,0.8)),
         rho_it = runif(1,0.2,0.8),
         rho_ven = runif(1,0.2,0.8),
         k = runif(1,0.01,2),
         tau = runif(1,1,2)
  )}

sample_stan_model = function(
  modelPath,
  n_chains ,
  n_warmups,
  n_iter ,
  n_thin ,
  pars ,
  seed_values,
  n_pop_it ,
  scale_time_step,
  n_recov_it,
  index_M_it ,
  index_A_it ,
  index_O_it ,
  index_Al_it ,
  n_pop_veneto,
  n_recov_veneto,
  index_M_veneto ,
  index_A_veneto ,
  index_O_veneto ,
  index_Al_veneto , 
  epsilon , 
  sigma ,
  gamma  ,
  mu  ,
  phi_PCR ,
  phi_Ag,
  start_date_it , 
  end_date_it ,
  time_intervention_it,
  time_seed_alpha_it,
  time_seed_M_it ,
  time_vac_it , 
  start_date_veneto,
  end_date_veneto,
  time_intervention_veneto, 
  time_seed_alpha_veneto,
  time_seed_M_veneto,
  time_vac_veneto,
  prev,
  adapt_delta
){
  
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  source("R/format_stan_data.R")
  
  list_of_inits = list()
  
  if (prev==F){
  for(i in 1:n_chains)  {
    list_of_inits[[i]] =  draw_init_values(seed = seed_values[i]
    )}
  } else {
    for(i in 1:n_chains)  {
      list_of_inits[[i]] =  draw_init_values_prev(seed = seed_values[i]
      )}
  }

  model = stan_model(paste(modelPath))
  
  
  stan_data = format_stan_data(
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
    time_vac_veneto= time_vac_veneto
  )
  
  time.start = Sys.time()
  
  fit = sampling(
    model,
    data = stan_data,
    init = list_of_inits,
    chains =n_chains ,
    warmup = n_warmups,
    iter =n_iter ,
    thin =n_thin,
    control = list(
      adapt_delta = adapt_delta, 
      max_treedepth = 15
    )
  ) 
  
  time.end = Sys.time()
  
  print(time.end - time.start)
  
  
  check_divergences(fit)

  return(fit)

  }
