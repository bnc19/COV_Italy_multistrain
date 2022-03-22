

replicate_rstan_fixed = function(posterior,
                                 n_difeq,
                                 n_pop ,
                                 n_recov ,
                                 n_months,
                                 gamma,
                                 sigma,
                                 phi_PCR,
                                 phi_Ag,
                                 n_days,
                                 daily_Ag_i,
                                 daily_PCR_i,
                                 average_daily_vaccination_i,
                                 x_i_data ,
                                 SEIR_model
                                )
                                {
  
  
  rho = as.numeric(posterior[9])
  beta = as.numeric(posterior[c(6,5,7,8)])
  omega = as.numeric(posterior[10:11])
  seed =  as.numeric(posterior[13:16])
  k  =  as.numeric(posterior[12])

  model_data = list(
    
    n_difeq = n_difeq,
    n_pop = n_pop,
    n_days = n_days,
    n_recov = n_recov,
    n_months = n_months,
    
    gamma = gamma,
    sigma = sigma,
    rho = rho,
    beta = beta, 
    seed= seed, 
    t0 = 0,
    ts = 1:n_days,
    phi_PCR = phi_PCR,
    phi_Ag = phi_Ag,
    
    Ag_daily  = daily_Ag_i,
    PCR_daily = daily_PCR_i,
    
    x_r_data = c(average_daily_vaccination_i, sigma, gamma, rho, beta, omega, seed,  phi_PCR,phi_Ag),
    x_i_data = x_i_data

  )


  
  
  SEIR__fit = sampling(
    SEIR_model,
    data = model_data,
    chains = 1,
    iter = 1,
    algorithm = "Fixed_param"
    
  )
  
  
  SEIR_fit_posts =  rstan::extract(SEIR__fit)
  
  
  fit_SEIR_M = SEIR_fit_posts$reported_incidence[,, 1]
  fit_SEIR_A = SEIR_fit_posts$reported_incidence[,, 2]
  fit_SEIR_O = SEIR_fit_posts$reported_incidence[,, 3]
  fit_SEIR_Al = SEIR_fit_posts$reported_incidence[,, 4]
  
  true_SEIR_M = SEIR_fit_posts$true_incidence[,, 1]
  true_SEIR_A = SEIR_fit_posts$true_incidence[,, 2]
  true_SEIR_O = SEIR_fit_posts$true_incidence[,, 3]
  true_SEIR_Al = SEIR_fit_posts$true_incidence[,, 4]
  
 total_incidence =  c(SEIR_fit_posts$total_incidence)
 total_reported_incidence =  c(SEIR_fit_posts$total_reported_incidence)

 pPCR =  c(SEIR_fit_posts$pPCR_daily)
  
 R0_M = SEIR_fit_posts$R_0[,, 1]
 R0_A = SEIR_fit_posts$R_0[, ,2]
 R0_O = SEIR_fit_posts$R_0[, ,3]
 R0_Al = SEIR_fit_posts$R_0[, ,4]
 

  out = as.matrix(data.frame(fit_SEIR_M, fit_SEIR_A, fit_SEIR_O, fit_SEIR_Al,
                             true_SEIR_M, true_SEIR_A, true_SEIR_O, true_SEIR_Al, 
                             total_incidence, total_reported_incidence, 
                             R0_M, R0_A, R0_O, R0_Al, pPCR
                             ))  # reported incidence, plot against data
  

  return(out)
  
}



################ change reporting rate to 1 e.g. mass antigen testing ################ 


change_rho = function(data.frame){
  data.frame$rho_chain = 1
  return(data.frame)
}
