# DIC functions

calc_dic = function (
  posterior_chains,
  log_lik,
  fixed_model_path,
  scale_time_step,
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
  phi_Ag =  0.643,
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
  italy_testing_in_veneto = F
){
  
  source("R/sample_posterior_chains.R")
  source("R/run_cf.R")
  library(rstan)
  
  # compile model
  model = stan_model(fixed_model_path)
  
  # calculate posterior mean 
  post_chains_sum = summarise_posterior_chains_DIC(posterior_chains = posterior_chains)
  
  # run model using posterior mean values 
  model_posts = replicate_rstan_fixed(model = model,
                                      posterior_sample_row = data.frame(t(post_chains_sum))[1,],
                                      scale_time_step = scale_time_step,
                                      n_pop_it = n_pop_it,
                                      n_recov_it =n_recov_it,
                                      index_M_it = index_M_it,
                                      index_A_it = index_A_it,
                                      index_O_it =  index_O_it,
                                      index_Al_it = index_Al_it,
                                      n_pop_veneto = n_pop_veneto,
                                      n_recov_veneto = n_recov_veneto,
                                      index_M_veneto = index_M_veneto,
                                      index_A_veneto = index_A_veneto,
                                      index_O_veneto = index_O_veneto,
                                      index_Al_veneto = index_Al_veneto,
                                      epsilon = epsilon,
                                      sigma = sigma,
                                      gamma = gamma,
                                      mu = mu ,
                                      phi_PCR = phi_PCR,
                                      phi_Ag =  phi_Ag,
                                      start_date_it =  start_date_it,
                                      end_date_it = end_date_it,
                                      time_intervention_it = time_intervention_it ,
                                      time_seed_alpha_it = time_seed_alpha_it,
                                      time_seed_M_it = time_seed_M_it ,
                                      time_vac_it = time_vac_it,
                                      start_date_veneto =  start_date_veneto,
                                      end_date_veneto = end_date_veneto,
                                      time_intervention_veneto = time_intervention_veneto,
                                      time_seed_alpha_veneto = time_seed_alpha_veneto,
                                      time_seed_M_veneto = time_seed_M_veneto,
                                      time_vac_veneto = time_vac_veneto,
                                      italy_testing_in_veneto = F)
  
  
  # Calculate model deviance 
  
  E_log_lik =  model_posts$sum_LL # mean post value 
  
  LL = mean(rowSums(log_lik)) 
  
  pDIC = 2 * (E_log_lik - mean(rowSums(log_lik)))  # sum across the data points to get total log likelihood then calculate the mean for all iterations 
  
  dic = -2 * E_log_lik  + 2 * pDIC
  
  return(c(dic, LL))
}
