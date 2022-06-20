# run multivariant models assuming fixed parameters ----------------------------
# returns a matrix -------------------------------------------------------------

replicate_rstan_fixed = function(model_path ,
                                 posterior_sample_row,
                                 scale_time_step = 5)
                                {
  
  
source("R/format_stan_cf_data.R")
  
  model = stan_model(model_path)

 data =  format_stan_cf_data(
    scale_time_step = scale_time_step,
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
    posterior_sample_row = posterior_sample_row
  )
    
fit = sampling(
    model,
    data = data,
    chains = 1,
    iter = 1,
    algorithm = "Fixed_param"
    
  )
  
  
posts =  rstan::extract(fit)


return(posts)
  
}

# extract reported and true incidence from model fit ---------------------------
# returns a matrix -------------------------------------------------------------

extract_fit_results = function(posts,
                               location){
  
  if(location== "Veneto"){
    rep_M = posts$daily_incidence_ven[,, 1]
    rep_A = posts$daily_incidence_ven[,, 2]
    rep_O = posts$daily_incidence_ven[,, 3]
    rep_Al = posts$daily_incidence_ven[,, 4]
    
    true_M = posts$daily_incidence_ven_t[,, 1]
    true_A = posts$daily_incidence_ven_t[,, 2]
    true_O = posts$daily_incidence_ven_t[,, 3]
    true_Al = posts$daily_incidence_ven_t[,, 4]
  } else if (location == "Italy") {
    rep_M = posts$daily_incidence_it[,, 1]
    rep_A = posts$daily_incidence_it[,, 2]
    rep_O = posts$daily_incidence_it[,, 3]
    rep_Al = posts$daily_incidence_it[,, 4]
    
    true_M = posts$daily_incidence_it_t[,, 1]
    true_A = posts$daily_incidence_it_t[,, 2]
    true_O = posts$daily_incidence_it_t[,, 3]
    true_Al = posts$daily_incidence_it_t[,, 4]
  }
  
  out = as.matrix(data.frame(rep_M, rep_A, rep_O, rep_Al,
                             true_M, true_A, true_O, true_Al))
  
  
  return(out)
}


################ change reporting rate to 1 e.g. mass antigen testing ################ 


change_rho = function(data.frame){
  data.frame$rho_chain = 1
  return(data.frame)
}
