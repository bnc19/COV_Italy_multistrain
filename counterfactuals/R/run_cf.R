# run multivariant models assuming fixed parameters ----------------------------
# returns a matrix -------------------------------------------------------------

replicate_rstan_fixed = function(model ,
                                 posterior_sample_row,
                                 scale_time_step = 2,
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
                                 italy_testing_in_veneto = F)
                          
                                {
  
  
source("R/format_stan_cf_data.R")
  

 data =  format_stan_cf_data(
    scale_time_step = scale_time_step,
    n_pop_it = n_pop_it,
    n_recov_it = n_recov_it,
    index_M_it = index_M_it,
    index_A_it = index_A_it,
    index_O_it =  index_O_it,
    index_Al_it =index_Al_it,
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
    start_date_it = start_date_it,
    end_date_it =end_date_it,
    time_intervention_it = time_intervention_it ,
    time_seed_alpha_it = time_seed_alpha_it,
    time_seed_M_it = time_seed_M_it ,
    time_vac_it = time_vac_it,
    start_date_veneto = start_date_veneto,
    end_date_veneto = end_date_veneto,
    time_intervention_veneto = time_intervention_veneto,
    time_seed_alpha_veneto = time_seed_alpha_veneto,
    time_seed_M_veneto =time_seed_M_veneto,
    time_vac_veneto = time_vac_veneto,
    posterior_sample_row = posterior_sample_row,
    italy_testing_in_veneto = italy_testing_in_veneto
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
                               location,
                               model = "cf"){
  
  if(location== "Veneto"){
    rep_M = posts$daily_incidence_ven[,, 1]
    rep_A = posts$daily_incidence_ven[,, 2]
    rep_O = posts$daily_incidence_ven[,, 3]
    rep_Al = posts$daily_incidence_ven[,, 4]
    
    true_M = posts$daily_incidence_ven_t[,, 1]
    true_A = posts$daily_incidence_ven_t[,, 2]
    true_O = posts$daily_incidence_ven_t[,, 3]
    true_Al = posts$daily_incidence_ven_t[,, 4]
    
    ratio_M = posts$ratio_ven[,, 1]
    ratio_A = posts$ratio_ven[,, 2]
    ratio_O = posts$ratio_ven[,, 3]
    ratio_Al = posts$ratio_ven[,, 4]
    
    
    
    R0_M = posts$R0_ven_daily[,, 1]
    R0_A = posts$R0_ven_daily[,, 2]
    R0_O = posts$R0_ven_daily[,, 3]
    R0_Al = posts$R0_ven_daily[,, 4]
    
    pMO = posts$pPCR_ven[,]
    
  } else if (location == "Italy") {
    rep_M = posts$daily_incidence_it[,, 1]
    rep_A = posts$daily_incidence_it[,, 2]
    rep_O = posts$daily_incidence_it[,, 3]
    rep_Al = posts$daily_incidence_it[,, 4]
    
    true_M = posts$daily_incidence_it_t[,, 1]
    true_A = posts$daily_incidence_it_t[,, 2]
    true_O = posts$daily_incidence_it_t[,, 3]
    true_Al = posts$daily_incidence_it_t[,, 4]
    
    ratio_M = posts$ratio_it[,, 1]
    ratio_A = posts$ratio_it[,, 2]
    ratio_O = posts$ratio_it[,, 3]
    ratio_Al = posts$ratio_it[,, 4]
    
    
    R0_M = posts$R0_it_daily[,, 1]
    R0_A = posts$R0_it_daily[,, 2]
    R0_O = posts$R0_it_daily[,, 3]
    R0_Al = posts$R0_it_daily[,, 4]
    
    pMO = posts$pPCR_it[,]
  }
  
  if(model == "baseline"){
  out = as.matrix(data.frame(rep_M, rep_A, rep_O, rep_Al,
                             true_M, true_A, true_O, true_Al, 
                             ratio_M,ratio_A,ratio_O,ratio_Al,
                             R0_M, R0_A, R0_O, R0_Al, pMO))
  } else if (model == "cf"){
    out = as.matrix(data.frame(rep_M, rep_A, rep_O, rep_Al,
                               true_M, true_A, true_O, true_Al))
  } else if (model == "SA"){
    out = as.matrix(data.frame(rep_M, rep_A, rep_O, rep_Al))
  }
  
  return(out)
}


################ change reporting rate to 1 e.g. mass antigen testing ################ 


change_rho = function(data.frame){
  data.frame$rho_chain = 1
  return(data.frame)
}




# summarise posterior across iterations and add date ---------------------------
# returns a data frame of reported and true inc over time ----------------------

summarise_results = function(posterior_results,
                             start_date,
                             end_date,
                             S0,
                             no.col = 8){
  
  start = as.Date.character(start_date, format = "%d-%m-%Y")
  end = as.Date.character(end_date, format = "%d-%m-%Y")
  
  all_dates = 
    seq.Date(from = start, to = end,  by = "days")
  
  out =  posterior_results %>% as.data.frame.table() %>%
    rename(time = Var1, variant = Var2,ni = Var3, value = Freq) %>%
    filter(!grepl("ratio",variant )) %>% 
    filter(!grepl("R0",variant )) %>%
    filter(!grepl("pMO", variant)) %>%  
    dplyr::mutate(ni = as.numeric(ni),
                  time = as.numeric(time),
                  value =  value / S0 * 100000) %>%
    group_by(time,variant) %>%
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)) %>% 
    ungroup( ) %>%  
    mutate(Date = rep(all_dates, each = no.col)) %>% 
    separate(variant, into = c("output", "variant")) %>% 
    pivot_wider(id_cols = c("Date", "variant"), names_from = output,
                values_from = c(lower,mean,upper))
  
  return(out)
  
}


# calculate ratio of incidence reported incidence ------------------------------
# returns a data frame of reported inc, true inc, ratio for each variant -------
# aggregated over the study period ---------------------------------------------

calculate_ratio_reported = function(posterior_results,
                                    S0){
  
  out = posterior_results %>% as.data.frame.table() %>%
    rename(time = Var1, variant = Var2,ni = Var3, value = Freq) %>%
    filter(grepl("ratio",variant )) %>%  
    dplyr::mutate(ni = as.numeric(ni),
                  time = as.numeric(time),
                  value = value ) %>%
    group_by(variant, ni) %>% 
    summarise(value = mean(value, na.rm=T)) %>%   # calculate cumulative incidence 
    group_by(variant) %>%
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) %>%  
    ungroup( )
  
  
  out2 = posterior_results %>% as.data.frame.table() %>%
    rename(time = Var1, variant = Var2,ni = Var3, value = Freq) %>%
    filter(!grepl("ratio",variant )) %>%  
    dplyr::mutate(ni = as.numeric(ni),
                  time = as.numeric(time),
                  value = value ) %>%
    group_by(variant, ni) %>% 
    summarise(value = sum(value, na.rm=T) / S0 * 100) %>%   # calculate cumulative incidence 
    pivot_wider(id_cols = ni, names_from = variant, values_from = value ) %>%  
    mutate(true_tot = true_M+true_A+true_O+true_Al, 
           rep_tot = rep_M+rep_A+rep_O+rep_Al) %>% 
    pivot_longer(cols= -ni, names_to = "variant") %>% 
    group_by(variant) %>%
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) %>%  
    ungroup( ) %>% 
    bind_rows(out)
  
  return(out2)
}



# 
# # calculate R0 of each variant over time ---------------------------------------
# 
# 
# calculate_R0= function(posterior_results,
#                        start_date,
#                        end_date
#                        ){
#   
#   start = as.Date.character(start_date, format = "%d-%m-%Y")
#   end = as.Date.character(end_date, format = "%d-%m-%Y")
#   
#   all_dates = 
#     seq.Date(from = start, to = end,  by = "days")
#   
#   out = posterior_results %>% as.data.frame.table() %>%
#     rename(time = Var1, variant = Var2,ni = Var3, value = Freq) %>%
#     filter(grepl("R0",variant )) %>%  
#     dplyr::mutate(ni = as.numeric(ni),
#                   time = as.numeric(time),
#                   value = value ) %>%
#     group_by(variant,time) %>%
#     summarise(
#       lower = quantile(value, 0.025),
#       mean = mean(value),
#       upper = quantile(value, 0.975)
#     ) %>%  
#     ungroup( ) %>%  
#     mutate(Date = rep(all_dates,  4)) 
#   
# 
#   
#   return(out)
# }
# 
# # 
# # 
# # # calculate probability of takking molecular diagnostic test over time ---------
# # 
# # 
# # calculate_pMO = function(posterior_results,
# #                        start_date,
# #                        end_date
# # ){
# #   
# #   start = as.Date.character(start_date, format = "%d-%m-%Y")
# #   end = as.Date.character(end_date, format = "%d-%m-%Y")
# #   
# #   all_dates = 
# #     seq.Date(from = start, to = end,  by = "days")
# #   
# #   out = posterior_results %>% as.data.frame.table() %>%
# #     rename(time = Var1, variant = Var2,ni = Var3, value = Freq) %>%
# #     filter(variant == "pMO") %>%  
# #     dplyr::mutate(ni = as.numeric(ni),
# #                   time = as.numeric(time),
# #                   value = value ) %>%
# #     group_by(time) %>%
# #     summarise(
# #       lower = quantile(value, 0.025),
# #       mean = mean(value),
# #       upper = quantile(value, 0.975)
# #     ) %>%  
# #     ungroup( ) %>%  
# #     mutate(Date = all_dates) 
# #   
# #   
# #   
# #   return(out)
# # }
