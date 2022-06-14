format_stan_data= function(
  n_pop_it ,
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
  scale_time_step
){
  

  # Import Italy data ------------------------------------------------------------
  A_data_it = read.csv("data/Dataset_italy_A_v5.csv")$Freq 
  M_data_it  = read.csv("data/Dataset_italy_M_v5.csv")$Freq 
  O_data_it  = read.csv("data/Dataset_italy_O_v1.csv")$Freq 
  Al_data_it  = read.csv("data/Dataset_italy_Alpha_v1.csv")$Freq 
  n_seq_it  = read.csv("data/Dataset_italy_A_v5.csv")$TotSeq
  avg_daily_rep_inc_it  = read.csv("data/Dataset_italy_A_v5.csv")$new_reported_cases_daily
  daily_PCR_it  = round(read.csv("data/Italy_daily_test_data.csv")$pcr_daily)
  daily_Ag_it = round(read.csv("data/Italy_daily_test_data.csv")$antigen_daily)
  average_daily_vaccination_it  =  read.csv("data/data_vac_italy_day.csv")$prop_vac
  
  # Import Veneto data -----------------------------------------------------------
  A_data_veneto = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq 
  M_data_veneto  = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq 
  O_data_veneto  = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq 
  Al_data_veneto  = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq 
  n_seq_veneto  = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq
  avg_daily_rep_inc_veneto  = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily
  daily_PCR_veneto  = round(read.csv("data/Veneto_daily_test_data.csv")$pcr_daily)
  daily_Ag_veneto = round(read.csv("data/Veneto_daily_test_data.csv")$antigen_daily)
  average_daily_vaccination_veneto  =  read.csv("data/data_vac_veneto_day.csv")$prop_vac
  
  
  # format dates and find index Italy --------------------------------------------
  
  
  start_date_d_it = as.Date.character(start_date_it, format = "%d-%m-%Y")
  end_date_d_it = as.Date.character(end_date_it, format = "%d-%m-%Y")
  
  all_dates_it = seq.Date(from = start_date_d_it, to = end_date_d_it ,  by = "days")
  
  index_switch1_it =  
    which(all_dates_it == as.Date.character(time_intervention_it[1], format = "%d-%m-%Y"))
  index_switch2_it =  
    which(all_dates_it == as.Date.character(time_intervention_it[2], format = "%d-%m-%Y"))
  index_seed_alpha_it =  
    which(all_dates_it == as.Date.character(time_seed_alpha_it, format = "%d-%m-%Y"))
  index_seed_M_it =  
    which(all_dates_it == as.Date.character(time_seed_M_it, format = "%d-%m-%Y"))
  time_vac_it =  
    which(all_dates_it == as.Date.character(time_vac_it, format = "%d-%m-%Y"))
  
  n_days_it = length(all_dates_it)
  n_months_it = round(n_days_it / 30)
  
  index_1st_month_it = 
    c(which(as.numeric(format(all_dates_it, "%d")) == 1), (tail(n_days_it) + 1)) 
  # index by 1st of each month and the final day of the last month
  
  
  # extract correct dates of data 
  
  
  daily_PCR_i_it = 
    daily_PCR_it[(length(daily_PCR_it) - n_days_it + 1):length(daily_PCR_it)] 
  daily_Ag_i_it = 
    daily_Ag_it[(length(daily_Ag_it) - n_days_it + 1):length(daily_Ag_it)]  
  
  average_daily_vaccination_i_it = rep(0, n_days_it)
  average_daily_vaccination_i_it[time_vac_it:n_days_it] = 
    average_daily_vaccination_it[1:length(time_vac_it:n_days_it)]
  
  n_seq_i_it = n_seq_it[(length(n_seq_it) - n_months_it + 1):length(n_seq_it)]
  avg_daily_rep_inc_i_it = 
    avg_daily_rep_inc_it[(length(avg_daily_rep_inc_it)-n_months_it+1): 
                           length(avg_daily_rep_inc_it)]
  
  # having extracted GISAID data from first date of fitting 
  # index variant specific data, removing NA values 
  
  index_M_i_it  = index_M_it - (length(M_data_it)- n_months_it)
  index_A_i_it  = index_A_it - (length(A_data_it)- n_months_it)
  index_O_i_it  = index_O_it - (length(O_data_it)- n_months_it)
  index_Al_i_it = index_Al_it - (length(Al_data_it)- n_months_it)
  
  
  # estimate variant specific reported incidence 
  n_reported_M_it  = 
    round(avg_daily_rep_inc_it * M_data_it / n_seq_it)[index_M_it]
  n_reported_A_it  = 
    round(avg_daily_rep_inc_it * A_data_it / n_seq_it)[index_A_it]
  n_reported_O_it  = 
    round(avg_daily_rep_inc_it * O_data_it / n_seq_it)[index_O_it]
  n_reported_Al_it = 
    round(avg_daily_rep_inc_it * Al_data_it / n_seq_it)[index_Al_it]
  
  # n data 
  
  n_data_it = c(length(index_M_it),length(index_A_it), 
                length(index_O_it),length(index_Al_it))
  
  
  
  
  # format dates and find index Veneto -------------------------------------------
  
  
  start_date_d_veneto = as.Date.character(start_date_veneto, format = "%d-%m-%Y")
  end_date_d_veneto = as.Date.character(end_date_veneto, format = "%d-%m-%Y")
  
  all_dates_veneto = 
    seq.Date(from = start_date_d_veneto, to = end_date_d_veneto,  by = "days")
  
  index_switch1_veneto =  
    which(all_dates_veneto == as.Date.character(time_intervention_veneto[1], format = "%d-%m-%Y"))
  index_switch2_veneto =  
    which(all_dates_veneto == as.Date.character(time_intervention_veneto[2], format = "%d-%m-%Y"))
  index_seed_alpha_veneto =  
    which(all_dates_veneto == as.Date.character(time_seed_alpha_veneto, format = "%d-%m-%Y"))
  index_seed_M_veneto =  
    which(all_dates_veneto == as.Date.character(time_seed_M_veneto, format = "%d-%m-%Y"))
  time_vac_veneto =  
    which(all_dates_veneto == as.Date.character(time_vac_veneto, format = "%d-%m-%Y"))
  
  n_days_veneto = length(all_dates_veneto)
  n_months_veneto = round(n_days_veneto/ 30)
  
  index_1st_month_veneto =
    c(which(as.numeric(format(all_dates_veneto, "%d")) == 1), (tail(n_days_veneto) + 1)) 
  # index by 1st of each month and the final day of the last month
  
  
  # extract correct dates of data ------------------------------------------------
  
  
  daily_PCR_i_veneto = 
    daily_PCR_veneto[(length(daily_PCR_veneto) - n_days_veneto+ 1):length(daily_PCR_veneto)] 
  daily_Ag_i_veneto = 
    daily_Ag_veneto[(length(daily_Ag_veneto) - n_days_veneto+ 1):length(daily_Ag_veneto)]  
  
  average_daily_vaccination_i_veneto = rep(0, n_days_veneto)
  average_daily_vaccination_i_veneto[time_vac_veneto:n_days_veneto] = average_daily_vaccination_veneto[1:length(time_vac_veneto:n_days_veneto)]
  
  n_seq_i_veneto = 
    n_seq_veneto[(length(n_seq_veneto) - n_months_veneto+ 1):length(n_seq_veneto)]
  
  avg_daily_rep_inc_i_veneto = 
    avg_daily_rep_inc_veneto[(length(avg_daily_rep_inc_veneto)-n_months_veneto+1): 
                               length(avg_daily_rep_inc_veneto)]
  
  # having extracted GISAID data from first date of fitting 
  # index variant specific data, removing NA values 
  
  index_M_i_veneto  = index_M_veneto- (length(M_data_veneto)- n_months_veneto)
  index_A_i_veneto  = index_A_veneto- (length(A_data_veneto)- n_months_veneto)
  index_O_i_veneto  = index_O_veneto- (length(O_data_veneto)- n_months_veneto)
  index_Al_i_veneto = index_Al_veneto- (length(Al_data_veneto)- n_months_veneto)
  
  
  # estimate variant specific reported incidence 
  
  
  
  n_reported_M_veneto  =
    round(avg_daily_rep_inc_veneto* M_data_veneto/ n_seq_veneto)[index_M_veneto]
  n_reported_A_veneto  = 
    round(avg_daily_rep_inc_veneto* A_data_veneto/ n_seq_veneto)[index_A_veneto]
  n_reported_O_veneto  = 
    round(avg_daily_rep_inc_veneto* O_data_veneto/ n_seq_veneto)[index_O_veneto]
  n_reported_Al_veneto = 
    round(avg_daily_rep_inc_veneto* Al_data_veneto/ n_seq_veneto)[index_Al_veneto]
  
  # n data 
  
  n_data_veneto = c(length(index_M_veneto),length(index_A_veneto), 
                    length(index_O_veneto),length(index_Al_veneto))
  
  
  # List of data for Rstan -------------------------------------------------------
  
  model_data_real = list(
    # Shared data 
    scale_time_step = scale_time_step,
    gamma = gamma ,
    sigma = sigma ,
    mu = mu, 
    epsilon = epsilon, 
    phi_Ag =  phi_Ag,
    phi_PCR = phi_PCR , 
    n_var = 4,
    
    # Italy data 
    
    n_data_it = n_data_it, 
    
    index_M_it  = index_M_i_it,
    index_A_it  = index_A_i_it,
    index_O_it  = index_O_i_it,
    index_Al_it = index_Al_i_it,
    
    n_pop_it    = n_pop_it,
    n_days_it   = length(all_dates_it),
    n_recov_it  = n_recov_it,
    n_months_it = n_months_it,
    n_ts_it     = (length(all_dates_it) * scale_time_step),
    
    y_M_it  = round(n_reported_M_it  / (n_pop_it-n_recov_it) * 100000),
    y_A_it  = round(n_reported_A_it  / (n_pop_it-n_recov_it) * 100000),
    y_O_it  = round(n_reported_O_it  /  (n_pop_it-n_recov_it) * 100000),
    y_Al_it = round(n_reported_Al_it /  (n_pop_it-n_recov_it) * 100000),
    
    Ag_daily_it  = round(rep(daily_Ag_i_it            / scale_time_step, each = scale_time_step)),
    PCR_daily_it = round(rep(daily_PCR_i_it           / scale_time_step, each = scale_time_step)),
    vac_it       = rep(average_daily_vaccination_i_it / scale_time_step, each = scale_time_step),
    
    time_switch1_it    = index_switch1_it    * scale_time_step,
    time_switch2_it    = index_switch2_it    * scale_time_step,
    time_seed_M_it     = index_seed_M_it     * scale_time_step,
    time_seed_alpha_it = index_seed_alpha_it * scale_time_step,
    
    month_index_it = index_1st_month_it ,
    
    # Veneto data 
    
    n_data_ven   = n_data_veneto, 
    
    index_M_ven  = index_M_i_veneto,
    index_A_ven  = index_A_i_veneto,
    index_O_ven  = index_O_i_veneto,
    index_Al_ven = index_Al_i_veneto,
    
    n_pop_ven    = n_pop_veneto,
    n_days_ven   = length(all_dates_veneto),
    n_recov_ven  = n_recov_veneto,
    n_months_ven = n_months_veneto,
    n_ts_ven     = (length(all_dates_veneto) * scale_time_step),
    
    y_M_ven  = round(n_reported_M_veneto  / (n_pop_veneto-n_recov_veneto) * 100000),
    y_A_ven  = round(n_reported_A_veneto  / (n_pop_veneto-n_recov_veneto) * 100000),
    y_O_ven  = round(n_reported_O_veneto  / (n_pop_veneto-n_recov_veneto) * 100000),
    y_Al_ven = round(n_reported_Al_veneto / (n_pop_veneto-n_recov_veneto) * 100000),
    
    Ag_daily_ven  = round(rep(daily_Ag_i_veneto            / scale_time_step, each = scale_time_step)),
    PCR_daily_ven = round(rep(daily_PCR_i_veneto           / scale_time_step, each = scale_time_step)),
    vac_ven       = rep(average_daily_vaccination_i_veneto / scale_time_step, each = scale_time_step), 
    
    time_switch1_ven    = index_switch1_veneto    * scale_time_step,
    time_switch2_ven    = index_switch2_veneto    * scale_time_step,
    time_seed_M_ven     = index_seed_M_veneto     * scale_time_step,
    time_seed_alpha_ven = index_seed_alpha_veneto * scale_time_step,
    
    month_index_ven = index_1st_month_veneto
  )
  
  return(model_data_real)
}
