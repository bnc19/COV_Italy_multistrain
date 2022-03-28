##################################################################################################################################
# Function to format posterior chains and genomic and epidemiological data to run multivariant models assuming fixed parameters.      
##################################################################################################################################

format_data_for_SEIQR_model = function(
                          daily_reported_incidence,
                          average_monthly_vaccination,
                          monthly_PCR,
                          monthly_Ag,
                          deaths,
                          daily_PCR ,
                          daily_Ag,
                          start_date,
                          end_date ,
                          time_intervention ,
                          time_seed_alpha,
                          time_seed_M             ) {
  
  

  # convert character date to date object 
  
  start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
  end_date_d =  as.Date.character(end_date, format = "%d-%m-%Y")
  
  # vector of all dates to run model over 
  daily_date = seq.Date(from = start_date_d, to = end_date_d ,  by = "days")

  # convert from date to index number in data
  
  time_switch1 =  which(daily_date == as.Date.character(time_intervention[1], format = "%d-%m-%Y"))
  time_switch2 =  which(daily_date == as.Date.character(time_intervention[2], format = "%d-%m-%Y"))
  seed_alpha =  which(daily_date == as.Date.character(time_seed_alpha, format = "%d-%m-%Y"))
  seed_M =  which(daily_date == as.Date.character(time_seed_M, format = "%d-%m-%Y"))
  
  
  # number of days to run model over 
  n_days = length(daily_date)

  n_months = round(n_days / 30)
  
  index_1st_month = c(which(as.numeric(format(daily_date, "%d")) == 1), (tail(n_days) + 1)) # index by 1st of each month and the final day of the last month
  
  
  # format monthly data 
  
  
  monthly_PCR_i = monthly_PCR[(length(monthly_PCR) - n_months + 1):length(monthly_PCR)]
  
  monthly_Ag_i = monthly_Ag[(length(monthly_Ag) - n_months + 1):length(monthly_Ag)]
  
  
  average_monthly_vaccination_i = average_monthly_vaccination[(length(average_monthly_vaccination) -
                                                             n_months + 1):length(average_monthly_vaccination)]
  
  
  # format daily data to match desired run time 
  
  
  daily_reported_incidence_i = daily_reported_incidence[(length(daily_reported_incidence) - n_days + 1):length(daily_reported_incidence)]
  
  daily_PCR_i = daily_PCR[(length(daily_PCR) - n_days + 1):length(daily_PCR)]
  
  daily_Ag_i = daily_Ag[(length(daily_Ag) - n_days + 1):length(daily_Ag)]
  
  deaths_i = deaths$cum_deaths[(length(deaths$cum_deaths) - n_days + 1):length(deaths$cum_deaths)]
  # 
  # 
  # average_daily_vaccination$date.of.vaccination = as.Date.character(average_daily_vaccination$date.of.vaccination, format = "%Y-%m-%d")
  # average_daily_vaccination2 = left_join(data.frame(date.of.vaccination = daily_date), average_daily_vaccination)
  # 
  # average_daily_vaccination2$prop_vac =  ifelse(is.na(average_daily_vaccination2$prop_vac),  0, average_daily_vaccination2$prop_vac  )
  # 
  # average_daily_vaccination_i = average_daily_vaccination2$prop_vac [(length(average_daily_vaccination2$prop_vac) -  n_days + 1):length(average_daily_vaccination2$prop_vac)]
  # 
  # 
  
  ### combine data for output
  
  external_data = data.frame(
    daily_date,
    deaths_i,
    daily_Ag_i,
    daily_PCR_i,
    daily_reported_incidence_i
  )
  model_times = data.frame(
    n_days,
    seed_M,
    seed_alpha,
    time_switch1,
    time_switch2,
    n_months
    
  )
  
  monthly_data = data.frame(
    monthly_PCR_i,
    monthly_Ag_i,
    average_monthly_vaccination_i
  )
  
out =  list(external_data,  model_times,monthly_data, index_1st_month)

  
  return(out)
  
  }
