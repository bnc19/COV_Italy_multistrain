run_SEIR_stan_Italy_sens = function(Al_data,
                               A_data  ,
                               M_data  ,
                               O_data  ,
                               n_seq,
                               n_difeq ,
                               n_pop ,
                               n_recov ,
                               index_M = NULL ,
                               index_A = NULL ,
                               index_O = NULL ,
                               index_Al = NULL ,
                               modelPath,
                               Location,
                               n_chains,
                               n_warmups,
                               n_iter ,
                               n_thin ,
                               pars ,
                               filePath ,
                               ini_1_SEIR,
                               average_daily_vaccination,
                               daily_reported_incidence,
                               average_daily_reported_incidence,
                               daily_PCR ,
                               daily_Ag,
                               monthly_PCR ,
                               monthly_Ag,
                               start_date,
                               end_date ,
                               time_intervention ,
                               time_seed_alpha,
                               time_seed_M,
                               time_vaccination,
                               sigma = 1 / 5.1,
                               gamma = 1 / 2.1 ,
                               phi_Ag ,
                               phi_PCR , 
                               prior_seed_mean ,
                               prior_seed_sd, 
                               scale_time_step) {

  ########## Load packages    ########## 
  
  library(dplyr)
  library(ggplot2)
  library(rstan)
  library(bayesplot)
  library(Hmisc)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  
  

  ##########  If index not given, no missing data   ########## 

  # index for GISAID data
  
  if (is.null(index_A)) {
    index_M = 1:length(M_data)
  }
  
  if (is.null(index_A)) {
    index_A = 1:length(A_data)
  }
  
  if (is.null(index_O)) {
    index_O = 1:length(O_data)
  }
  
  if (is.null(index_Al)) {
    index_Al = 1:length(Al_data)
  }
  

  ##########  data wrangling to work with dates   ########## 

  
  
  start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
  end_date_d = as.Date.character(end_date, format = "%d-%m-%Y")
  
  daily_date = seq.Date(from = start_date_d, to = end_date_d ,  by = "days")
  
  index_switch1 =  which(daily_date == as.Date.character(time_intervention[1], format = "%d-%m-%Y")) * scale_time_step
  index_switch2 =  which(daily_date == as.Date.character(time_intervention[2], format = "%d-%m-%Y")) * scale_time_step
  index_seed_alpha =  which(daily_date == as.Date.character(time_seed_alpha, format = "%d-%m-%Y")) * scale_time_step
  index_seed_M =  which(daily_date == as.Date.character(time_seed_M, format = "%d-%m-%Y")) * scale_time_step

  
  n_days = length(daily_date)
  n_months = round(n_days / 30)
  
  index_1st_month = c(which(as.numeric(format(daily_date, "%d")) == 1), (tail(n_days) + 1)) # index by 1st of each month and the final day of the last month
  

  
  # format daily data to length of seed plus data fitting 

  daily_PCR_i = daily_PCR[(length(daily_PCR) - n_days + 1):length(daily_PCR)] / scale_time_step
  daily_PCR_i2 = rep(daily_PCR_i, each = scale_time_step)
  
  daily_Ag_i = daily_Ag[(length(daily_Ag) - n_days + 1):length(daily_Ag)]  / scale_time_step
  daily_Ag_i2 = rep(daily_Ag_i, each = scale_time_step)
  
  
  ## format vaccination 
  
  time_vac =  which(daily_date == as.Date.character("26-12-2020", format = "%d-%m-%Y"))
  
  average_daily_vaccination_i = rep(0, (time_vac+length(average_daily_vaccination)))

  average_daily_vaccination_i[(time_vac +1):(time_vac  + length(average_daily_vaccination))] = average_daily_vaccination / scale_time_step
  
  average_daily_vaccination_i2 = rep(average_daily_vaccination_i[1:n_days], each = scale_time_step) 
  
  
  
  # # format monthly data to length of seed plus data fitting 
  
  average_daily_reported_incidence_i = average_daily_reported_incidence[(length(average_daily_reported_incidence)-n_months+1): length(average_daily_reported_incidence)]
  
  
  n_seq_i = n_seq[(length(n_seq) - n_months + 1):length(n_seq)]
  
  
  ### having extracted GISAID data from first time of fitting (first occurrence of A variant), this allows to index variant specific data, removing NA values 
  
  index_M_i = index_M - (length(M_data)- n_months)
  index_A_i = index_A - (length(A_data)- n_months)
  index_O_i = index_O - (length(O_data)- n_months)
  index_Al_i  = index_Al - (length(Al_data)- n_months)
  
  
  ##########  data as list for Stan model  ########## 

  
  model_data = list(
    n_data_A  =  length(index_A),
    n_data_M =  length(index_M),
    n_data_O =  length(index_O),
    n_data_Al =  length(index_Al),
    
    
    n_difeq = n_difeq,
    n_pop = n_pop,
    n_days = n_days,
    n_recov = n_recov,
    n_months = n_months,
    n_variants = 4,
    n_ts =( n_days*scale_time_step),
    
    
    n_seq_M = n_seq[index_M],
    n_seq_A = n_seq[index_A],
    n_seq_O = n_seq[index_O],
    n_seq_Al = n_seq[index_Al],
    
    n_reported_M = average_daily_reported_incidence[index_M],
    n_reported_A = average_daily_reported_incidence[index_A],
    n_reported_O = average_daily_reported_incidence[index_O],
    n_reported_Al = average_daily_reported_incidence[index_Al],
    
    y_M = round(M_data[index_M]),
    y_A = round(A_data[index_A]),
    y_O = round(O_data[index_O]),
    y_Al = round(Al_data[index_Al]),
    
    index_M = index_M_i,
    index_A = index_A_i,
    index_O = index_O_i,
    index_Al = index_Al_i,
    
    gamma = gamma /  scale_time_step,
    sigma = sigma / scale_time_step,
    phi_Ag =  phi_Ag,
    phi_PCR = phi_PCR , 
  
    
    Ag_daily  = round(daily_Ag_i2),
    PCR_daily = round(daily_PCR_i2),
    vac = average_daily_vaccination_i2,
    
    time_switch1 = index_switch1,
    time_switch2 =index_switch2,
    
    time_seed_M = index_seed_M,
    time_seed_alpha = index_seed_alpha,

    month_index = index_1st_month,
    
    seed_mean = prior_seed_mean ,
    seed_sd = prior_seed_sd,
    scale_time_step = scale_time_step
  )
  
  
  ##########  set up model   ########## 

  
  SEIR_model = stan_model(modelPath)
  
  

  ##########  run model    ########## 
  time.start <- Sys.time()
  
    SEIR__fit = sampling(
      SEIR_model,
      data = model_data,
      init = ini_1_SEIR,
      chains = n_chains,
      warmup = n_warmups,
      iter = n_iter,
      thin = n_thin, 
      control = list(
        adapt_delta = 0.90, 
        max_treedepth = 15
      )) 
    
    time.end <- Sys.time()
    time.end - time.start 
  
  SEIR__fit_summary = summary(SEIR__fit, pars = pars)$summary
  
  
  write.csv(SEIR__fit_summary, file = paste0(filePath, "/posterior.csv"))
  
  

  ########## diagnostics   ########## 

  SEIR_fit_posterior = as.array(SEIR__fit)
  color_scheme_set("viridis")
  
  # Markov chain traceplots
  
  param_trace = mcmc_trace(SEIR_fit_posterior, pars = c(pars))
  ggsave(param_trace, file = paste0(filePath, "/param_trace.jpg"), width = 20, height = 20)
  
  
  param_pairs = mcmc_pairs(SEIR_fit_posterior, pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho"   ))
  ggsave(param_pairs, file = paste0(filePath, "/param_pairs.jpg"))
  
  

  ##########  plot model and cases   ########## 
  
  # extract data from posterior
  SEIR_fit_posts =  rstan::extract(SEIR__fit)
  
  
  fit_SEIR_M = SEIR_fit_posts$reported_incidence[,, 1]
  fit_SEIR_A = SEIR_fit_posts$reported_incidence[,, 2]
  fit_SEIR_O = SEIR_fit_posts$reported_incidence[,, 3]
  fit_SEIR_Al = SEIR_fit_posts$reported_incidence[,, 4]
  

  # list all variants
  
  reportedList = list(fit_SEIR_M, fit_SEIR_A, fit_SEIR_O, fit_SEIR_Al)  # reported incidence, plot against data

  dataList = list(M_data[index_M], A_data[index_A], O_data[index_O], Al_data[index_Al])  # data
  indexList = list(index_M_i, index_A_i, index_O_i, index_Al_i)  # index data by month
  
  # empty lists
  Median = Lower = Upper  =  dataConf = dataDay = list()
 
  # calculate median and 95% CrI. Time = columns, iterations = rows.
  
  for (i in 1:length(dataList)) {
    # find the median of each column
    
    Median[[i]] = apply(reportedList[[i]], 2, median)
    Lower[[i]] = apply(reportedList[[i]], 2, quantile, probs = c(0.025))
    Upper[[i]] = apply(reportedList[[i]], 2, quantile, probs = c(0.975))
  }
  
  # account for seeding and aggregating from day to month
  
  monthly_date = daily_date[ which(as.numeric(format(daily_date, "%d")) == 1) ]
  

  # calculate binomial confidence intervals
  
  for (i in 1:length(dataList)) {
    dataConf[[i]] = data.frame(
      Date = monthly_date[indexList[[i]]],
      binconf(dataList[[i]], n_seq_i[indexList[[i]]]) * average_daily_reported_incidence_i[indexList[[i]]]
    )
    
  }
  
  dataDay = left_join(data.frame(Date = daily_date), dataConf[[1]])
  for (i in 2:length(dataConf)) {
    dataDay = left_join(dataDay, dataConf[[i]], by = "Date")
  }
  
  
  reportedIncidence = data.frame(dataDay,
                                 Median, Lower, Upper)
  
  
  colnames(reportedIncidence) = c(
    "Date",
    "M_Variant_M",
    "M_Variant_L",
    "M_Variant_U",
    "A_Variant_M",
    "A_Variant_L",
    "A_Variant_U",
    "O_Variant_M",
    "O_Variant_L",
    "O_Variant_U",
    "Al_Variant_M",
    "Al_Variant_L",
    "Al_Variant_U",
    "M_fit_M",
    "A_fit_M" ,
    "O_fit_M",
    "Al_fit_M",
    "M_fit_L",
    "A_fit_L" ,
    "O_fit_L" ,
    "Al_fit_L",
    "M_fit_U",
    "A_fit_U" ,
    "O_fit_U",
    "Al_fit_U"
  )
  
  

  
  ########## plot reported incidence against data  ########## 
  
  reportedIncidencePlot = ggplot(reportedIncidence, aes(x = Date, y = M_Variant_M)) +
    geom_point(shape = 19, size = 3, (aes(color = "M234I-A376T"))) +
    geom_errorbar(aes(ymin = M_Variant_L, ymax = M_Variant_U)) +
    geom_ribbon(aes(ymin = M_fit_L, ymax = M_fit_U),
                fill = "orange",
                alpha = 0.5) +
    geom_line(aes(y = M_fit_M, color = "M234I-A376T model"), size = 1) +
    geom_point(shape = 19, size = 3, (aes(y = A_Variant_M, color = "A220V"))) +
    geom_errorbar(color = "blue", aes(ymin = A_Variant_L, ymax = A_Variant_U)) +
    geom_ribbon(aes(ymin = A_fit_L, ymax = A_fit_U),
                fill = "sky blue",
                alpha = 0.5) +
    geom_line(aes(y = A_fit_M, color = "A220V model"), size = 1) +
    geom_point(shape = 19, size = 3, (aes(y = O_Variant_M, color = "Other"))) +
    geom_errorbar(color = "red", aes(ymin = O_Variant_L, ymax = O_Variant_U)) +
    geom_ribbon(aes(ymin = O_fit_L, ymax = O_fit_U),
                fill = "pink",
                alpha = 0.5) +
    geom_line(aes(y = O_fit_M, color = "Other model"), size = 1) +
    geom_point(shape = 19, size = 3, (aes(y = Al_Variant_M, color = "Alpha"))) +
    geom_errorbar(color = "mediumpurple4", aes(ymin = Al_Variant_L, ymax = Al_Variant_U)) +
    geom_ribbon(aes(ymin = Al_fit_L, ymax = Al_fit_U),
                fill = "lavender",
                alpha = 0.5) +
    geom_line(aes(y = Al_fit_M, color = "Alpha model"), size = 1) +
    scale_colour_manual(
      name = '',
      values = c(
        'M234I-A376T' = 'black',
        'M234I-A376T model' = 'orange',
        "A220V" = "blue",
        'A220V model' = 'sky blue',
        'Other model' = 'pink',
        'Other' = 'red',
        'Alpha model' = 'lavender',
        'Alpha' = 'mediumpurple4'
      )
    ) +
    labs(x = "Date", y = paste0("Incidence")) +
    theme_bw() + theme(text = element_text(size = 16)) +
    scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
    ggtitle(Location)
  
  
  
  
  ########## compile data  ########## 
  
 
  reportedIncidence[, -1] =  round(reportedIncidence[, -1])
  
  
  beta_chain_M =  SEIR_fit_posts$beta[,1]
  beta_chain_A =  SEIR_fit_posts$beta[,2]
  beta_chain_O =  SEIR_fit_posts$beta[,3]
  beta_chain_Al =  SEIR_fit_posts$beta[,4]
  
  rho_chain =  SEIR_fit_posts$rho
  omega1_chain =  SEIR_fit_posts$omega[,1]
  omega2_chain =  SEIR_fit_posts$omega[,2]
  
  I0_chain_M =  SEIR_fit_posts$I0[,1]
  I0_chain_A =  SEIR_fit_posts$I0[,2]
  I0_chain_O =  SEIR_fit_posts$I0[,3]
  I0_chain_Al =  SEIR_fit_posts$I0[,4]
  k_chain = SEIR_fit_posts$k
  
  
  
  posteriorChains = data.frame( beta_chain_A, beta_chain_M, beta_chain_O, beta_chain_Al,
                               rho_chain, omega1_chain, omega2_chain,k_chain, I0_chain_M, I0_chain_A, I0_chain_O, I0_chain_Al )
  
  ########## Save plot and Data  ########## 
  ggsave(
    reportedIncidencePlot,
    file = paste0(filePath, "/reportedIncidencePlot.jpg"),
    height = 20,
    width = 30,
    unit = "cm",
    dpi = 720
  )

  write.csv(posteriorChains, file =  paste0(filePath, "/posteriorChains.csv"))
 
  write.csv(reportedIncidence,
            file = paste0(filePath,  "/reportedIncidence.csv"))

  return(SEIR__fit)
  
}
