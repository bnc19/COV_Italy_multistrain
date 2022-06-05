rm(list = ls())
setwd("Q:/COV_Italy_multistrain/model_fitting")
modelPath = "model/new_neg_bin_sens_estimate_new_likelihood_presymp_both.stan"

# Load packages ----------------------------------------------------------------
library(bayesplot)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(Hmisc)
library(cowplot)

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

# Define variables Italy -------------------------------------------------------


# initial conditions
n_pop_it = 59257566 - 4847026
n_recov_it = 1482377 - 93401

# index data Italy 
index_M_it = 6:14
index_A_it = 3:14
index_O_it =  3:9
index_Al_it = 9:14

# Define variables Veneto ------------------------------------------------------

# initial conditions
n_pop_veneto= 4847026
n_recov_veneto= 93401

# index data Veneto 
index_M_veneto = 8:14
index_A_veneto = c(5,7:14)
index_O_veneto = c(5,7:9)
index_Al_veneto = 9:14

# Shared parameters ------------------------------------------------------------


#5.1 days = incubation period (no symptoms) 
epsilon = 1/ 2.8 ## mean of values given in Vo table S5 assuming R0 = 2.7

# 1/ sigma = 5.1 - 2.8 

sigma = 1 / (5.1 - 2.8 )
gamma = 1 / 2.1 
mu = 0.59 # probability symptomatic 

phi_PCR = 0.920
phi_Ag =  0.643 # change antigen sensitivity 


# Dates Italy ------------------------------------------------------------------
start_date_it =  "01-05-2020" 
end_date_it = "31-05-2021"
time_intervention_it= c("15-11-2020" , "15-03-2021") 

time_seed_alpha_it = "01-11-2020"
time_seed_M_it= "01-08-2020" 
time_vac_it = "27-12-2020" # date first vaccine (first datapoint)


# Dates Veneto -----------------------------------------------------------------
start_date_veneto=  "01-07-2020" 
end_date_veneto= "31-05-2021"
time_intervention_veneto= c("15-11-2020" , "15-03-2021") 

time_seed_alpha_veneto= "01-11-2020"
time_seed_M_veneto= "01-10-2020"
time_vac_veneto= "27-12-2020" # date first vaccine (first datapoint)


# Stan parameters -------------------------------------------------------------- 
n_chains =4
n_warmups =2000
n_iter = 4000
n_thin = 1
pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", 
         "rho_it" , "rho_ven" ,
         "omega[1]"  , "omega[2]" #, "omega[3]"  , "omega[4]"
         , "I0_it[1]", "I0_it[2]", "I0_it[3]", "I0_it[4]",
         "I0_ven[1]", "I0_ven[2]", "I0_ven[3]", "I0_ven[4]",
         "R0_it", "R0_ven")

ini_1_SEIR = function(){
  list(  beta = replicate(4,runif(1,0,3)),
         I0_it = replicate(4, runif(1, 1,20)),
         I0_ven = replicate(4, runif(1, 1,20)),
         omega = replicate(2,runif(1,0.2,0.8)),
         rho_it = runif(1,0.2,0.8),
         rho_ven = runif(1,0.2,0.8),
         k = runif(1,0.01,2)
  )}



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
  
  gamma = gamma ,
  sigma = sigma ,
  mu = mu, 
  epsilon = epsilon, 
  phi_Ag =  phi_Ag,
  phi_PCR = phi_PCR , 
  n_var = 4,
  
  # Italy data 
  
  n_data_it = n_data_it, 
  index_M_it = index_M_i_it,
  index_A_it = index_A_i_it,
  index_O_it = index_O_i_it,
  index_Al_it = index_Al_i_it,
  
  n_pop_it = n_pop_it,
  n_days_it = length(all_dates_it),
  n_recov_it = n_recov_it,
  n_months_it = n_months_it,
  n_ts_it =length(all_dates_it),
  
  y_M_it = round(n_reported_M_it / (n_pop_it-n_recov_it) * 100000),
  y_A_it = round(n_reported_A_it / (n_pop_it-n_recov_it) * 100000),
  y_O_it = round(n_reported_O_it /  (n_pop_it-n_recov_it) * 100000),
  y_Al_it = round(n_reported_Al_it/  (n_pop_it-n_recov_it) * 100000),
  
  Ag_daily_it  = daily_Ag_i_it,
  PCR_daily_it = daily_PCR_i_it,
  vac_it = average_daily_vaccination_i_it,
  
  time_switch1_it = index_switch1_it,
  time_switch2_it =index_switch2_it,
  
  time_seed_M_it = index_seed_M_it,
  time_seed_alpha_it = index_seed_alpha_it,
  
  month_index_it = index_1st_month_it,
  
  # Veneto data 
  n_data_ven= n_data_veneto, 
  index_M_ven= index_M_i_veneto,
  index_A_ven= index_A_i_veneto,
  index_O_ven= index_O_i_veneto,
  index_Al_ven= index_Al_i_veneto,
  
  n_pop_ven= n_pop_veneto,
  n_days_ven= length(all_dates_veneto),
  n_recov_ven= n_recov_veneto,
  n_months_ven= n_months_veneto,
  n_ts_ven=length(all_dates_veneto),
  
  y_M_ven= round(n_reported_M_veneto/ (n_pop_veneto-n_recov_veneto) * 100000),
  y_A_ven= round(n_reported_A_veneto/ (n_pop_veneto-n_recov_veneto) * 100000),
  y_O_ven= round(n_reported_O_veneto/  (n_pop_veneto-n_recov_veneto) * 100000),
  y_Al_ven= round(n_reported_Al_veneto/  (n_pop_veneto-n_recov_veneto) * 100000),
  
  Ag_daily_ven = daily_Ag_i_veneto,
  PCR_daily_ven= daily_PCR_i_veneto,
  vac_ven= average_daily_vaccination_i_veneto,
  
  time_switch1_ven= index_switch1_veneto,
  time_switch2_ven=index_switch2_veneto,
  
  time_seed_M_ven= index_seed_M_veneto,
  time_seed_alpha_ven= index_seed_alpha_veneto,
  
  month_index_ven= index_1st_month_veneto
)


# Run model --------------------------------------------------------------------


SEIR_model = stan_model(paste(modelPath))


time.start = Sys.time()

SEIR__fit = sampling(
  SEIR_model,
  data = model_data_real,
  init = ini_1_SEIR,
  control = list(
    adapt_delta = 0.85, 
    max_treedepth = 12
  )
) 

time.end = Sys.time()
time.end - time.start 


# Diagnostics ------------------------------------------------------------------

SEIR__fit_summary = summary(SEIR__fit, pars = pars)$summary
write.csv(SEIR__fit_summary, file="SEIR__fit_summary_10_L_location_omega.csv")

SEIR_fit_posterior = as.array(SEIR__fit)

# Markov chain traceplots

mcmc_trace(SEIR_fit_posterior, pars =pars)
mcmc_pairs(SEIR_fit_posterior, 
           pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho[1"   ))
mcmc_pairs(SEIR_fit_posterior, 
           pars = c( "beta[1]" , "I0[1]" , "rho", "inv_epsilon", "k", "omega[1]"   ))





#-------------------------------------------------------------------------------

# extract posterior estimates 
SEIR__fit_ext = rstan::extract(SEIR__fit)

varaints = c("M234I-A376T", "A220V", "Other", "Alpha")



# Plot model Italy -------------------------------------------------------------

dataConf_it = dataDay_it  = list()
# account for seeding and aggregating from day to month

monthly_date_it  = all_dates_it [ which(as.numeric(format(all_dates_it , "%d")) == 1) ]

dataList_it  = list(M_data_it[index_M_it ], A_data_it[index_A_it ], 
                    O_data_it[index_O_it ], Al_data_it[index_Al_it ])  # data
indexList_it  = list(index_M_i_it , index_A_i_it , index_O_i_it , index_Al_i_it )  # index data by month

# calculate binomial confidence intervals

for (i in 1:length(dataList_it )) {
  dataConf_it [[i]] = data.frame(
    Date = monthly_date_it [indexList_it [[i]]],
    binconf(dataList_it [[i]], n_seq_i_it[indexList_it [[i]]]) * 
      (avg_daily_rep_inc_i_it [indexList_it [[i]]] / (n_pop_it  -n_recov_it ) * 100000) ,
    variant = varaints[i]
  )
  
}

dataConfDf_it  = bind_rows(dataConf_it )
dataDay_it  = left_join(data.frame(Date = all_dates_it ),dataConfDf_it )


Italy_fit  = SEIR__fit_ext$incidence_it %>% as.data.frame.table() %>%
  rename(ni = iterations, time = Var2, variant = Var3, value = Freq) %>%
  dplyr::mutate(ni = as.numeric(ni),
                time = as.numeric(time)) %>%
  group_by(time,variant) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>% 
  ungroup( ) %>%  
  mutate(Date = rep(all_dates_it , each = 4)) %>% 
  left_join(dataDay_it , by = "Date") %>%
  mutate(Variant = ifelse(variant.x == "A", "M234I-A376T",
                          ifelse(variant.x == "B", "A220V",
                                 ifelse(variant.x=="C", "Other", "Alpha")))) %>% 
  ggplot(aes(x = Date , y = mean)) +
  geom_line(aes(color = Variant)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst ,color = variant.y)) +
  geom_errorbar(aes(ymin=Lower, ymax = Upper, color = variant.y)) +
  labs(y = paste0("")) +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "none") +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Italy")+
  ylim(c(0,150))





# Plot model Veneto -------------------------------------------------------------------
dataConf_veneto = dataDay_veneto  = list()

monthly_date_veneto = # account for seeding and aggregating from day to month
  all_dates_veneto [ which(as.numeric(format(all_dates_veneto, "%d")) == 1) ]

dataList_veneto = # data list 
  list(M_data_veneto[index_M_veneto], A_data_veneto[index_A_veneto], 
                O_data_veneto[index_O_veneto], Al_data_veneto[index_Al_veneto])  

indexList_veneto = # index data by month
  list(index_M_i_veneto , index_A_i_veneto, index_O_i_veneto, index_Al_i_veneto )  


# calculate binomial confidence intervals

for (i in 1:length(dataList_veneto )) {
  dataConf_veneto [[i]] = data.frame(
    Date = monthly_date_veneto [indexList_veneto [[i]]],
    binconf(dataList_veneto [[i]], n_seq_i_veneto[indexList_veneto [[i]]]) * 
      (avg_daily_rep_inc_i_veneto[indexList_veneto [[i]]] / (n_pop_veneto  -n_recov_veneto ) * 100000) ,
    variant = varaints[i]
  )
  
}

dataConfDf_veneto  = bind_rows(dataConf_veneto )
dataDay_veneto  = left_join(data.frame(Date = all_dates_veneto ),dataConfDf_veneto )


veneto_fit  = SEIR__fit_ext$incidence_ven %>% as.data.frame.table() %>%
  rename(ni = iterations, time = Var2, variant = Var3, value = Freq) %>%
  dplyr::mutate(ni = as.numeric(ni),
                time = as.numeric(time)) %>%
  group_by(time,variant) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>% 
  ungroup( ) %>%  
  mutate(Date = rep(all_dates_veneto, each = 4)) %>% 
  left_join(dataDay_veneto, by = "Date") %>%
  mutate(Variant = ifelse(variant.x == "A", "M234I-A376T",
                          ifelse(variant.x == "B", "A220V",
                                 ifelse(variant.x=="C", "Other", "Alpha")))) %>% 
  ggplot(aes(x = Date , y = mean)) +
  geom_line(aes(color = Variant)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst ,color = variant.y)) +
  geom_errorbar(aes(ymin=Lower, ymax = Upper, color = variant.y))+
  labs(y = paste0("Reported Incidence per 100,000 population")) +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c(0.12,0.86)) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Veneto") + 
  ylim(c(0,150))

#  Save Veneto and Italy fits --------------------------------------------------

rep_inc_plot = plot_grid(veneto_fit, Italy_fit)


ggsave(
  rep_inc_plot,
  file =  "reportedIncidencePlot_10L_location_omega.jpg",
  height = 20,
  width = 50,
  unit = "cm",
  dpi = 720
)
