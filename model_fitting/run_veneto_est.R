rm(list = ls())
setwd("Q:/COV_Italy_multistrain/model_fitting")
modelPath = "model/new_neg_bin_sens_estimate_new_likelihood_presymp.stan"

# Load packages ----------------------------------------------------------------
library(bayesplot)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(Hmisc)


# Import data - ----------------------------------------------------------------
A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq 
M_data = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq 
O_data = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq 
Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq 
n_seq = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq
average_daily_reported_incidence = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily
daily_PCR = round(read.csv("data/Veneto_daily_test_data.csv")$pcr_daily)
daily_Ag= round(read.csv("data/Veneto_daily_test_data.csv")$antigen_daily)
average_daily_vaccination =  read.csv("data/data_vac_veneto_day.csv")$prop_vac

# Define variables -------------------------------------------------------------

Location = "Veneto"

# initial conditions
n_pop = 4847026
n_recov = 93401

# index data 
index_M = 8:14
index_A = c(5,7:14)
index_O = c(5,7:9)
index_Al = 9:14


# parameters 


#5.1 days = incubation period (no symptoms) 
epsilon = 1/ 2.8 ## mean of values given in Vo table S5 assuming R0 = 2.7

# 1/ sigma = 5.1 - 2.8 

sigma = 1 / (5.1 - 2.8 )
gamma = 1 / 2.1 
mu = 0.59 # probability symptomatic 

phi_PCR = 0.920
phi_Ag =  0.643 # change antigen sensitivity 


# dates 
start_date = "01-07-2020"
end_date = "31-05-2021"
time_intervention = c("15-11-2020" , "15-03-2021")
time_seed_alpha = "01-11-2020"
time_seed_M = "01-10-2020"
time_vac = "27-12-2020" # date first vaccine (first datapoint)


# stan stuff 
n_chains =4
n_warmups =2000
n_iter = 4000
n_thin = 1
pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", 
         "rho"  , "omega[1]"  , "omega[2]"
         , "I0[1]", "I0[2]", "I0[3]", "I0[4]")

ini_1_SEIR = function(){
  list(  beta = replicate(4,runif(1,0,3)),
         I0 = replicate(4, runif(1, 1,200)),
         omega = replicate(2,runif(1,0.2,0.8)),
         rho = runif(1,0.2,0.8),
         k = runif(1,0.01,2)
  )}



# format dates and find index --------------------------------------------------


start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
end_date_d = as.Date.character(end_date, format = "%d-%m-%Y")

all_dates = seq.Date(from = start_date_d, to = end_date_d ,  by = "days")

index_switch1 =  which(all_dates == as.Date.character(time_intervention[1], format = "%d-%m-%Y"))
index_switch2 =  which(all_dates == as.Date.character(time_intervention[2], format = "%d-%m-%Y"))
index_seed_alpha =  which(all_dates == as.Date.character(time_seed_alpha, format = "%d-%m-%Y"))
index_seed_M =  which(all_dates == as.Date.character(time_seed_M, format = "%d-%m-%Y"))
time_vac =  which(all_dates == as.Date.character(time_vac, format = "%d-%m-%Y"))

n_days = length(all_dates)
n_months = round(n_days / 30)

index_1st_month = c(which(as.numeric(format(all_dates, "%d")) == 1), (tail(n_days) + 1)) 
# index by 1st of each month and the final day of the last month


# extract correct dates of data ------------------------------------------------


daily_PCR_i = daily_PCR[(length(daily_PCR) - n_days + 1):length(daily_PCR)] 
daily_Ag_i = daily_Ag[(length(daily_Ag) - n_days + 1):length(daily_Ag)]  

average_daily_vaccination_i = rep(0, n_days)
average_daily_vaccination_i[time_vac:n_days] = average_daily_vaccination[1:length(time_vac:n_days)]

n_seq_i = n_seq[(length(n_seq) - n_months + 1):length(n_seq)]
average_daily_reported_incidence_i = average_daily_reported_incidence[(length(average_daily_reported_incidence)-n_months+1): length(average_daily_reported_incidence)]




# having extracted GISAID data from first date of fitting ----------------------
# index variant specific data, removing NA values ------------------------------

index_M_i = index_M - (length(M_data)- n_months)
index_A_i = index_A - (length(A_data)- n_months)
index_O_i = index_O - (length(O_data)- n_months)
index_Al_i  = index_Al - (length(Al_data)- n_months)


# estimate variant specific reported incidence ---------------------------------



n_reported_M = round(average_daily_reported_incidence * M_data / n_seq)[index_M]
n_reported_A = round(average_daily_reported_incidence * A_data / n_seq)[index_A]
n_reported_O = round(average_daily_reported_incidence * O_data / n_seq)[index_O]
n_reported_Al = round(average_daily_reported_incidence * Al_data / n_seq)[index_Al]


model_data_real = list(
  n_data_A  =  length(index_A),
  n_data_M =  length(index_M),
  n_data_O =  length(index_O),
  n_data_Al =  length(index_Al),
  
  index_M = index_M_i,
  index_A = index_A_i,
  index_O = index_O_i,
  index_Al = index_Al_i,
  
  n_pop = n_pop,
  n_days = length(all_dates),
  n_recov = n_recov,
  n_months = n_months,
  n_variants = 4,
  n_ts =length(all_dates),
  
  
  y_M = round(n_reported_M / (n_pop-n_recov) * 100000),
  y_A = round(n_reported_A / (n_pop-n_recov) * 100000),
  y_O = round(n_reported_O /  (n_pop-n_recov) * 100000),
  y_Al = round(n_reported_Al/  (n_pop-n_recov) * 100000),
  
  
  gamma = gamma ,
  sigma = sigma ,
  mu = mu, 
  epsilon = epsilon, 
  phi_Ag =  phi_Ag,
  phi_PCR = phi_PCR , 
  
  Ag_daily  = daily_Ag_i,
  PCR_daily = daily_PCR_i,
  vac = average_daily_vaccination_i,
  
  time_switch1 = index_switch1,
  time_switch2 =index_switch2,
  
  time_seed_M = index_seed_M,
  time_seed_alpha = index_seed_alpha,
  
  month_index = index_1st_month
)


#-------------------------------------------------------------------------------


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


#-------------------------------------------------------------------------------

SEIR__fit_summary = summary(SEIR__fit, pars = pars)$summary



SEIR_fit_posterior = as.array(SEIR__fit)

# Markov chain traceplots

mcmc_trace(SEIR_fit_posterior, pars =pars)
mcmc_pairs(SEIR_fit_posterior, pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho"   ))
mcmc_pairs(SEIR_fit_posterior, pars = c( "beta[1]" , "I0[1]" , "rho", "inv_epsilon", "k", "omega[1]"   ))
#-------------------------------------------------------------------------------
dataConf = dataDay = list()
# account for seeding and aggregating from day to month

monthly_date = all_dates[ which(as.numeric(format(all_dates, "%d")) == 1) ]

dataList = list(M_data[index_M], A_data[index_A], O_data[index_O], Al_data[index_Al])  # data
indexList = list(index_M_i, index_A_i, index_O_i, index_Al_i)  # index data by month
varaints = c("A", "B", "C", "D")
# calculate binomial confidence intervals

for (i in 1:length(dataList)) {
  dataConf[[i]] = data.frame(
    Date = monthly_date[indexList[[i]]],
    binconf(dataList[[i]], n_seq_i[indexList[[i]]]) * (average_daily_reported_incidence_i[indexList[[i]]] / (n_pop -n_recov) * 100000) ,
    variant = varaints[i]
  )
  
}

dataConfDf = bind_rows(dataConf)
dataDay = left_join(data.frame(Date = all_dates),dataConfDf)




SEIR__fit_ext = rstan::extract(SEIR__fit)

veneto_fit  = SEIR__fit_ext$incidence %>% as.data.frame.table() %>%
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
  mutate(Date = rep(all_dates, each = 4)) %>% 
  left_join(dataDay, by = "Date") %>%
  ggplot(aes(x = Date , y = mean)) +
  geom_line(aes(color = variant.x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = variant.x), alpha = 0.4) +
  geom_point(aes(y = PointEst ,color = variant.y)) +
  geom_errorbar(aes(ymin=Lower, ymax = Upper, color = variant.y))

