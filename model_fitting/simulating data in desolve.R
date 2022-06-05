rm(list = ls())
setwd("Q:/COV_Italy_multistrain/model_fitting")


########## Load packages    ########## 
library(bayesplot)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(Hmisc)



# fake data
SEIQR = function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    
    if(time < time_intervention1 ) {
      foiM = betaM*(IM/S0)
      foiA = betaA*(IA/S0)
      foiO = betaO*(IO/S0)
      foiAl = betaAl*(IAl/S0)
    }  else if (time >= time_intervention1 & time < time_intervention2) {
      foiM = (1-omega1) * betaM*(IM/S0)
      foiA = (1-omega1) * betaA*(IA/S0)
      foiO = (1-omega1) * betaO*(IO/S0)
      foiAl = (1-omega1) * betaAl*(IAl/S0)
    } else if (time > time_intervention2){
      foiM = (1-omega2) * betaM*(IM/S0)
      foiA = (1-omega2) * betaA*(IA/S0)
      foiO = (1-omega2) * betaO*(IO/S0)
      foiAl = (1-omega2) * betaAl*(IAl/S0)
    }
    
    
    dS <- - S * (foiM + foiA + foiO + foiAl)
    
    dEM <- foiM * S - sigma*EM
    dIM <- (1-rho)*sigma*EM - gamma*IM 
    
    dEA <- foiA*S - sigma*EA
    dIA <- (1-rho)*sigma*EA - gamma*IA 
    
    dEO <- foiO*S - sigma*EO
    dIO <- (1-rho)*sigma*EO - gamma*IO 
    
    dEAl <- foiAl*S - sigma * EAl
    dIAl <- (1-rho) * sigma * EAl - gamma*IAl 
    
    
    dQ <- rho * sigma * (EO + EM + EA + EAl) - gamma * Q
    dR <- gamma*(IM+IA+IO+IAl+Q)
    
    return(list(c(dS, dEM, dIM,dEA, dIA,dEO, dIO,dEAl,dIAl, dQ, dR)))
  })
}



# load packages and source user defined functions -----------------------------#       

library(deSolve)
library(tidyverse)


# define global variables -----------------------------------------------------#  

# start and end data of model fitting 
start_date =  as.Date.character("01-07-2020", format = "%d-%m-%Y") 
end_date = as.Date.character("31-05-2021", format = "%d-%m-%Y")

# sequence of dates of model fitting 
all_dates = seq.Date(from = start_date, to = end_date ,  by = "days")

# date to seed m 
seed_m =  which(all_dates == as.Date.character("01-10-2020", format = "%d-%m-%Y")) 


# date to seed alpha 
seed_al =  which(all_dates == as.Date.character("01-11-2020", format = "%d-%m-%Y")) 


# timing of interventions 

time_intervention1 = which(all_dates == as.Date.character("15-11-2020", format = "%d-%m-%Y")) 

time_intervention2 = which(all_dates == as.Date.character("15-03-2021", format = "%d-%m-%Y")) 


# model times 
ts = 1:length(all_dates)

# model parameters 

R0M = 1.4 # reproduction number 
R0A = 1.5 # reproduction number 
R0O = 1.4 # reproduction number 
R0Al = 2.8

sigma = 1/5.1   # rate of progression
gamma = 1/2.1     # rate of recovery 
rho = 0.6         # reporting probability 
omega1 = 0.5
omega2 = 0.4
betaM =( R0M * gamma )/ (1-rho)    # betaM, the transmission rate = R0M / (1-rho) * gamma 
betaA =( R0A * gamma )/ (1-rho)    # betaA, the transmission rate = R0A / (1-rho) * gamma 
betaO =( R0O* gamma )/ (1-rho)    # betaA, the transmission rate = R0A / (1-rho) * gamma 
betaAl =( R0Al * gamma )/ (1-rho)    # betaAl, the transmission rate = R0A / (1-rho) * gamma 

# initial states 

n_pop = 4847026
n_recov = 93401
S0 = n_pop - n_recov
n_infA = n_infO = 10 # 1 initial infection 
n_infM = 500 # 1 initial infection 
n_infAl = 100
  
initial_state = c(S= S0 - n_infA - n_infO , EM = 0 , IM=0,EA = 0, IA = n_infA, EO = 0 , IO = n_infO,EAl = 0, IAl = 0, Q=0, R=n_recov)

params = c(betaM,betaA, betaO,betaAl, sigma, gamma,rho,S0, time_intervention1, time_intervention2, omega1, omega2)
# seed later 

seed_event =  data.frame(var = c("IM", "S", "IAl", "S"),
                           time = c(seed_m,seed_m,seed_al,seed_al), 
                           value = c(n_infM,-n_infM, n_infAl,-n_infAl),
                           method = c("add" ,"add", "add", "add"))




# Solve the model using deSolve -----------------------------------------------#   

model = ode(initial_state, ts, SEIQR, params, events = list(data = seed_event))

out.df  = round(as.data.frame(model))

# plot the model states to check they make sense ------------------------------#  

for(i in 2:ncol(out.df)) plot(out.df[,i])


# plot epidemic curve ---------------------------------------------------------#  

sim_data = data.frame(time = ts,
                      rep_incM = round(rho *out.df$EM* sigma),
                      rep_incA = round(rho * out.df$EA * sigma),
                      rep_incO = round(rho * out.df$EO * sigma),
                      rep_incAl = round(rho * out.df$EAl * sigma))

ggplot(sim_data, aes(x= time , y = rep_incA)) +
  geom_line() +
  geom_line(aes(y = rep_incM), color = "red") +
  geom_line(aes(y = rep_incO), color = "blue") +
  geom_line(aes(y = rep_incAl), color = "pink")
# + 
#   geom_point(data = incidencedf, aes(y = rep_incA)) +
#   geom_point(data = incidencedf,aes(y = rep_incM), color = "red") +
#   geom_point(data = incidencedf,aes(y = rep_incO), color = "blue") +
#   geom_point(data = incidencedf,aes(y = rep_incAl), color = "pink") 



# aggregate to months ---------------------------------------------------------#

# index by 1st of each month and the final day of the last month
index_1st_month = c(which(as.numeric(format(all_dates, "%d")) == 1), 
                    (tail(length(all_dates)) + 1)) 


n_variants = 4
n_months = round(length(all_dates)/30)
monthly_incidence = matrix(0,n_months, n_variants)

index = 1; 

for (m in 1:n_months){ 
  ind = index + 1;
  for (i in 1:4){
    monthly_incidence[m,i] =  mean( sim_data[index_1st_month[index]:(index_1st_month[ind]-1),i+1] );    }
  index = index + 1;
}

monthly_incidencedf = data.frame(
  time = 1:n_months,
  rep_incM = monthly_incidence[,1],
  rep_incA = monthly_incidence[,2],
  rep_incO = monthly_incidence[,3],
  rep_incAl = monthly_incidence[,4]
)
ggplot(monthly_incidencedf, aes(x= time , y = rep_incA)) +
  geom_line() +
  geom_line(aes(y = rep_incM), color = "red") +
  geom_line(aes(y = rep_incO), color = "blue") +
  geom_line(aes(y = rep_incAl), color = "pink") 

################################################################################
# Fit model to simulated data assuming vaccination = 0 and pPCR = 1
################################################################################


seedA = 1
seedO = 1
seedM = 3
seedAl =4

phi_PCR = 0.920
phi_Ag =  0.643 

model_data = list(
  n_data_A  =  (n_months-seedA),
  n_data_M =  (n_months - seedM),
  n_data_O =   (n_months-seedO),
  n_data_Al =  (n_months - seedAl),
  
  n_pop = n_pop,
  n_days = length(all_dates),
  n_recov = n_recov,
  n_months = n_months,
  n_variants = 4,
  n_ts =length(all_dates),
  

  y_M = round(monthly_incidence[(seedM+1):n_months,1]),
  y_A = round(monthly_incidence[(seedA+1):n_months,2]),
  y_O = round(monthly_incidence[(seedO+1):n_months,3]),
  y_Al = round(monthly_incidence[(seedAl+1):n_months,4]),

  
  gamma = gamma ,
  sigma = sigma ,
  phi_Ag =  phi_Ag,
  phi_PCR = phi_PCR , 
  
  
  Ag_daily  = rep(0, length(all_dates)),
  PCR_daily = rep(1,length(all_dates)),
  vac = rep(0, length(all_dates)),
  
  time_switch1 = time_intervention1,
  time_switch2 =time_intervention2,
  
  time_seed_M = seed_m,
  time_seed_alpha = seed_al,
  
  month_index = index_1st_month
)



##########  set up model   ########## 


SEIR_model = stan_model("model/new_neg_bin_sens_estimate_new_likelihood.stan")


ini_1_SEIR = function(){
  list(  beta = replicate(4,runif(1,1,4)),
         I0 = replicate(4, runif(1, 1,500)),
         omega = replicate(2,runif(1,0.2,0.6)),
         rho = runif(1,0.2,0.6),
         k = runif(1,0.01,1)
  )}

##########  run model    ########## 
time.start <- Sys.time()

SEIR__fit = sampling(
  SEIR_model,
  data = model_data,
  init = ini_1_SEIR,
  chains = 3,
  warmup = 1000,
  iter = 2000,
  thin = 1 , 
   control = list(
    # adapt_delta = 0.8, 
     max_treedepth = 12
   )
  ) 

time.end <- Sys.time()
time.end - time.start 


pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
         , "I0[1]", "I0[2]", "I0[3]", "I0[4]")

SEIR__fit_summary = summary(SEIR__fit, pars = pars)$summary



SEIR_fit_posterior = as.array(SEIR__fit)

# Markov chain traceplots

mcmc_trace(SEIR_fit_posterior, pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho"  ))
mcmc_pairs(SEIR_fit_posterior, pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho"   ))


SEIR__fit_ext = rstan::extract(SEIR__fit)

monthly_incidencedf_long =monthly_incidencedf %>% 
  pivot_longer(!time, names_to = "variant", values_to="data") %>% 
  mutate(variant = rep(c("A", "B", "C", "D"), n_months))

plot_m1_solv_results_real  = SEIR__fit_ext$monthly_incidence %>% as.data.frame.table() %>%
  rename(ni = iterations, time = Var2, variant = Var3, value = Freq) %>%
  dplyr::mutate(ni = as.numeric(ni),
                time = as.numeric(time)) %>%
  group_by(time,variant) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>% 
  left_join(monthly_incidencedf_long) %>%
  ggplot(aes(x = time , y = mean)) +
  geom_line(aes(color = variant)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = variant), alpha = 0.4) +
  geom_point(aes(y = data,color = variant))



# -----------------------------------------------------------------------------
A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq 
M_data = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq 
O_data = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq 
Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq 
n_seq = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq

index_M = 8:14
index_A = c(5,7:14)
index_O = c(5,7:9)
index_Al = 9:14  

average_daily_reported_incidence = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily

### having extracted GISAID data from first time of fitting (first occurrence of A variant), this allows to index variant specific data, removing NA values 

index_M_i = index_M - (length(M_data)- n_months)
index_A_i = index_A - (length(A_data)- n_months)
index_O_i = index_O - (length(O_data)- n_months)
index_Al_i  = index_Al - (length(Al_data)- n_months)



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
  phi_Ag =  phi_Ag,
  phi_PCR = phi_PCR , 
  
  
  
  
  Ag_daily  = round(daily_Ag_i2),
  PCR_daily = round(daily_PCR_i2),
  vac = average_daily_vaccination_i2,
  
  time_switch1 = time_intervention1,
  time_switch2 =time_intervention2,
  
  time_seed_M = seed_m,
  time_seed_alpha = seed_al,
  
  month_index = index_1st_month
)


#-------------------------------------------------------------------------------


SEIR_model = stan_model("model/new_neg_bin_sens_estimate_new_likelihood.stan")


time.start <- Sys.time()

SEIR__fit = sampling(
  SEIR_model,
  data = model_data_real,
  init = ini_1_SEIR,
  chains = 3,
  warmup = 1000,
  iter = 2000,
  thin = 1 , 
  control = list(
    # adapt_delta = 0.8, 
    max_treedepth = 12
  )
) 

time.end <- Sys.time()
time.end - time.start 


#-------------------------------------------------------------------------------
pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
         , "I0[1]", "I0[2]", "I0[3]", "I0[4]")

SEIR__fit_summary = summary(SEIR__fit, pars = pars)$summary



SEIR_fit_posterior = as.array(SEIR__fit)

# Markov chain traceplots

mcmc_trace(SEIR_fit_posterior, pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho"   ))
mcmc_pairs(SEIR_fit_posterior, pars = c( "beta[1]" , "beta[2]" , "beta[3]"  ,"beta[4]" , "rho"   ))

#-------------------------------------------------------------------------------
dataConf = dataDay = list()
# account for seeding and aggregating from day to month

monthly_date = daily_date[ which(as.numeric(format(daily_date, "%d")) == 1) ]

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
dataDay = left_join(data.frame(Date = daily_date),dataConfDf)




SEIR__fit_ext = rstan::extract(SEIR__fit)

plot_m1_solv_results_real  = SEIR__fit_ext$incidence %>% as.data.frame.table() %>%
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
  mutate(Date = rep(daily_date, each = 4)) %>% 
  left_join(dataDay, by = "Date") %>%
  ggplot(aes(x = Date , y = mean)) +
  geom_line(aes(color = variant.x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = variant.x), alpha = 0.4) +
  geom_point(aes(y = PointEst ,color = variant.y)) +
  geom_errorbar(aes(ymin=Lower, ymax = Upper, color = variant.y))

