
################  Run all models and counterfactuals ################  



################   Set up  ################  


# rm(list = ls())

library(tidyverse)
library(ggplot2)
library(rstan)
library(bayesplot)
library(ggpubr)
library(Hmisc)
library("ggsci")
library(cowplot)
library("wesanderson")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd(
  "C:/Users/bnc19/OneDrive - Imperial College London/SEIR italy cov/COV_Italy_multistrain_2/counterfactuals_V2/run_fixed_rstan"
)

source("R/replicate_rstan_fit.R")


################   Define values  ################  

number_of_samples =100

n_difeq = 11
gamma = 1/2.1
sigma = 1/5.1
phi_PCR = 0.92
phi_Ag = 0.643

n_pop = 4847026
n_recov = 93401
n_pop_italy = (59257566 - 4847026 )
n_recov_italy = (1482377 - 93401)


R0_scales = c(".8", "1", "1.2","1.4", "1.6")


################   Import data  ################  


daily_Ag_i_italy = read.csv("Italy/daily_Ag_i_italy.csv")[, 2]
daily_PCR_i_italy = read.csv("Italy/daily_PCR_i_italy.csv")[, 2]
average_daily_vaccination_i_italy  = read.csv("Italy/average_daily_vaccination_i_italy.csv")[, 2]
x_i_data_italy  = read.csv("Italy/x_i_data_italy.csv")[, 2]

daily_Ag_i = read.csv("Veneto/daily_Ag_i.csv")[, 2]
daily_PCR_i = read.csv("Veneto/daily_PCR_i.csv")[, 2]
average_daily_vaccination_i = read.csv("Veneto/average_daily_vaccination_i.csv")[, 2]
x_i_data = read.csv("Veneto/x_i_data.csv")[, 2]
x_i_data_italy_test = read.csv("Veneto/x_i_data_Veneto_italy_test.csv")[, 2]
daily_Ag_i_Itest = read.csv("Veneto/daily_Ag_i_Itest.csv")[, 2]
daily_PCR_i_Itest = read.csv("Veneto/daily_PCR_i.csv")[, 2]


n_months = x_i_data[1]
n_days = x_i_data[(3 * n_months + 5)]

n_months_italy  = x_i_data_italy[1]
n_days_italy  =  x_i_data_italy[(3 * n_months_italy + 5)]

# rstan model
SEIR_model = stan_model("Models/stan_model.stan")
SEIR_model_pPCR_0 = stan_model("Models/stan_model_pPCR_0.stan")
SEIR_model_PCR_after_ag = stan_model("Models/stan_model_PCR_after_ag.stan")


# posterior samples 

posterior_samples = read.csv( "Veneto/100_posterior_samples_veneto.csv")[,-1]
posterior_samples_italy = read.csv("Italy/100_posterior_samples_italy.csv")[,-1]



################   Define column names for output  ################  


col_names = c(
  "M_fit_M",
  "A_fit_M" ,
  "O_fit_M",
  "Al_fit_M",
  "T_M_fit_M",
  "T_A_fit_M" ,
  "T_O_fit_M",
  "T_Al_fit_M",
  "total_M",
  "total_rep_M",

  
 "M_R0_M",
 "A_R0_M" ,
 "O_R0_M",
 "Al_R0_M",
 "pPCR_M",
  
  "M_fit_L",
  "A_fit_L" ,
  "O_fit_L",
  "Al_fit_L",
  "T_M_fit_L",
  "T_A_fit_L" ,
  "T_O_fit_L",
  "T_Al_fit_L",
  "total_L",
 "total_rep_L",


 
 "M_R0_L",
 "A_R0_L" ,
 "O_R0_L",
 "Al_R0_L",
 "pPCR_L",
  
  "M_fit_U",
  "A_fit_U" ,
  "O_fit_U",
  "Al_fit_U",
  "T_M_fit_U",
  "T_A_fit_U" ,
  "T_O_fit_U",
  "T_Al_fit_U",
  "total_U",
 "total_rep_U",

 "M_R0_U",
 "A_R0_U" ,
 "O_R0_U",
 "Al_R0_U",
 "pPCR_U"
  
)




################## Model 1: baseline scenario ################## 

#### apply function across 1000 bootstrap samples of posterior estimates


##  TAKES 1-2 MINUTES TO RUN ##

modelfit1  = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior_samples[i, ],
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
    x_i_data,
    SEIR_model
  )
}, simplify = "array")

## model fit 
  
mean1  = data.frame(apply(modelfit1, c(1, 2),  mean))
lower1 = data.frame(apply(modelfit1, c(1, 2), quantile, probs = 0.025))
upper1 = data.frame(apply(modelfit1, c(1, 2), quantile, probs = 0.975))

modelfit1_df = data.frame(mean1, lower1, upper1)
colnames(modelfit1_df) = col_names


write.csv(modelfit1_df, "Results/modelfit1_df_veneto.csv")

## cumulative incidence 
cuminc1  = data.frame(apply (modelfit1, 2, colSums))



mean(9538  /  cuminc1$total_incidence   * 100) # calculate IFR from baseline scenario as internal check  


cuminc1 = cuminc1 %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)

mean_cum1 = data.frame(sapply(cuminc1, mean))
lower_cum1 = data.frame(sapply(cuminc1, quantile, probs = 0.025))
upper_cum1 = data.frame(sapply(cuminc1, quantile, probs = 0.975))


cuminc1_df = data.frame(t(cbind(mean_cum1, lower_cum1, upper_cum1)))
rownames(cuminc1_df) = c("mean", "lower", "upper")

write.csv(cuminc1_df, "Results/cuminc1_veneto.csv")



##################  Model 2: baseline testing, scale B_M ##################



posterior_samples.8 = posterior_samples
posterior_samples.1.2 = posterior_samples
posterior_samples.1.4 = posterior_samples
posterior_samples.1.6 = posterior_samples

posterior_samples.8$beta_chain_M = (posterior_samples$beta_chain_M * .8)
posterior_samples.1.2$beta_chain_M = (posterior_samples$beta_chain_M * 1.2)
posterior_samples.1.4$beta_chain_M = (posterior_samples$beta_chain_M * 1.4)
posterior_samples.1.6$beta_chain_M = (posterior_samples$beta_chain_M * 1.6)

beta_scale_posterior_samples_list = list(
  posterior_samples.8,
  posterior_samples,
  posterior_samples.1.2,
  posterior_samples.1.4,
  posterior_samples.1.6
)

#### apply function across 1000 bootstrap samples of posterior estimates


## TAKES 5 MINUTES TO RUN ##

modelfit2 = lapply(1:length(beta_scale_posterior_samples_list), function(x) {
  sapply(1:number_of_samples, function(i) {
    replicate_rstan_fixed(
      beta_scale_posterior_samples_list[[x]][i, ],
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
      x_i_data,
      SEIR_model
    )
  }, simplify = "array")
})




cuminc2  = lapply(1:length(modelfit2), function(i) {
  data.frame(apply (modelfit2[[i]], 2 , colSums))
})

mean_cum2 = lapply(1:length(modelfit2), function(i) {
  data.frame(sapply(cuminc2[[i]], mean))
})
lower_cum2 = lapply(1:length(modelfit2), function(i) {
  data.frame(sapply(cuminc2[[i]], quantile, probs = 0.025))
})
upper_cum2 = lapply(1:length(modelfit2), function(i) {
  data.frame(sapply(cuminc2[[i]], quantile, probs = 0.975))
})


for(i in 1:5){
  rownames(mean_cum2[[i]]) = col_names[1:15]
  rownames(lower_cum2[[i]]) = col_names[16:30]
  rownames(upper_cum2[[i]]) = col_names[31:45]
  
  colnames(mean_cum2[[i]]) = R0_scales[i]
  colnames(lower_cum2[[i]]) = R0_scales[i]
  colnames(upper_cum2[[i]]) = R0_scales[i]
  
 }
cum2list = list(mean_cum2, lower_cum2, upper_cum2)
cum2list2 = list()
for (i in 1:3){
  cum2list2[[i]]  = bind_cols(cum2list[[i]])
}

cuminc2_df = bind_rows(cum2list2)
write.csv(cuminc2_df, "Results/cuminc2_veneto.csv")


##################  Model 3: No PCR, baseline M ################## 

################## Only antigen testing, estimate rho, 0.643 test sens ################## ###################

modelfit3  = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples[i,],
    n_difeq = n_difeq,
    n_pop = n_pop ,
    n_recov = n_recov ,
    n_months =  n_months,
    gamma = gamma,
    sigma = sigma,
    phi_PCR,
    phi_Ag,
    n_days = n_days,
    daily_Ag_i =  daily_Ag_i,
    daily_PCR_i = daily_PCR_i,
    average_daily_vaccination_i = average_daily_vaccination_i,
    x_i_data =  x_i_data,
    SEIR_model = SEIR_model_pPCR_0
  )
}, simplify = "array")

cuminc3  = data.frame(apply (modelfit3, 2 , colSums))
mean_cum3 = data.frame(sapply(cuminc3, mean))
lower_cum3 = data.frame(sapply(cuminc3, quantile, probs = 0.025))
upper_cum3 = data.frame(sapply(cuminc3, quantile, probs = 0.975))
################## Only antigen testing, rho = 1, 0.643 test sens      ##################

posterior_samples_2_rho = posterior_samples

posterior_samples_2_rho$rho_chain = 1

modelfit3_rho = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples_2_rho[i,],
    n_difeq = n_difeq,
    n_pop = n_pop ,
    n_recov = n_recov ,
    n_months =  n_months,
    gamma = gamma,
    sigma = sigma,
    phi_PCR,
    phi_Ag,
    n_days = n_days,
    daily_Ag_i =  daily_Ag_i,
    daily_PCR_i = daily_PCR_i,
    average_daily_vaccination_i = average_daily_vaccination_i,
    x_i_data =  x_i_data,
    SEIR_model = SEIR_model_pPCR_0
  )
}, simplify = "array")



# cumulative data 

cuminc3_rho  = data.frame(apply (modelfit3_rho, 2 , colSums))

cuminc3_rho = cuminc3_rho %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)

mean_cum3_rho = data.frame(sapply(cuminc3_rho, mean))
lower_cum3_rho = data.frame(sapply(cuminc3_rho, quantile, probs = 0.025))
upper_cum3_rho = data.frame(sapply(cuminc3_rho, quantile, probs = 0.975))

cuminc3_df = data.frame(t(cbind(mean_cum3_rho, lower_cum3_rho, upper_cum3_rho)))
rownames(cuminc3_df) = c("mean", "lower", "upper")


write.csv(cuminc3_df, "Results/cuminc3_df_veneto.csv")

################## Only antigen testing, rho = 1, 0.875 test sens  ##################     ##################

modelfit3_rho_Ag_sens = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples_2_rho[i,],
    n_difeq = n_difeq,
    n_pop = n_pop ,
    n_recov = n_recov ,
    n_months =  n_months,
    gamma = gamma,
    sigma = sigma,
    phi_PCR,
    0.875,
    n_days = n_days,
    daily_Ag_i =  daily_Ag_i,
    daily_PCR_i = daily_PCR_i,
    average_daily_vaccination_i = average_daily_vaccination_i,
    x_i_data =  x_i_data,
    SEIR_model = SEIR_model_pPCR_0
  )
}, simplify = "array")

cuminc3_rho_Ag_sens  = data.frame(apply (modelfit3_rho_Ag_sens, 2 , colSums))




cuminc3_rho_Ag_sens = cuminc3_rho_Ag_sens %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)

mean_cum3_rho_Ag_sens = data.frame(sapply(cuminc3_rho_Ag_sens, mean))
lower_cum3_rho_Ag_sens = data.frame(sapply(cuminc3_rho_Ag_sens, quantile, probs = 0.025))
upper_cum3_rho_Ag_sens = data.frame(sapply(cuminc3_rho_Ag_sens, quantile, probs = 0.975))



cuminc3_df_sens = data.frame(t(cbind(mean_cum3_rho_Ag_sens, lower_cum3_rho_Ag_sens, upper_cum3_rho_Ag_sens)))
rownames(cuminc3_df_sens) = c("mean", "lower", "upper")

write.csv(cuminc3_df_sens, "Results/cuminc3_df_sens_veneto.csv")




################## Model 4: No PCR, scale B_M, mass testing  ################  



#### apply function across 1000 bootstrap samples of posterior estimates

beta_scale_posterior_samples_list_rho = lapply(beta_scale_posterior_samples_list, change_rho)

# ##  TAKES 1-2 MINUTES TO RUN ##

modelfit4 = lapply(1:length(beta_scale_posterior_samples_list), function(x) {
  sapply(1:number_of_samples, function(i) {
    replicate_rstan_fixed(
      beta_scale_posterior_samples_list_rho[[x]][i, ],
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
      x_i_data,
      SEIR_model = SEIR_model_pPCR_0
    )
  }, simplify = "array")
})




cuminc4  = lapply(1:length(modelfit4), function(i) {
  data.frame(apply (modelfit4[[i]], c(2) , colSums))
})


mean_cum4 = lapply(1:length(modelfit4), function(i) {
  data.frame(sapply(cuminc4[[i]], mean))
})
lower_cum4 = lapply(1:length(modelfit4), function(i) {
  data.frame(sapply(cuminc4[[i]], quantile, probs = 0.025))
})
upper_cum4 = lapply(1:length(modelfit4), function(i) {
  data.frame(sapply(cuminc4[[i]], quantile, probs = 0.975))
})




for(i in 1:5){
  rownames(mean_cum4[[i]]) = col_names[1:15]
  rownames(lower_cum4[[i]]) = col_names[16:30]
  rownames(upper_cum4[[i]]) = col_names[31:45]
  
  
  colnames(mean_cum4[[i]]) = R0_scales[i]
  colnames(lower_cum4[[i]]) = R0_scales[i]
  colnames(upper_cum4[[i]]) = R0_scales[i]
  
}
cum4list = list(mean_cum4, lower_cum4, upper_cum4)
cum4list2 = list()
for (i in 1:3){
  cum4list2[[i]]  = bind_cols(cum4list[[i]])
}

cuminc4_df = bind_rows(cum4list2)
write.csv(cuminc4_df, "Results/cuminc4_df_veneto.csv")

 
################## Model 5: PCR after neg Ag, baseline M  ################  



#### apply function across 1000 bootstrap samples of posterior estimates

##  TAKES 1-2 MINUTES TO RUN ##

modelfit5  = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples[i,],
    n_difeq = n_difeq,
    n_pop = n_pop ,
    n_recov = n_recov ,
    n_months =  n_months,
    gamma = gamma,
    sigma = sigma,
    phi_PCR,
    phi_Ag,
    n_days = n_days,
    daily_Ag_i =  daily_Ag_i,
    daily_PCR_i = daily_PCR_i,
    average_daily_vaccination_i = average_daily_vaccination_i,
    x_i_data =  x_i_data,
    SEIR_model = SEIR_model_PCR_after_ag
  )
}, simplify = "array")

cuminc5  = data.frame(apply (modelfit5, 2 , colSums))


cuminc5 = cuminc5 %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)


mean_cum5 = data.frame(sapply(cuminc5, mean))
lower_cum5 = data.frame(sapply(cuminc5, quantile, probs = 0.025))
upper_cum5 = data.frame(sapply(cuminc5, quantile, probs = 0.975))


cuminc5_df = data.frame(t(cbind(mean_cum5, lower_cum5, upper_cum5)))
rownames(cuminc5_df) = c("mean", "lower", "upper")


write.csv(cuminc5_df, "Results/cuminc5_df_veneto.csv")

################## Model 6: Italy testing, scale B_M ################## 



#### apply function across 1000 bootstrap samples of posterior estimates

##  TAKES 1-2 MINUTES TO RUN ##

modelfit6 = lapply(1:length(beta_scale_posterior_samples_list), function(x) {
  sapply(1:number_of_samples, function(i) {
    replicate_rstan_fixed(
      beta_scale_posterior_samples_list[[x]][i, ],
      n_difeq,
      n_pop ,
      n_recov ,
      n_months,
      gamma,
      sigma,
      phi_PCR,
      phi_Ag,
      n_days,
      daily_Ag_i_Itest,
      daily_PCR_i_Itest,
      average_daily_vaccination_i,
      x_i_data_italy_test,
      SEIR_model = SEIR_model
    )
  }, simplify = "array")
})




cuminc6  = lapply(1:length(modelfit6), function(i) {
  data.frame(apply (modelfit6[[i]], c(2) , colSums))
})


mean_cum6 = lapply(1:length(modelfit6), function(i) {
  data.frame(sapply(cuminc6[[i]], mean))
})
lower_cum6 = lapply(1:length(modelfit6), function(i) {
  data.frame(sapply(cuminc6[[i]], quantile, probs = 0.025))
})
upper_cum6 = lapply(1:length(modelfit6), function(i) {
  data.frame(sapply(cuminc6[[i]], quantile, probs = 0.975))
})



for(i in 1:5){
  rownames(mean_cum6[[i]]) = col_names[1:15]
  rownames(lower_cum6[[i]]) = col_names[16:30]
  rownames(upper_cum6[[i]]) = col_names[31:45]
  
  
  colnames(mean_cum6[[i]]) = R0_scales[i]
  colnames(lower_cum6[[i]]) = R0_scales[i]
  colnames(upper_cum6[[i]]) = R0_scales[i]
  
}
cum6list = list(mean_cum6, lower_cum6, upper_cum6)
cum6list2 = list()
for (i in 1:3){
  cum6list2[[i]]  = bind_cols(cum6list[[i]])
}

cuminc6_df = bind_rows(cum6list2)
write.csv(cuminc6_df, "Results/cuminc6_df_veneto.csv")



################ Model 7: Italy testing in Veneto ################   



#### apply function across 1000 bootstrap samples of posterior estimates

##  TAKES 1-2 MINUTES TO RUN ##


modelfit7  = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples[i, ],
    n_difeq = n_difeq,
    n_pop = n_pop ,
    n_recov =  n_recov ,
    n_months =  n_months,
    gamma =   gamma,
    sigma =   sigma,
    phi_PCR,
    phi_Ag,
    n_days =  n_days,
    daily_Ag_i =   daily_Ag_i_Itest,
    daily_PCR_i =   daily_PCR_i_Itest,
    average_daily_vaccination_i =   average_daily_vaccination_i,
    x_i_data =     x_i_data_italy_test,
    SEIR_model =    SEIR_model
  )
}, simplify = "array")


cuminc7  = data.frame(apply (modelfit7, 2 , colSums))



cuminc7 = cuminc7 %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)


mean_cum7 = data.frame(sapply(cuminc7, mean))
lower_cum7 = data.frame(sapply(cuminc7, quantile, probs = 0.025))
upper_cum7 = data.frame(sapply(cuminc7, quantile, probs = 0.975))


cuminc7_df = data.frame(t(cbind(mean_cum7, lower_cum7, upper_cum7)))
rownames(cuminc7_df) = c("mean", "lower", "upper")

write.csv(cuminc7_df, "Results/cuminc7_df_veneto.csv")


################  Model 1 Italy : baseline scenario ################ 


#### apply function across 1000 bootstrap samples of posterior estimates


##  TAKES 1-2 MINUTES TO RUN ##

modelfit1_italy  = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples_italy[i, ],
    n_difeq = n_difeq,
    n_pop = n_pop_italy ,
    n_recov =  n_recov_italy ,
    n_months =  n_months_italy,
    gamma =   gamma,
    sigma =   sigma,
    phi_PCR,
    phi_Ag,
    n_days =  n_days_italy,
    daily_Ag_i =   daily_Ag_i_italy,
    daily_PCR_i =   daily_PCR_i_italy,
    average_daily_vaccination_i =   average_daily_vaccination_i_italy,
    x_i_data =     x_i_data_italy,
    SEIR_model =    SEIR_model
  )
}, simplify = "array")


## model fit 

mean1_italy  = data.frame(apply (modelfit1_italy, c(1, 2),  mean))
lower1_italy = data.frame(apply(modelfit1_italy, c(1, 2), quantile, probs = 0.025))
upper1_italy = data.frame(apply(modelfit1_italy, c(1, 2), quantile, probs = 0.975))

modelfit1_df_italy = data.frame(mean1_italy, lower1_italy, upper1_italy)
colnames(modelfit1_df_italy) = col_names

write.csv(modelfit1_df_italy, "Results/modelfit1_df_italy.csv")

## cumulative incidence 

cuminc1_italy  = data.frame(apply (modelfit1_italy, 2, colSums))


mean(87811  /  cuminc1_italy$total_incidence   * 100) # calculate IFR from baseline scenario as internal chec




cuminc1_italy = cuminc1_italy %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)

mean_cum1_italy = data.frame(sapply(cuminc1_italy, mean))
lower_cum1_italy = data.frame(sapply(cuminc1_italy, quantile, probs = 0.025))
upper_cum1_italy = data.frame(sapply(cuminc1_italy, quantile, probs = 0.975))


cuminc1_df_italy = data.frame(t(cbind(mean_cum1_italy, lower_cum1_italy, upper_cum1_italy)))
rownames(cuminc1_df_italy) = c("mean", "lower", "upper")


write.csv(cuminc1_df_italy, "Results/cuminc1_df_italy.csv")
