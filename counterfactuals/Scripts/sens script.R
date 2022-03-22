
################    Run all models and counterfactuals ################   




################    Set up ################   


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


################    Define values ################   


# data to run rstan model
set.seed(23222)
number_of_samples = 100

n_difeq = 11
gamma = 1/2.1
sigma = 1/5.1
phi_PCR = 0.92
phi_Ag_sens = 0.875

n_pop = 4847026
n_recov = 93401
n_pop_italy  = 59257566 - 4847026
n_recov_italy  = 1482377 - 93401




################    Import data ################   


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

variant_data = read.csv("Veneto/Veneto_variant_data.csv")[, -1]
variant_data$Date = as.Date.character(variant_data$Date, format = "%Y-%m-%d")

variant_data_italy = read.csv("Italy/Italy_variant_data.csv")[,-1]
variant_data_italy$Date = as.Date.character(variant_data_italy$Date, format = "%Y-%m-%d")

# posterior chains Veneto
posterior_chains_sens = read.csv("Veneto/posteriorChains_sens.csv")[-1]


# posterior chains Italy
posterior_chains_italy_sens = read.csv("Italy/posteriorChains_sens.csv")[-1]



################    Sample from posterior chains ################   

set.seed(23122023)

# sample for Italy and Veneto
posterior_index_sens = sample(size = number_of_samples,
                         1:dim(posterior_chains_sens)[1],
                         replace = T)
set.seed(23122022)

posterior_index_italy = sample(size = number_of_samples,
                               1:dim(posterior_chains_italy_sens)[1],
                               replace = T)

bootdata_sens = list()
for (i in 1:number_of_samples) {
  bootdata_sens[[i]] = posterior_chains_sens[posterior_index_sens[i],]
}


bootdata_italy_sens= list()
for (i in 1:number_of_samples) {
  bootdata_italy_sens[[i]] = posterior_chains_italy_sens[posterior_index_italy[i],]
}


posterior_samples_sens = bind_rows(bootdata_sens)
posterior_samples_italy_sens= bind_rows(bootdata_italy_sens)



write.csv(posterior_samples_sens, "100_posterior_samples_veneto_sens.csv")
write.csv(posterior_samples_italy_sens, "100_posterior_samples_italy_sens.csv")


################    Calculate mean and 95% CrI posterior estiamtes ################   


Veneto_post_est_sens  = cbind(
  apply (posterior_samples_sens, 2,  mean),
  apply(posterior_samples_sens, 2, quantile, probs = 0.025),
  apply(posterior_samples_sens, 2, quantile, probs = 0.975)
)
Italy_post_est_sens = cbind(
  apply (posterior_samples_italy_sens, 2,  mean),
  apply(posterior_samples_italy_sens, 2, quantile, probs = 0.025),
  apply(posterior_samples_italy_sens, 2, quantile, probs = 0.975)
)



################    Define column names for output ################   


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



################   Model 1: baseline scenario ################   


#### apply function across 1000 bootstrap samples of posterior estimates


##  TAKES 1-2 MINUTES TO RUN ##

modelfit1_sens = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior_samples_sens[i, ],
    n_difeq,
    n_pop ,
    n_recov ,
    n_months,
    gamma,
    sigma,
    phi_PCR,
    phi_Ag_sens,
    n_days,
    daily_Ag_i,
    daily_PCR_i,
    average_daily_vaccination_i,
    x_i_data,
    SEIR_model
  )
}, simplify = "array")



mean1_sens  = data.frame(apply(modelfit1_sens, c(1, 2),  mean))
lower1_sens = data.frame(apply(modelfit1_sens, c(1, 2), quantile, probs = 0.025))
upper1_sens = data.frame(apply(modelfit1_sens, c(1, 2), quantile, probs = 0.975))

modelfit1_df_sens= data.frame(mean1_sens, lower1_sens, upper1_sens)
colnames(modelfit1_df_sens) = col_names

### plot Veneto 

modelfit1_df_sens2 = modelfit1_df_sens %>%  
  select(! c(M_R0_M,A_R0_M,O_R0_M,Al_R0_M,pPCR_M,M_R0_L,A_R0_L,O_R0_L,Al_R0_L,pPCR_L,M_R0_U,A_R0_U,O_R0_U,Al_R0_U,pPCR_U)) 



modelfit1_df_veneto_inc = modelfit1_df_sens2 / (n_pop - n_recov) * 100000


variant_data2_sens = variant_data
variant_data2_sens[, 2:13] = variant_data2_sens[, 2:13] / (n_pop - n_recov) * 100000
modelfit1_vd_sens = cbind(variant_data2_sens, modelfit1_df_veneto_inc)


modelfit1_plot_sens = ggplot(modelfit1_vd_sens, aes(x = Date, y = M_Variant_M)) +
  geom_ribbon(aes(ymin = M_fit_L, ymax = M_fit_U),
              fill = "orange",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = A_fit_L, ymax = A_fit_U),
              fill = "sky blue",
              alpha = 0.3) +
  geom_ribbon(aes(ymin = O_fit_L, ymax = O_fit_U),
              fill = "pink",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = Al_fit_L, ymax = Al_fit_U),
              fill = "plum3",
              alpha = 0.5, 
              linetype = 2) +
  geom_line(aes(y = M_fit_M, color = "M234I-A376T model"), size = 1.5) +
  geom_line(aes(y = A_fit_M, color = "A220V model"), size = 1.5) +
  geom_line(aes(y = O_fit_M, color = "Other model"), size = 1.5) +
  geom_line(aes(y = Al_fit_M, color = "Alpha model"), size = 1.5) +
  
  geom_point(shape = 19, size = 3, (aes(y = A_Variant_M, color = "A220V data"))) +
  geom_errorbar(color = "blue", aes(ymin = A_Variant_L, ymax = A_Variant_U), size =1) +
  
  
  geom_point(shape = 19, size = 3, (aes(y = O_Variant_M, color = "Other data"))) +
  geom_errorbar(color = "red", aes(ymin = O_Variant_L, ymax = O_Variant_U), size =1) +
  
  
  geom_point(shape = 19, size = 3, (aes(y = Al_Variant_M, color = "Alpha data"))) +
  geom_errorbar(color = "mediumpurple4", aes(ymin = Al_Variant_L, ymax = Al_Variant_U), size =1) +
  
  geom_point(shape = 19, size = 3, (aes(color = "M234I-A376T data"))) +
  geom_errorbar(aes(ymin = M_Variant_L, ymax = M_Variant_U), size =1) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T data' = 'black',
      'M234I-A376T model' = 'orange',
      "A220V data" = "blue",
      'A220V model' = 'sky blue',
      'Other model' = 'pink',
      'Other data' = 'red',
      'Alpha model' = 'plum3',
      'Alpha data' = 'mediumpurple4'
    )
  ) +
  labs(x = " ",
       y = paste0("Reported incidence per 100,000 population")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    legend.position = c(0.2, 0.83),
    legend.key.height = unit(.5, 'cm'),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14),
    axis.text.x = element_text(face="bold")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  ggtitle("Veneto") +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 20))
modelfit1_plot_sens





## cumulative incidence 
cuminc1_sens  = data.frame(apply (modelfit1_sens, 2, colSums))

cuminc1_sens = cuminc1_sens %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)

mean_cum1_sens = data.frame(sapply(cuminc1_sens, mean))
lower_cum1_sens = data.frame(sapply(cuminc1_sens, quantile, probs = 0.025))
upper_cum1_sens = data.frame(sapply(cuminc1_sens, quantile, probs = 0.975))


cuminc1_df_sens = data.frame(t(cbind(mean_cum1_sens, lower_cum1_sens, upper_cum1_sens)))
rownames(cuminc1_df_sens) = c("mean", "lower", "upper")

################    Model 1 Italy : baseline scenario ################   




#### apply function across 1000 bootstrap samples of posterior estimates


##  TAKES 1-2 MINUTES TO RUN ##

modelfit1_italy_sens  = sapply(1:number_of_samples, function(i) {
  replicate_rstan_fixed(
    posterior = posterior_samples_italy_sens[i, ],
    n_difeq = n_difeq,
    n_pop = n_pop_italy ,
    n_recov =  n_recov_italy ,
    n_months =  n_months_italy,
    gamma =   gamma,
    sigma =   sigma,
    phi_PCR,
    phi_Ag_sens,
    n_days =  n_days_italy,
    daily_Ag_i =   daily_Ag_i_italy,
    daily_PCR_i =   daily_PCR_i_italy,
    average_daily_vaccination_i =   average_daily_vaccination_i_italy,
    x_i_data =     x_i_data_italy,
    SEIR_model =    SEIR_model
  )
}, simplify = "array")



mean1_italy_sens  = data.frame(apply (modelfit1_italy_sens, c(1, 2),  mean))
lower1_italy_sens = data.frame(apply(modelfit1_italy_sens, c(1, 2), quantile, probs = 0.025))
upper1_italy_sens = data.frame(apply(modelfit1_italy_sens, c(1, 2), quantile, probs = 0.975))

modelfit1_df_italy_sens = data.frame(mean1_italy_sens, lower1_italy_sens, upper1_italy_sens)
colnames(modelfit1_df_italy_sens) = col_names


modelfit1_df_italy_sens2 = modelfit1_df_italy_sens %>%  
  select(! c(M_R0_M,A_R0_M,O_R0_M,Al_R0_M,pPCR_M,M_R0_L,A_R0_L,O_R0_L,Al_R0_L,pPCR_L,M_R0_U,A_R0_U,O_R0_U,Al_R0_U,pPCR_U)) 



modelfit1_df_italy_inc = modelfit1_df_italy_sens2 / (n_pop_italy - n_recov_italy) * 100000

variant_data_italy2_Sens = variant_data_italy
variant_data_italy2_Sens[, 2:13] = variant_data_italy[, 2:13] / (n_pop_italy - n_recov_italy) * 100000


modelfit1_vd_italy_Sens = cbind(variant_data_italy2_Sens, modelfit1_df_italy_inc)

## plot Italy 

modelfit1_plot_italy_sens = ggplot(modelfit1_vd_italy_Sens, aes(x = Date, y = M_Variant_M)) +
  geom_ribbon(aes(ymin = M_fit_L, ymax = M_fit_U),
              fill = "orange",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = A_fit_L, ymax = A_fit_U),
              fill = "sky blue",
              alpha = 0.3) +
  geom_ribbon(aes(ymin = O_fit_L, ymax = O_fit_U),
              fill = "pink",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = Al_fit_L, ymax = Al_fit_U),
              fill = "plum3",
              alpha = 0.5, 
              linetype = 2) +
  geom_line(aes(y = M_fit_M, color = "M234I-A376T model"), size = 1.5) +
  geom_line(aes(y = A_fit_M, color = "A220V model"), size = 1.5) +
  geom_line(aes(y = O_fit_M, color = "Other model"), size = 1.5) +
  geom_line(aes(y = Al_fit_M, color = "Alpha model"), size = 1.5) +
  
  geom_point(shape = 19, size = 3, (aes(y = A_Variant_M, color = "A220V data"))) +
  geom_errorbar(color = "blue", aes(ymin = A_Variant_L, ymax = A_Variant_U), size =1) +
  
  
  geom_point(shape = 19, size = 3, (aes(y = O_Variant_M, color = "Other data"))) +
  geom_errorbar(color = "red", aes(ymin = O_Variant_L, ymax = O_Variant_U), size =1) +
  
  
  geom_point(shape = 19, size = 3, (aes(y = Al_Variant_M, color = "Alpha data"))) +
  geom_errorbar(color = "mediumpurple4", aes(ymin = Al_Variant_L, ymax = Al_Variant_U), size =1) +
  
  geom_point(shape = 19, size = 3, (aes(color = "M234I-A376T data"))) +
  geom_errorbar(aes(ymin = M_Variant_L, ymax = M_Variant_U), size =1) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T data' = 'black',
      'M234I-A376T model' = 'orange',
      "A220V data" = "blue",
      'A220V model' = 'sky blue',
      'Other model' = 'pink',
      'Other data' = 'red',
      'Alpha model' = 'plum3',
      'Alpha data' = 'mediumpurple4'
    )
  ) +
  labs(x = "", y = paste0(" ")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    axis.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face="bold")
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  ggtitle("Rest of Italy") +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 20))

modelfit1_plot_italy_sens

## cumulative incidence 
cuminc1_it_sens  = data.frame(apply (modelfit1_italy_sens, 2, colSums))

cuminc1_it_sens = cuminc1_it_sens %>%  
  mutate(Rep_M = fit_SEIR_M/  true_SEIR_M , 
         Rep_A = fit_SEIR_A/  true_SEIR_A , 
         Rep_O = fit_SEIR_O/  true_SEIR_O , 
         Rep_Al = fit_SEIR_Al/  true_SEIR_Al ,
         Rep_tot = total_reported_incidence/total_incidence)

mean_cum1_sens_it = data.frame(sapply(cuminc1_it_sens, mean))
lower_cum1_sens_it = data.frame(sapply(cuminc1_it_sens, quantile, probs = 0.025))
upper_cum1_sens_it = data.frame(sapply(cuminc1_it_sens, quantile, probs = 0.975))


cuminc1_df_sens_it = data.frame(t(cbind(mean_cum1_sens_it, lower_cum1_sens_it, upper_cum1_sens_it)))
rownames(cuminc1_df_sens_it) = c("mean", "lower", "upper")



################  Plot proportion of incidence reported in Veneto and Italy ################  

proportion_reported_plot_sens = cuminc1_df_sens %>%  
  bind_rows(cuminc1_df_sens_it) %>%  
  mutate(X = rep(rownames(cuminc1_df_sens),2)) %>% 
  select(X, Rep_M,     Rep_A,     Rep_O ,   Rep_Al  , Rep_tot) %>% 
  rename(Stat = X) %>% 
  mutate(Country = factor(rep(c("Veneto", "Rest of Italy"), each = 3), levels = c("Veneto", "Rest of Italy"))) %>%  
  pivot_longer(cols= -c(Country,Stat) ) %>%  
  mutate(Variant = factor (rep ( c( "M234I-A376T", "A220V", "Other", "Alpha", "Total model"), 6)),
         value = value * 100) %>% 
  pivot_wider(id_cols = c(Country,Variant), names_from = Stat,values_from = value)  %>%  
  ggplot(aes(x = Country, y = mean)) +
  geom_point( shape = 19, size = 4,aes(color = rep(c("M234I-A376T", "A220V", "Other", "Alpha", "Total model"),2)), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = rep(c("M234I-A376T", "A220V", "Other", "Alpha", "Total model"),2)),
                position = position_dodge(width = 0.5), width =  0.4, size = 1) +
  labs(x = " ", y = paste0("Percentage of cumulative incidence reported")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c(0.12,0.92),
    axis.text.x = element_text(face="bold")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T' = 'orange',
      'A220V' = 'sky blue',
      'Other' = 'pink',
      'Alpha' = 'plum3',
      "Total model" = " darkgrey"),
    breaks = c('Total model'))






################   Save output ################  





incidence__sens_fig = plot_grid(
  modelfit1_plot_sens,
  modelfit1_plot_italy_sens,
  proportion_reported_plot_sens,
  ncol = 3,
  labels = c("a", "b", "c"),
  align = "h"
)




ggsave(
  plot = incidence__sens_fig,
  filename = "FigS10.jpg",
  height = 20,
  width = 60,
  units = "cm",
  dpi = 300
)




