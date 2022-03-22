
#### Set up 
# rm(list = ls())
root <- "Q:/Antigen_testing"
setwd("Q:/Antigen_testing")



ctx <- context::context_save(root, packages=c("dplyr", "ggplot2", "rstan", "bayesplot", "Hmisc", "BH", "RcppEigen"), 
                             sources=c("R/run_testing_rstan_sens.R"))




# Set up to run with multiple cores on clusters
config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")


# Create a queue within the context (/environment)
obj <- didehpc::queue_didehpc(ctx, config) 



obj$task_get("8f1858559470fd583bf9c9d916c5cd6a")$log()


#########################################################
###### Run NegBin Italy on cluster Sensitivity  #########
#########################################################

I_V1 =  obj$enqueue(
  run_SEIR_stan_Italy_sens (
    A_data = read.csv("data/Dataset_Italy_A_v5.csv")$Freq_new,
    M_data = read.csv("data/Dataset_Italy_M_v5.csv")$Freq_new,
    O_data = read.csv("data/Dataset_Italy_O_v1.csv")$Freq_new,
    Al_data = read.csv("data/Dataset_Italy_Alpha_v1.csv") $Freq_new,
    n_seq = read.csv("data/Dataset_Italy_A_v5.csv")$TotSeq_new,
    n_difeq = 11,
    n_pop = (59257566 - 4847026 ),
    n_recov = (1482377 - 93401),
    index_M = 6:14,
    index_A = 3:14,
    index_O =  3:9,
    index_Al = 9:14,
    modelPath = "new_neg_bin_sens.stan",
    Location = "Italy",
    n_chains =3,
    n_warmups =1000,
    n_iter = 2000,
    n_thin = 1,
    pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
             , "I0[1]", "I0[2]", "I0[3]", "I0[4]" , "k"
             
    ),
    
    filePath = "Results/Italy/0.875/I_V2/LAD_99",
    ini_1_SEIR = function(){
      list(   beta = replicate(4,runif(1,2,4)),
              I0 = replicate(4, runif(1, 1,10000)),
              omega = replicate(2,runif(1,0.5,.8)),
              rho = runif(1,0.5,.8),
              k = runif(1,0.01,0.5)
      )},
    average_daily_vaccination =  c(0,0,0,0, 0,0,0,0,0,0.0003,0.0005, 0.0009 ,0.0014, 0.0029),
    
    average_daily_reported_incidence = read.csv("data/Dataset_Italy_A_v5.csv")$new_reported_cases_daily_new,
    daily_reported_incidence = read.csv("data/dailyReportedIncidence_italy.csv")$new_case ,
    
    daily_PCR = round(read.csv("data/Italy_daily_test_data.csv")$pcr_daily_average ) ,
    daily_Ag= round(read.csv("data/Italy_daily_test_data.csv")$antigen_daily_average       ) ,
    monthly_PCR = round(read.csv("data/Italy_monthly_test_data.csv")$pcr_daily_average),
    monthly_Ag = round(read.csv("data/Italy_monthly_test_data.csv")$antigen_daily_average) ,
    
    
    start_date = "01-05-2020",
    end_date = "31-05-2021",
    time_intervention = c("07-11-2020" , "15-03-2021"),
    time_seed_alpha = "01-11-2020",
    time_seed_M = "01-08-2020",
    
    sigma = 1 / 5.1,
    gamma = 1 / 2.1 ,
    phi_PCR = 0.920,
    phi_Ag = 0.875,
    NB = TRUE,
    AD = 0.99 , 
    prior_seed_mean = 1000,
    prior_seed_sd = 1000,
    rho_a = 2,
    rho_b = 2,
    deaths = read.csv("data/Italy_deaths.csv")
  )
)

obj$task_get("aaf62da4421f270e39693e94e8626583")$log()


# 2000, higher beta AD 0.99 b00382747d4bf1fd5d47a58d7f389e85
# same but 0.95 6760ed9627166d6c0b1e731df4707bf8
#  0.99 3df7646f1d27ef4af8a4b02babcdf37e

###################

x = x$result()

library(loo)
library(rstanarm)
library(bayesplot )
log_lik <- extract_log_lik(x, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik), cores = 2)
loov2 <- loo(log_lik, r_eff = r_eff, cores = 2)

if (any(pareto_k_values(loov2) > 0.7)) {
  loo2 <- loo(x, save_psis = TRUE, k_threshold = 0.7)
}

#########################################################
######## Run NegBin Veneto on cluster sensitivity #######
#########################################################



V_V1 = obj$enqueue(
  run_SEIR_stan_Italy_sens (
    A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq ,
    M_data = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq ,
    O_data = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq ,
    Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq ,
    n_seq = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq,
    
    n_difeq = 11,
    n_pop = 4847026,
    n_recov = 93401,
    
    
    modelPath = "new_neg_bin_sens.stan",
    Location = "Veneto",
    n_chains =3,
    n_warmups =2000,
    n_iter = 4000,
    n_thin = 1,
    
    pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
             , "I0[1]", "I0[2]", "I0[3]", "I0[4]"),
    
    
    filePath = "Results/Veneto2/far_init" ,
    ini_1_SEIR = function(){
      list(  beta = replicate(4,runif(1,0,10)),
             I0 = replicate(4, runif(1, 1,10000)),
             omega = replicate(2,runif(1,0,1)),
             rho = runif(1,0,1),
             k = runif(1,0.01,10)
      )},
    average_daily_reported_incidence = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily,
    daily_reported_incidence = read.csv("data/dailyReportedIncidence.csv")$new_case ,
    daily_PCR = round(read.csv("data/Veneto_daily_test_data.csv")$pcr_daily),
    daily_Ag= round(read.csv("data/Veneto_daily_test_data.csv")$antigen_daily),
    monthly_PCR = round(read.csv("data/Veneto_monthly_test_data.csv")$pcr_daily_average),
    monthly_Ag = round(read.csv("data/Veneto_monthly_test_data.csv")$antigen_daily_average),
    average_daily_vaccination =  c(0,0,0,0, 0,0,0,0,0,0.00043, 0.00039,0.00096,0.0018, 0.0028),
    index_M = 8:14,
    index_A = c(5,7:14),
    index_O = c(5,7:9),
    index_Al = 9:14 ,
    start_date = "01-07-2020",
    end_date = "31-05-2021",
    time_intervention = c("15-11-2020" , "15-03-2021"),
    time_seed_alpha = "01-11-2020",
    time_seed_M = "01-10-2020",
    sigma = 1 / 5.1,
    gamma = 1 / 2.1 ,
    phi_PCR = 0.920,
    phi_Ag = 0.643,
    NB = TRUE,
    AD = 0.95,
    prior_seed_mean = 1,
    prior_seed_sd = 1000,
    rho_a = 2,
    rho_b = 2,
    deaths = read.csv("data/Veneto_deaths.csv")
    
  )
  
)


obj$task_get("a8d7d346d52885b7a3f768f00e184a64")$log()



########### 0.875 closer start values ###############


# 0.9 AD   aa76da48165b1de0b1bf107dc199604c

# 0.95 AD  a8d7d346d52885b7a3f768f00e184a64

# 0.99 AD  2d929bdac3ae207851f59a63dbb65a2c

################ far away initial values #################
# aaf62da4421f270e39693e94e8626583



# 0.643 4000 1eeadafed23649a4a733fac370073001



library(loo)
log_lik_2 <- extract_log_lik(V_V1L, merge_chains = FALSE)
r_eff_2 <- relative_eff(exp(log_lik_2))
loo_2 <- loo(log_lik_2, r_eff = r_eff_2, cores = 2)
print(loo_2)

