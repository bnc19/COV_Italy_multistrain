
################################         READ ME         ################################   
# 
# This script fits a multivariant model to epidemiological and genomic data from Veneto 
# or the rest of Italy. The model can be run on a high performance cluster (HPC) or
# locally. Currently, each model runs 3 chains for  200 iterations, discarding the first 
# 100 as burnin. Initial values are close to the posterior means. This is to demo the
# models ~ quickly (12 hours on a normal computer). In order to obtain the results 
# provided in the manuscript, 3 chains with more diffuse starting points were run for
# 2000 iterations, discarding the first 1000 iterations as burnin. This takes ~48-72 
# hours to run on the HPC. 
#
################################                         ################################   

# create folders 

dir.create("model_fitting/Results")
dir.create("model_fitting/Results/Italy")
dir.create("model_fitting/Results/Veneto")


################################    Set up to run on HPC ################################   
# # rm(list = ls())
# # root = "Q:/COV_Italy_multistrain/model_fitting"
# # setwd("Q:/COV_Italy_multistrain/model_fitting")
# 
# 
# ctx <- context::context_save(root, packages=c("dplyr", "ggplot2", "rstan", "bayesplot", "Hmisc", "BH", "RcppEigen"),
#                              sources=c("R/run_testing_rstan_sens.R"))
# 
# # Set up to run with multiple cores on clusters
# config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")
# 
# 
# # Create a queue within the context (/environment)
# obj <- didehpc::queue_didehpc(ctx, config)

################################    Run NegBin Italy on cluster  ################################    


# I_V1 =  obj$enqueue(
#   run_SEIR_stan_Italy_sens (
#     A_data = read.csv("data/Dataset_Italy_A_v5.csv")$Freq_new,  # Genomic data 
#     M_data = read.csv("data/Dataset_Italy_M_v5.csv")$Freq_new,
#     O_data = read.csv("data/Dataset_Italy_O_v1.csv")$Freq_new,
#     Al_data = read.csv("data/Dataset_Italy_Alpha_v1.csv") $Freq_new,
#     n_seq = read.csv("data/Dataset_Italy_A_v5.csv")$TotSeq_new,
#     n_difeq = 11,   
#     n_pop = (59257566 - 4847026 ),
#     n_recov = (1482377 - 93401),
#     index_M = 6:14,
#     index_A = 3:14,
#     index_O =  3:9,
#     index_Al = 9:14,
#     modelPath = "model/new_neg_bin_sens.stan",  
#     Location = "Italy",
#     n_chains =3,
#     n_warmups =100,
#     n_iter = 200,
#     n_thin = 1,
#     pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
#              , "I0[1]", "I0[2]", "I0[3]", "I0[4]" , "k" ),
#     
#     filePath = "Results/Italy",
#     ini_1_SEIR = function(){
#       list(  beta = replicate(4,runif(1,1,3)),   
#              I0 = replicate(4, runif(1, 1,10000)),
#              omega = replicate(2,runif(1,0.1,0.9)),
#              rho = runif(1,0.1,0.9),
#              k = runif(1,0.01,5)
#       )},
#     average_daily_vaccination = read.csv("data/Vac_Italy_For_Month.csv")$prop_vac,
#     average_daily_reported_incidence = read.csv("data/Dataset_Italy_A_v5.csv")$new_reported_cases_daily_new,
#     daily_reported_incidence = read.csv("data/dailyReportedIncidence_italy.csv")$new_case ,
#     daily_PCR = round(read.csv("data/Italy_daily_test_data.csv")$pcr_daily_average ) ,
#     daily_Ag= round(read.csv("data/Italy_daily_test_data.csv")$antigen_daily_average       ) ,
#     monthly_PCR = round(read.csv("data/Italy_monthly_test_data.csv")$pcr_daily_average),
#     monthly_Ag = round(read.csv("data/Italy_monthly_test_data.csv")$antigen_daily_average) ,
#     start_date = "01-05-2020",
#     end_date = "31-05-2021",
#     time_intervention = c("07-11-2020" , "15-03-2021"),
#     time_seed_alpha = "01-11-2020",
#     time_seed_M = "01-08-2020",
#     sigma = 1 / 5.1,   
#     gamma = 1 / 2.1 ,
#     phi_PCR = 0.920,
#     phi_Ag = 0.875, # 0.643 - change antigen sensitivity 
#     prior_seed_mean = 1000,
#     prior_seed_sd = 1000
#   )
# )
# 
# 
# ################  Run NegBin Veneto on cluster ################  
# 
# V_V1 = obj$enqueue(
#   run_SEIR_stan_Italy_sens (
#     A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq ,
#     M_data = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq ,
#     O_data = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq ,
#     Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq ,
#     n_seq = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq,
#     
#     n_difeq = 11,
#     n_pop = 4847026,
#     n_recov = 93401,
#     
#     
#     modelPath = "model/new_neg_bin_sens.stan",
#     Location = "Veneto",
#     n_chains =3,
#     n_warmups =100,
#     n_iter = 200,
#     n_thin = 1,
#     
#     pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
#              , "I0[1]", "I0[2]", "I0[3]", "I0[4]"),
#     
#     filePath = "Results/Veneto" ,
#     ini_1_SEIR = function(){
#       list(  beta = replicate(4,runif(1,1,3)),
#              I0 = replicate(4, runif(1, 1,1000)),
#              omega = replicate(2,runif(1,0.1,0.9)),
#              rho = runif(1,0.1,0.9),
#              k = runif(1,0.01,2)
#       )},
#     average_daily_reported_incidence = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily,
#     daily_reported_incidence = read.csv("data/dailyReportedIncidence_veneto.csv")$new_case ,
#     daily_PCR = round(read.csv("data/Veneto_daily_test_data.csv")$pcr_daily),
#     daily_Ag= round(read.csv("data/Veneto_daily_test_data.csv")$antigen_daily),
#     monthly_PCR = round(read.csv("data/Veneto_monthly_test_data.csv")$pcr_daily_average),
#     monthly_Ag = round(read.csv("data/Veneto_monthly_test_data.csv")$antigen_daily_average),
#     average_daily_vaccination =  read.csv("data/Vac_Veneto_For_Month.csv")$prop_vac,
#     index_M = 8:14,
#     index_A = c(5,7:14),
#     index_O = c(5,7:9),
#     index_Al = 9:14 ,
#     start_date = "01-07-2020",
#     end_date = "31-05-2021",
#     time_intervention = c("15-11-2020" , "15-03-2021"),
#     time_seed_alpha = "01-11-2020",
#     time_seed_M = "01-10-2020",
#     sigma = 1 / 5.1,
#     gamma = 1 / 2.1 ,
#     phi_PCR = 0.920,
#     phi_Ag = 0.875,  # 0.643 - change antigen sensitivity 
#     prior_seed_mean = 1,
#     prior_seed_sd = 1000
#   )
# )
# 
# 



################  Run NegBin Italy Locally  ################  


I_V1 = run_SEIR_stan_Italy_sens (
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
  modelPath = "model/new_neg_bin_sens.stan",
  Location = "Italy",
  n_chains =3,
  n_warmups =100,
  n_iter = 200,
  n_thin = 1,
  pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
           , "I0[1]", "I0[2]", "I0[3]", "I0[4]" , "k" ),
  
  filePath = "Results/Italy",
  ini_1_SEIR = function(){
    list(  beta = replicate(4,runif(1,1,3)),
           I0 = replicate(4, runif(1, 1,10000)),
           omega = replicate(2,runif(1,0.1,0.9)),
           rho = runif(1,0.1,0.9),
           k = runif(1,0.01,5)
    )},
  average_daily_vaccination = read.csv("data/Vac_Italy_For_Month.csv")$prop_vac,
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
  phi_Ag = 0.875,  # 0.643 - change antigen sensitivity 
  prior_seed_mean = 1000,
  prior_seed_sd = 1000
)


################  Run NegBin Veneto Locally  ################  

V_V1 = run_SEIR_stan_Italy_sens (
  A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq ,
  M_data = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq ,
  O_data = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq ,
  Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq ,
  n_seq = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq,
  
  n_difeq = 11,
  n_pop = 4847026,
  n_recov = 93401,
  
  
  modelPath = "model/new_neg_bin_sens.stan",
  Location = "Veneto",
  n_chains =3,
  n_warmups =100,
  n_iter = 200,
  n_thin = 1,
  
  pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
           , "I0[1]", "I0[2]", "I0[3]", "I0[4]"),
  
  filePath = "Results/Veneto" ,
  ini_1_SEIR = function(){
    list(  beta = replicate(4,runif(1,1,3)),
           I0 = replicate(4, runif(1, 1,1000)),
           omega = replicate(2,runif(1,0.1,0.9)),
           rho = runif(1,0.1,0.9),
           k = runif(1,0.01,2)
    )},
  average_daily_reported_incidence = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily,
  daily_reported_incidence = read.csv("data/dailyReportedIncidence_veneto.csv")$new_case ,
  daily_PCR = round(read.csv("data/Veneto_daily_test_data.csv")$pcr_daily),
  daily_Ag= round(read.csv("data/Veneto_daily_test_data.csv")$antigen_daily),
  monthly_PCR = round(read.csv("data/Veneto_monthly_test_data.csv")$pcr_daily_average),
  monthly_Ag = round(read.csv("data/Veneto_monthly_test_data.csv")$antigen_daily_average),
  average_daily_vaccination =  read.csv("data/Vac_Veneto_For_Month.csv")$prop_vac,
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
  phi_Ag = 0.875,  # 0.643 - change antigen sensitivity 
  prior_seed_mean = 1,
  prior_seed_sd = 1000
)


