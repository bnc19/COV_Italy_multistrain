
################  format data to run counterfactuals ################  

################  Set up ################  

# rm(list = ls())
#  setwd("Q:/COV_Italy_multistrain/model_fitting")


library(tidyverse)
library(Hmisc)
source("counterfactuals/R/Format_data_for_SEIQR_model.R")  # function to format external data based on seeding 
source("counterfactuals/R/Format_variant_data.R")



################   Import and define external data for Italy ################  

susceptible_italy =  59257566 - 4847026
recovered_italy =  1482377 - 93401


# external data on:

# per capita vaccination
# daily reported incidence
# daily number of PCR tests
# daily number of Ag tests
# daily number of deaths 


# convert dates into date format and format daily data to match 

formated_data_italy  = format_data_for_SEIQR_model(
  read.csv("model_fitting/data/dailyReportedIncidence_italy.csv")$new_case,
  c(0,0,0,0, 0,0,0,0,0,0.0003,0.0005, 0.0009 ,0.0014, 0.002),
  round(read.csv("model_fitting/data/Italy_monthly_test_data.csv")$pcr_daily_average),
  round(read.csv("model_fitting/data/Italy_monthly_test_data.csv")$antigen_daily_average),
  round(read.csv("model_fitting/data/Italy_daily_test_data.csv")$pcr_daily) ,
  round(read.csv("model_fitting/data/Italy_daily_test_data.csv")$antigen_daily),
  start_date = "01-05-2020",
  end_date = "31-05-2021",
  time_intervention = c("07-11-2020" , "15-03-2021"),
  time_seed_alpha = "01-11-2020",
  time_seed_M = "01-08-2020" )


Italy_variant_data = Format_variant_data(
  index_M = 6:14,
  index_A = 3:14,
  index_O =  3:9,
  index_Al = 9:14,
  start_date = "01-05-2020",
  end_date = "31-05-2021",
  A_data = read.csv("model_fitting/data/Dataset_Italy_A_v5.csv")$Freq_new,
  M_data = read.csv("model_fitting/data/Dataset_Italy_M_v5.csv")$Freq_new,
  O_data = read.csv("model_fitting/data/Dataset_Italy_O_v1.csv")$Freq_new,
  Al_data = read.csv("model_fitting/data/Dataset_Italy_Alpha_v1.csv") $Freq_new,
  n_seq = read.csv("model_fitting/data/Dataset_Italy_A_v5.csv")$TotSeq_new,
  average_daily_reported_incidence = read.csv("model_fitting/data/Dataset_Italy_A_v5.csv")$new_reported_cases_daily_new)

x_i_data_italy = c(
  formated_data_italy[[2]]$n_months,
  formated_data_italy[[4]],
  formated_data_italy[[3]]$monthly_PCR_i,
  formated_data_italy[[3]]$monthly_Ag_i,
  susceptible_italy - recovered_italy,
  recovered_italy,
  formated_data_italy[[2]]$n_days,
  formated_data_italy[[2]]$time_switch1,
  formated_data_italy[[2]]$time_switch2,
  formated_data_italy[[2]]$seed_alpha, 
  formated_data_italy[[2]]$seed_M
)



write.csv(Italy_variant_data, "counterfactuals/Italy/Italy_variant_data.csv")
write.csv(x_i_data_italy, "counterfactuals/Italy/x_i_data_italy.csv")
write.csv(formated_data_italy[[3]]$average_monthly_vaccination_i, "counterfactuals/Italy/average_daily_vaccination_i_italy.csv")
write.csv(formated_data_italy[[1]]$daily_Ag_i, "counterfactuals/Italy/daily_Ag_i_italy.csv")
write.csv(formated_data_italy[[1]]$daily_PCR_i, "counterfactuals/Italy/daily_PCR_i_italy.csv")




#############################################################################
# Import and define external data for Veneto 
#############################################################################

susceptible_ven =  4847026
recovered_ven =  93401


# external data on:

# per capita vaccination
# daily reported incidence
# daily number of PCR tests
# daily number of Ag tests
# daily number of deaths 



# convert dates into date format and format daily data to match 

formated_data_Veneto  = format_data_for_SEIQR_model(
  read.csv("model_fitting/data/dailyReportedIncidence_veneto.csv")$new_case,
  c(0,0,0,0, 0,0,0,0,0,0.00043, 0.00039,0.00096,0.0018, 0.0028),
  round(read.csv("model_fitting/data/Veneto_monthly_test_data.csv")$pcr_daily_average),
  round(read.csv("model_fitting/data/Veneto_monthly_test_data.csv")$antigen_daily_average),
  round(read.csv("model_fitting/data/Veneto_daily_test_data.csv")$pcr_daily) ,
  round(read.csv("model_fitting/data/Veneto_daily_test_data.csv")$antigen_daily),
  start_date = "01-07-2020",
  end_date = "31-05-2021",
  time_intervention = c("15-11-2020" , "15-03-2021"),
  time_seed_alpha = "01-11-2020",
  time_seed_M = "01-10-2020" )


Veneto_variant_data = Format_variant_data(
index_M = 8:14,
index_A = c(5,7:14),
index_O = c(5,7:9),
index_Al = 9:14 ,
start_date = "01-07-2020",
end_date = "31-05-2021",
A_data = read.csv("model_fitting/data/Dataset_Veneto_A_v5.csv")$Freq ,
M_data = read.csv("model_fitting/data/Dataset_Veneto_M_v5.csv")$Freq ,
O_data = read.csv("model_fitting/data/Dataset_Veneto_O_v1.csv")$Freq ,
Al_data = read.csv("model_fitting/data/Dataset_Veneto_Alpha_v1.csv")$Freq ,
n_seq = read.csv("model_fitting/data/Dataset_Veneto_A_v5.csv")$TotSeq,
average_daily_reported_incidence = read.csv("model_fitting/data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily)

x_i_data_Veneto = c(
  formated_data_Veneto[[2]]$n_months,
  formated_data_Veneto[[4]],
  formated_data_Veneto[[3]]$monthly_PCR_i,
  formated_data_Veneto[[3]]$monthly_Ag_i,
  susceptible_ven - recovered_ven,
  recovered_ven,
  formated_data_Veneto[[2]]$n_days,
  formated_data_Veneto[[2]]$time_switch1,
  formated_data_Veneto[[2]]$time_switch2,
  formated_data_Veneto[[2]]$seed_alpha, 
  formated_data_Veneto[[2]]$seed_M
)


################ italy testing in veneto #######################

formated_data_Veneto_I_test  = format_data_for_SEIQR_model(
  read.csv("model_fitting/data/dailyReportedIncidence_veneto.csv")$new_case,
  c(0,0,0,0, 0,0,0,0,0,0.00043, 0.00039,0.00096,0.0018, 0.0028),
  round(read.csv("model_fitting/data/Italy_monthly_test_data.csv")$pcr_daily_average),
  round(read.csv("model_fitting/data/Italy_monthly_test_data.csv")$antigen_daily_average),
  round(read.csv("model_fitting/data/Italy_daily_test_data.csv")$pcr_daily) ,
  round(read.csv("model_fitting/data/Italy_daily_test_data.csv")$antigen_daily),
  start_date = "01-07-2020",
  end_date = "31-05-2021",
  time_intervention = c("15-11-2020" , "15-03-2021"),
  time_seed_alpha = "01-11-2020",
  time_seed_M = "01-10-2020" )


x_i_data_Veneto_italy_test = c(
  formated_data_Veneto[[2]]$n_months,
  formated_data_Veneto[[4]],
  formated_data_Veneto_I_test[[3]]$monthly_PCR_i,
  formated_data_Veneto_I_test[[3]]$monthly_Ag_i,
  susceptible_ven - recovered_ven,
  recovered_ven,
  formated_data_Veneto[[2]]$n_days,
  formated_data_Veneto[[2]]$time_switch1,
  formated_data_Veneto[[2]]$time_switch2,
  formated_data_Veneto[[2]]$seed_alpha, 
  formated_data_Veneto[[2]]$seed_M
)



write.csv(Veneto_variant_data, "counterfactuals/Veneto/Veneto_variant_data.csv")
write.csv(x_i_data_Veneto, "counterfactuals/Veneto/x_i_data.csv")
write.csv(formated_data_Veneto[[3]]$average_monthly_vaccination_i, "counterfactuals/Veneto/average_daily_vaccination_i.csv")
write.csv(formated_data_Veneto[[1]]$daily_Ag_i, "counterfactuals/Veneto/daily_Ag_i.csv")
write.csv(formated_data_Veneto[[1]]$daily_PCR_i, "counterfactuals/Veneto/daily_PCR_i.csv")
write.csv(x_i_data_Veneto_italy_test, "counterfactuals/Veneto/x_i_data_Veneto_italy_test.csv")
write.csv(formated_data_Veneto_I_test[[1]]$daily_Ag_i, "counterfactuals/Veneto/daily_Ag_i_Itest.csv")
write.csv(formated_data_Veneto_I_test[[1]]$daily_PCR_i, "counterfactuals/Veneto/daily_PCR_i_Itest.csv")

