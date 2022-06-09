source("R/run_model_fitting.R")

#  Set up to run on HPC  -------------------------------------------------------
rm(list = ls())
root = "Q:/COV_Italy_multistrain/model_fitting"
setwd("Q:/COV_Italy_multistrain/model_fitting")


ctx <- context::context_save(root, packages=c("tidyverse", "rstan", "bayesplot", "Hmisc", 
                                              "BH", "RcppEigen", "png", "knitr", "rstudioapi",
                                              "StanHeaders", "cowplot"),
                             sources=c("R/run_model_fitting.R"))

# Set up to run with multiple cores on clusters
config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")


# Create a queue within the context (/environment)
obj <- didehpc::queue_didehpc(ctx, config)

##################### Interventions middle of month ############################


# Run model with testing of symptomatic only -----------------------------------
test_symp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/symp_test",
    model_path = "model/est_test_symp.stan"
  )
)

test_symp$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_test",
    model_path = "model/est_test_asymp.stan"
  )
)

test_asymp$log()

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test",
    model_path = "model/est_test_asymp_and_symp.stan"
  )
)

test_asymp_and_symp$log()


##################### Interventions first of month ############################

# Run model with testing of symptomatic only -----------------------------------
test_symp_int_1st =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/symp_test_int_1st",
    model_path = "model/est_test_symp.stan",
    time_intervention_it = c("01-11-2020", "01-03-2021"),
    time_intervention_veneto  = c("01-11-2020", "01-03-2021")
  )
)
test_symp_int_1st$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp_int_1st =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/test_asymp_int_1st",
    model_path = "model/est_test_asymp.stan",
    time_intervention_it = c("01-11-2020", "01-03-2021"),
    time_intervention_veneto  = c("01-11-2020", "01-03-2021")
    
  )
)

test_asymp_int_1st$log()

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_int_1st =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test_int_1st",
    model_path = "model/est_test_asymp_and_symp.stan",
    time_intervention_it = c("01-11-2020", "01-03-2021"),
    time_intervention_veneto  = c("01-11-2020", "01-03-2021")
    
  )
)

test_asymp_and_symp_int_1st$log()




##################### Interventions last of month ############################

# Run model with testing of symptomatic only -----------------------------------
test_symp_int_last =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/symp_test_int_last",
    model_path = "model/est_test_symp.stan",
    time_intervention_it = c("30-11-2020", "31-03-2021"),
    time_intervention_veneto  = c("30-11-2020", "31-03-2021")
  )
)
test_symp_int_last$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp_int_last =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/test_asymp_int_last",
    model_path = "model/est_test_asymp.stan",
    time_intervention_it = c("30-11-2020", "31-03-2021"),
    time_intervention_veneto  = c("30-11-2020", "31-03-2021")
    
  )
)

test_asymp_int_last$log()

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_int_last =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test_int_last",
    model_path = "model/est_test_asymp_and_symp.stan",
    time_intervention_it = c("30-11-2020", "31-03-2021"),
    time_intervention_veneto  = c("30-11-2020", "31-03-2021")
    
  )
)

test_asymp_and_symp_int_last$log()
