

#  Set up to run on HPC  -------------------------------------------------------
rm(list = ls())
root = "Q:/COV_Italy_multistrain/model_fitting"
setwd("Q:/COV_Italy_multistrain/model_fitting")


ctx <- context::context_save(root, packages=c("tidyverse", "rstan", "bayesplot", "Hmisc", 
                                              "BH", "RcppEigen", "png", "knitr", "rstudioapi",
                                              "StanHeaders", "cowplot", "loo"),
                             sources=c("R/run_model_fitting.R"))

# Set up to run with multiple cores on clusters
config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")


# Create a queue within the context (/environment)
obj <- didehpc::queue_didehpc(ctx, config)

# 
source("R/run_model_fitting.R")

# Global parameters ------------------------------------------------------------

n_iter = 1000
n_warmups = round(n_iter/2)


############################# Fit models to incidence ##########################


# Run model with testing of symptomatic only -----------------------------------
test_symp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results",
    model_path = "model/est_test_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups
  )
)

test_symp$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_test",
    model_path = "model/est_test_asymp.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups
  )
)


test_asymp$log()


# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test",
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups
  )
)


test_asymp_and_symp$log()



############################# Fit models to prevalence ##########################

# Add Tau 
pars = c("lp__",
         "beta[1]",
         "beta[2]",
         "beta[3]",
         "beta[4]",
         "rho_it" ,
         "rho_ven" ,
         "omega[1]",
         "omega[2]",
         "I0_ven[1]",
         "I0_ven[2]",
         "I0_ven[3]",
         "I0_ven[4]",
         "I0_it[1]",
         "I0_it[2]",
         "I0_it[3]",
         "I0_it[4]",
         "k",
         "tau")

# Run model with testing of symptomatic only -----------------------------------

test_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/symp_test_prev",
    model_path = "model/est_test_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars

  )
)

test_symp_prev$log()


# Run model with testing of asymptomatic plus all symptomatic isolate ----------

test_asymp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_test_prev",
    model_path = "model/est_test_asymp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    scale_time_step = 2
    
  )
)

test_asymp_prev$log() 

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test_prev",
    model_path = "model/est_test_asymp_and_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars
  )
  )

test_asymp_and_symp_prev$log() 


