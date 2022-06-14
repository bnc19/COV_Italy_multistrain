

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

n_iter = 2000
n_warmups = round(n_iter/2)


############################# Fit models to incidence ##########################


# Run model with testing of symptomatic only -----------------------------------
test_symp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/symp_test/scale2",
    model_path = "model/est_test_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2
  )
)

test_symp$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_test/scale2",
    model_path = "model/est_test_asymp.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2
  )
)


test_asymp$log()


# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test/scale2",
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2
  )
)


test_asymp_and_symp$log()



############################# Fit models to prevalence ##########################

# Add Tau 
pars = c("lp__",
         "beta",
         "rho_it" ,
         "rho_ven" ,
         "omega",
         "I0_ven",
         "I0_it",
         "k",
         "tau")

# Run model with testing of symptomatic only -----------------------------------

test_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/symp_test_prev/scale2",
    model_path = "model/est_test_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = TRUE,
    scale_time_step = 2
  )
)

test_symp_prev$log()

# 2bccc9fa3e255a10bfa406b1a1b3aae1 99 adapt delta, scale * 1
# 3e655b735f38af627fd4ef143511bdcb 99 adapt delta, scale * 2 


# Run model with testing of asymptomatic plus all symptomatic isolate ----------

test_asymp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_test_prev/scale2",
    model_path = "model/est_test_asymp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    scale_time_step = 5,
    prev = T
    
  )
)
test_asymp_prev$log()

# cec584a021a2a12ac51d2a02c873240a 99 adapt delta, scale * 2 
# 4d10d61ec7aa9681322cd88e8fd59de9 99 adapt delta, scale * 5 

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= "Results/asymp_and_symp_test_prev/scale2",
    model_path = "model/est_test_asymp_and_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = T,
    scale_time_step = 2
  )
  )

test_asymp_and_symp_prev$log() 


# 4b2975fe6ebd07a3d5767d6471790236 99 adapt delta, scale * 1 
# 00ebe31dfc72258723f4028fdc0aa041 99 adapt delta, scale * 2 

