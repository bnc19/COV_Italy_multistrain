

#  Set up to run on HPC  -------------------------------------------------------
rm(list = ls())
root = "Q:/COV_Italy_multistrain/model_fitting"
setwd("Q:/COV_Italy_multistrain/model_fitting")


ctx <- context::context_save(root, packages=c("tidyverse", "rstan", "bayesplot", "Hmisc", 
                                              "BH", "RcppEigen", "png", "knitr", "rstudioapi",
                                              "StanHeaders", "cowplot", "loo"),
                             sources=c("R/run_model_fitting.R", "R/diagnose_stan_fit.R"))

# Set up to run with multiple cores on clusters
config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")



# Create a queue within the context (/environment)
obj <- didehpc::queue_didehpc(ctx, config)


source("R/run_model_fitting.R")
source("R/diagnose_stan_fit.R")
obj$task_get("8ba38475681fb6d61b3fb9fa560e369f")$log()


# Global parameters ------------------------------------------------------------

n_iter = 2000
n_warmups = round(n_iter/2)
scale_time_step = 5


dir.create(paste0("Results" ,"/inc"))
dir.create(paste0("Results", "/prev"))


############################# Fit models to incidence ##########################
file_path_inc = paste0("Results","/inc")

# Run model with testing of symptomatic only -----------------------------------
test_symp =  obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test"),
    model_path = "model/est_test_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = scale_time_step
  )
)

test_symp$log()

test_symp_fit=test_symp$result()

symp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =test_symp_fit ,
    pars = pars,
    file_path = paste0(file_path_inc,"/symp_test")
  ))

symp_diag$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test"),
    model_path = "model/est_test_asymp.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = scale_time_step
  )
)


test_asymp$log()
test_asymp_fit=test_asymp$result()

asymp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =test_asymp_fit ,
    pars = pars,
    file_path = paste0(file_path_inc,"/symp_test")
  ))

asymp_diag$log()


# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_and_symp_test"),
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = scale_time_step,
    adapt_delta = 0.85
  )
)

test_asymp_and_symp$log()  # 10 div
############################# Fit models to prevalence ##########################


# Add Tau 
pars = c("lp__", 
         "beta[1]","beta[2]","beta[3]","beta[4]",
         "rho_it" , "rho_ven" ,
         "omega[1]", "omega[2]","omega[3]","omega[4]",
         "I0_it[1]", "I0_it[2]","I0_it[3]","I0_it[4]",
         "I0_ven[1]", "I0_ven[2]","I0_ven[3]","I0_ven[4]",
         "k", "tau")

file_path_prev  = paste0("Results", "/prev")



# Run model with testing of symptomatic only -----------------------------------

test_symp_prev=  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/symp_test"),
    model_path = "model/est_test_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = TRUE,
    scale_time_step = scale_time_step,
    adapt_delta = 0.99
  )
)

test_symp_prev$log()

## diag 

fit_symp_prev = obj$task_get("b42507cb0fd2ea7eed154f41f389742e")$result()

symp_prev_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_symp_prev ,
    pars = pars,
    file_path = paste0(file_path_prev,"/symp_test")
  ))

symp_diag$log()


# Run model with testing of asymptomatic plus all symptomatic isolate ----------

test_asymp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_test"),
    model_path = "model/est_test_asymp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    scale_time_step = scale_time_step,
    prev = TRUE,
    adapt_delta = 0.99
    
  )
)

test_asymp_prev$log()

## diag  


fit_asymp_prev =  obj$task_get("cff8091ccbad7869c9bc95a0d5df284a")$result()

asymp_prev_diag = 
  obj$enqueue(diagnose_stan_fit(
  fit =fit_asymp_prev  ,
  pars = pars,
  file_path = paste0(file_path_prev, "/asymp_test")
))

asymp_prev_diag$log()


# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_and_symp_test"),
    model_path = "model/est_test_asymp_and_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = T,
    scale_time_step = scale_time_step,
    adapt_delta = 0.99
  )
  )


test_asymp_and_symp_prev$log()


## diag 


fit_asymp_and_symp_prev  =obj$task_get("3f13a3f65869f8f508e9a7a37a2c89f5")$result()

asymp_symp_prev_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_asymp_and_symp_prev  ,
    pars = pars,
    file_path = paste0(file_path_prev, "/asymp_and_symp_test")
  ))

asymp_symp_prev_diag$log()
