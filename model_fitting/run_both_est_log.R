

#  Set up to run on HPC  -------------------------------------------------------
rm(list = ls())
root = "Q:/COV_Italy_multistrain2/model_fitting"
setwd("Q:/COV_Italy_multistrain2/model_fitting")


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
obj$task_get("dc7c8d55baf40bc205d9b9a3849daa2e")$log()


# Global parameters ------------------------------------------------------------

n_iter = 1000
n_warmups = round(n_iter/2)
scale_time_step = 6


dir.create(paste0("Results/scale", scale_time_step))
dir.create(paste0("Results/scale", scale_time_step,"/inc"))
dir.create(paste0("Results/scale", scale_time_step,"/prev"))


############################# Fit models to incidence ##########################
file_path_inc = paste0("Results/scale", scale_time_step, "/inc")

# Run model with testing of symptomatic only -----------------------------------
test_symp2 =  obj$enqueue(
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
# f97022895f5d62f9139c71ac39ae279a

test_symp2 $log()
# b7f60fceb3d5d394ba52c03e1214eb41

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp2 =  obj$enqueue(
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
# 49ea207c56135704f5df31f8c4a89263

test_asymp2$log()

# d57a018949494a880ae5e0cec30421a0

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp2 =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_and_symp_test"),
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = scale_time_step
  )
)


test_asymp_and_symp$log()
# d864226d68bb0504b3a1565a2ce71387


test_asymp_and_symp2$log()
# a089e577aa8e7e6730a7af89ce08c835


############################# Fit models to prevalence ##########################


# Add Tau 
pars = c("lp__", 
         "beta[1]","beta[2]","beta[3]","beta[4]",
         "rho_it" , "rho_ven" ,
         "omega[1]", "omega[2]","omega[3]","omega[4]",
         "I0_it[1]", "I0_it[2]","I0_it[3]","I0_it[4]",
         "I0_ven[1]", "I0_ven[2]","I0_ven[3]","I0_ven[4]",
         "k", "tau")

file_path_prev  = paste0("Results/scale", scale_time_step, "/prev")



# Run model with testing of symptomatic only -----------------------------------

test_symp_prev=  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/symp_test"),
    model_path = "model/log/est_test_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = TRUE,
    scale_time_step = scale_time_step,
    adapt_delta = 0.99
  )
)

test_symp_prev$log()

# 8ba38475681fb6d61b3fb9fa560e369f bad priors?

## diag 

fit_symp = obj$task_get("8ba38475681fb6d61b3fb9fa560e369f")$result()

symp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_symp ,
    pars = pars,
    file_path = paste0(file_path_prev,"/symp_test")
  ))

symp_diag$log()


# Run model with testing of asymptomatic plus all symptomatic isolate ----------

test_asymp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_test"),
    model_path = "model/log/est_test_asymp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    scale_time_step = scale_time_step,
    prev = TRUE,
    adapt_delta = 0.99
    
  )
)


test_asymp_prev$log()

# 3fb8138800fd7d9c6493496f047ad126 bad priors?

## diag  


fit_asymp =  obj$task_get("fba3228c9bfd82a3787b82627a74a878")$result()

asymp_diag = 
  obj$enqueue(diagnose_stan_fit(
  fit =fit_asymp ,
  pars = pars,
  file_path = paste0(file_path_prev, "LN_tau/asymp_test")
))

asymp_diag$log()


# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_and_symp_test"),
    model_path = "model/log/est_test_asymp_and_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = T,
    scale_time_step = scale_time_step,
    adapt_delta = 0.99
  )
  )


test_asymp_and_symp_prev$log()

# 3ad50ef56b49609752266a40672cc55c bad priors?

obj$task_get("341c2494c602d68b14ae40e4eb2d42f3")$log()
## diag 
rstan::check_divergences(fit_asymp_and_symp)

fit_asymp_and_symp =test_asymp_and_symp_prev$result()

asymp_symp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_asymp_and_symp ,
    pars = pars,
    file_path = paste0(file_path_prev, "/asymp_and_symp_test")
  ))

asymp_symp_diag$log()
