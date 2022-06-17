

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
obj$task_get("630462737707f802dc22c6da7607982e")$log()


# Global parameters ------------------------------------------------------------


dir.create(paste0("Results/scale", scale_time_step))
dir.create(paste0("Results/scale", scale_time_step,"/inc"))
dir.create(paste0("Results/scale", scale_time_step,"/prev"))

# 2 = adapt delta 0.8
# otherwise adpat delta  0.99

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

n_iter = 2000
n_warmups = round(n_iter/2)
scale_time_step = 5

# Add Tau 
pars = c("lp__", 
         "beta[1]","beta[2]","beta[3]","beta[4]",
         "rho_it" , "rho_ven" ,
         "omega[1]", "omega[2]","omega[3]","omega[4]",
         "I0_it[1]", "I0_it[2]","I0_it[3]","I0_it[4]",
         "I0_ven[1]", "I0_ven[2]","I0_ven[3]","I0_ven[4]",
         "k", "tau")

dir.create(paste0("Results/scale", scale_time_step,"/prev"))
dir.create(paste0("Results/scale", scale_time_step,"/prev/inform_priors"))
file_path_prev  = paste0("Results/scale", scale_time_step, "/prev")



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

# a49c08fd033a69208d5f930bcf27d5ee og priors, adapt delta .99  scale 5 2000 it
# 7b3107528afcc711cd899f884931f3ed improve priors adapt delta .99  scale 5 2000 it
# 8cc06b0b9a66c040f0d3f92d68205d74 improve priors adapt delta .8 scale 3
# 812846e254be4ea78a57d6e4fe02f96e improve priors adapt delta .99  scale 3

fit_symp = obj$task_get("2255f358c398afa8a77302b81d228041")$result()

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

# 2ef09d4a4ca4364aec7c2f1e23b649c2 og priors, adapt delta .99  scale 5 2000 it
# b61bdbfda4ce5f346edb1d6647062af2 improve priors adapt delta .99  scale 5 2000 it
# b93f54b8d03b6207f21cf120d7b72dbb improve priors adapt delta .8 scale 3   - done 
# 6408bfe35e7c3fbcfc9e0154ac3187a9 improve priors adapt delta .99 scale 3


fit_asymp =  obj$task_get("a6a01d1e24adf5bf859384e583806419")$result()

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


obj$task_get("ab2c12ce996b41ce60917905e8436150")$log()

# ef0a7aa38b0b920f972f8a24464c3853 og priors, adapt delta .99  scale 5 2000 it
# be87b990ccc7b41b3dc5854598fbdba9 improve priors adapt delta .99  scale 5 2000 it
# 37f6aa0a15b9aa5798aa5ccdf157f005  improve priors adapt delta .8 scale 3
# 17b45cdd783aaa1879d2dc00f1cf923d  improve priors adapt delta .99 scale 3

rstan::check_divergences(fit_asymp_and_symp)

fit_asymp_and_symp = obj$task_get("630462737707f802dc22c6da7607982e")$result()

asymp_symp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_asymp_and_symp ,
    pars = pars,
    file_path = paste0(file_path_prev, "/LN_tau/asymp_and_symp_test")
  ))

asymp_symp_diag$log()
