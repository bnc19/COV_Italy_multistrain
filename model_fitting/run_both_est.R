

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

obj$task_get("83a4ee0fa9540d43d77bd871e5516ca6")$log()


# Global parameters ------------------------------------------------------------

n_iter = 1000
n_warmups = round(n_iter/2)
scale_time_step = 2


dir.create(paste0("Results/scale", scale_time_step))
dir.create(paste0("Results/scale", scale_time_step,"/inc"))
dir.create(paste0("Results/scale", scale_time_step,"/prev"))

# 2 = adapt delta 0.8
# otherwise adpat delta  0.99

############################# Fit models to incidence ##########################
file_path_inc = paste0("Results/scale", scale_time_step, "/inc")

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
# 9310d258dbdbd72505c635531d4fd86f

test_symp2 $log()
# 71a2e95e5cf7b6eb8db5f16f9181e49d

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
# 7e63ee4b0dbf8ba616a2950029e78335

test_asymp2$log()

# 84467ed3a1a44dc99c3c83898509396c

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
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
# 0d80c1439437e55211fceac07178d832


test_asymp_and_symp2$log()
# 6bd50e7be79a1b1f256b5827f9b7e329


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

test_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/symp_test"),
    model_path = "model/est_test_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = TRUE,
    scale_time_step = scale_time_step
  )
)

test_symp_prev$log()
#  9d29bd36c1c80cf6cd427b21e746c486


test_symp_prev2$log()
#  70b532c86ad975626d2d1b3454b9f233


# Run model with testing of asymptomatic plus all symptomatic isolate ----------

test_asymp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_test"),
    model_path = "model/est_test_asymp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    scale_time_step = scale_time_step,
    prev = TRUE
    
  )
)

test_asymp_prev$log()
# 157b812d61c524a2e6a6a6f8be8b1207

test_asymp_prev2$log()

# 24be56c409253d91212d35a55f99c788

# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp_prev =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_and_symp_test"),
    model_path = "model/est_test_asymp_and_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars = pars,
    prev = T,
    scale_time_step = scale_time_step
  )
  )

test_asymp_and_symp_prev$log() 

# 6d26c448b7ada5c69e7b4c656fa2557f


test_asymp_and_symp_prev2$log()

#  6dd9783506377d89255ec32c5d295a98