

#  Set up to run on HPC  -------------------------------------------------------
rm(list = ls())
root = "Q:/COV_Italy_multistrain/model_fitting"
setwd("Q:/COV_Italy_multistrain/model_fitting")


ctx <- context::context_save(root, packages=c("tidyverse", "rstan", "bayesplot", "Hmisc", 
                                              "BH", "RcppEigen", "png", "knitr", "rstudioapi",
                                              "StanHeaders", "cowplot"),
                             sources=c("R/model_fitting_on_cluster.R"))

# Set up to run with multiple cores on clusters
config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")


# Create a queue within the context (/environment)
obj <- didehpc::queue_didehpc(ctx, config)

# Run model with testing of symptomatic only -----------------------------------
test_symp =  obj$enqueue(
model_fittting_on_cluster(
  file_path= "Results/symp_test",
  model_path = "model/est_test_symp.stan"
)
)

test_symp$log()

# Run model with testing of asymptomatic plus all symptomatic isolate ----------
test_asymp =  obj$enqueue(
  model_fittting_on_cluster(
    file_path= "Results/asymp_test",
    model_path = "model/est_test_asymp.stan"
  )
)



# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
  model_fittting_on_cluster(
    file_path= "Results/asymp_and_symp_test",
    model_path = "model/est_test_asymp_and_symp.stan"
  )
)
