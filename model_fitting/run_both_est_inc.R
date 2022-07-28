

#  Set up to run on HPC  -------------------------------------------------------
rm(list = ls())
root = "Q:/COV_Italy_multistrain_2/model_fitting"
setwd("Q:/COV_Italy_multistrain_2/model_fitting")


ctx <- context::context_save(root, packages=c("tidyverse", "rstan", "bayesplot", "Hmisc", 
                                              "BH", "RcppEigen", "png", "knitr", "rstudioapi",
                                              "StanHeaders", "cowplot", "loo"),
                             sources=c("R/run_model_fitting.R", "R/diagnose_stan_fit.R"))

# Set up to run with multiple cores on clusters
config <- didehpc::didehpc_config(cores = 8,  parallel =T, cluster = "dideclusthn")


# Create a queue within the context (/environment)
obj <- didehpc::queue_didehpc(ctx, config)

source("R/summarise_posterior_chains.R")
source("R/run_model_fitting.R")
source("R/diagnose_stan_fit.R")
obj$task_get("092ab9c282d1e194ede113de72e7c344")$log()


# Global parameters ------------------------------------------------------------

n_iter = 2000
n_warmups = round(n_iter/2)



dir.create(paste0("Results" ,"/inc"))
dir.create(paste0("Results", "/prev"))


############################# Fit models to incidence ##########################
file_path_inc = paste0("Results","/inc")

pars_inc = c("lp__", 
             "beta[1]","beta[2]","beta[3]","beta[4]",
             "rho_it" , "rho_ven" ,
             "omega[1]", "omega[2]","omega[3]","omega[4]",
             "I0_it[1]", "I0_it[2]","I0_it[3]","I0_it[4]",
             "I0_ven[1]", "I0_ven[2]","I0_ven[3]","I0_ven[4]",
             "k")

# Run model with testing of symptomatic only -----------------------------------
test_symp_m2 =  obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/2"),
    model_path = "model/est_test_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"),
    adapt_delta = 0.88,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020"
  )
)

test_symp$log()  # bound Italy seed, rho and same omega  
test_symp_m2$log()

summarise_posterior_chains(read.csv(paste0(file_path_inc,"/symp_test/2/posterior_chains.csv")))

# Run model with testing of asymptomatic plus all symptomatic test and isolate ----------

test_asymp_m2 =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test/2"),
    model_path = "model/est_test_asymp.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"),
    adapt_delta = 0.80,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020"
  )
)

test_asymp$log()  # bound Italy seed, rho and same omega  
test_asymp_m2$log() 


summarise_posterior_chains(read.csv(paste0(file_path_inc,"/asymp_test/2/posterior_chains.csv")))

# Run model with testing of asymptomatic plus all symptomatic test and maybe isolate ----------

test_asymp_m2 =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test_2/2"),
    model_path = "model/est_test_asymp2.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"),
    adapt_delta =0.85,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020"
  ),

)

test_asymp_m2$log() 

summarise_posterior_chains(read.csv(paste0(file_path_inc,"/asymp_test_2/1/posterior_chains.csv")))




# Run model with same level of Asymp and Symp testing --------------------------
test_asymp_and_symp_m2 =  obj$enqueue(
run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_and_symp_test/2"),
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 2,
    adapt_delta = 0.85,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020"
  )
)

test_asymp_and_symp$log()   # bound Italy seed, rho and same omega
test_asymp_and_symp_m2$log()
