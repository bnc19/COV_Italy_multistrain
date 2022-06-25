

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
obj$task_get("ef76a48dcbe4660b208316eafff53cea")$log()


# Global parameters ------------------------------------------------------------

n_iter = 1000
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
test_symp =  obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 5,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"),
    adapt_delta = 0.85
  )
)

test_symp$log()  # bound Italy seed, rho and same omega  


# diag ------------------------------   
test_symp_fit=test_symp$result()

symp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =test_symp_fit ,
    pars = pars_inc,
    file_path = paste0(file_path_inc,"/symp_test")
  )
  )

symp_diag$log()



# Run model with testing of asymptomatic plus all symptomatic test and maybe isolate ----------

test_asymp_2_quick =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test_2/same_omega_bound_seed/same_rho"),
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
              "k")
  )
)

test_asymp_2$log() # # bound Italy seed / same omega and rho   
test_asymp_2_quick$log() ## nound italy seed / diff omega / same rho

# diag ------------------------------ 
test_asymp_2_fit=test_asymp_2$result()

asymp_2_diag = 
  obj$enqueue(
    diagnose_stan_fit(
      fit =test_asymp_2_fit ,
      pars = pars_inc,
      file_path = paste0(file_path_inc,"/asymp_test")
    )
  )

asymp_2_diag$log()



# Run model with testing of asymptomatic plus all symptomatic test and isolate ----------

test_asymp =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_asymp.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 5,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"),
    adapt_delta = 0.85
  )
)

test_asymp$log()  # bound Italy seed, rho and same omega  

# diag ------------------------------ 
test_asymp_fit=test_asymp$result()

asymp_diag = 
  obj$enqueue(
    diagnose_stan_fit(
    fit =test_asymp_fit ,
    pars = pars_inc,
    file_path = paste0(file_path_inc,"/asymp_test")
  )
  )

asymp_diag$log()


# Run model with same level of asymp and symp testing --------------------------
test_asymp_and_symp =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_and_symp_test/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = F,
    scale_time_step = 5,
    adapt_delta = 0.85,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k")
  )
)

test_asymp_and_symp$log()   # bound Italy seed, rho and same omega  

# diag ------------------------------ 

asymp_symp_fit=test_asymp_and_symp$result()

asymp_symp_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =asymp_symp_fit ,
    pars = pars_inc,
    file_path = paste0(file_path_inc,"/symp_test")
  ))

asymp_symp_diag$log()



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

test_symp_prev_rho2=  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/symp_test/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    pars= c("lp__", 
                      "beta[1]","beta[2]","beta[3]","beta[4]",
                      "rho_it" , "rho_ven" ,
                      "omega[1]", "omega[2]",
                      "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
                      "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
                      "k", "tau"),
    prev = TRUE,
    scale_time_step = 5,
    adapt_delta = 0.85
  )
)

test_symp_prev$log()  # bound Italy seed and same omega  / finished  
test_symp_prev_rho$log( )  # bound Italy seed, rho and same omega 

# diag ------------------------------  

fit_symp_prev = test_symp_prev$result()

symp_prev_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_symp_prev ,
    pars = pars,
    file_path = paste0(file_path_prev,"/symp_test")
  ))

symp_diag$log()


# Run model with testing of asymptomatic plus all symptomatic test and maybe isolate ----------

test_asymp_prev_2_rho =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_test_2/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_asymp_prev2.stan", 
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = T,
    scale_time_step = 5,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k", "tau"),
    adapt_delta = 0.90
  )
)

test_asymp_prev_2$log()  # bound Italy seed and same omega  # adapt delta was 0.8 but no divergent transitions! 
test_asymp_prev_2_rho$log() # bound Italy seed, rho and same omega 

# diag ------------------------------ 
test_asymp_2_fit=test_asymp_prev_2$result()

asymp_prev_2_diag = 
  obj$enqueue(
    diagnose_stan_fit(
      fit =test_asymp_2_fit ,
      pars = pars_inc,
      file_path = paste0(file_path_inc,"/asymp_test")
    )
  )

asymp_diag$log()



# Run model with testing of asymptomatic plus all symptomatic test and isolate ----------

test_asymp_prev_rho =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_test/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_asymp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    scale_time_step = 5,
    prev = TRUE,
    adapt_delta = 0.90,
    pars= c("lp__", 
      "beta[1]","beta[2]","beta[3]","beta[4]",
      "rho_it" , "rho_ven" ,
      "omega[1]", "omega[2]",
      "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
      "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
      "k", "tau")
    
  )
)


test_asymp_prev$log()  # bound Italy seed and same omega  // running 
test_asymp_prev_rho$log() #  
# diag ------------------------------ 


fit_asymp_prev =  test_asymp_prev$result()

asymp_prev_diag = 
  obj$enqueue(diagnose_stan_fit(
  fit =fit_asymp_prev  ,
  pars = pars,
  file_path = paste0(file_path_prev, "/asymp_test")
))

asymp_prev_diag$log()


# Run model with same level of asymp and symp testing --------------------------

test_asymp_and_symp_prev_rho =  obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_prev,"/asymp_and_symp_test/same_omega_bound_seed/bound_rho"),
    model_path = "model/est_test_asymp_and_symp_prev.stan",
    n_iter = n_iter,
    n_warmups = n_warmups,
    prev = T,
    scale_time_step = 5,
    adapt_delta = 0.90,
    pars= c("lp__", 
            "beta[1]","beta[2]","beta[3]","beta[4]",
            "rho_it" , "rho_ven" ,
            "omega[1]", "omega[2]",
            "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
            "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
            "k", "tau")
  )
  )



test_asymp_and_symp_prev$log()   # bound Italy seed and same omega / running 
test_asymp_and_symp_prev_rho$log() #  

# diag ------------------------------ 


fit_asymp_and_symp_prev  = test_asymp_and_symp_prev$result()

asymp_symp_prev_diag = 
  obj$enqueue(diagnose_stan_fit(
    fit =fit_asymp_and_symp_prev  ,
    pars = pars,
    file_path = paste0(file_path_prev, "/asymp_and_symp_test")
  ))

asymp_symp_prev_diag$log()
