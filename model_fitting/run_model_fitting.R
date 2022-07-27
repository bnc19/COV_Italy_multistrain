################################         READ ME         ################################   
# 
# This script fits a multivariant model to epidemiological and genomic data from Veneto 
# or the rest of Italy. The model can be run on a high performance cluster (HPC) or
# locally. Currently, each model runs 3 chains for  200 iterations, discarding the first 
# 100 as burnin. Initial values are close to the posterior means. This is to demo the
# models ~ quickly (12 hours on a normal computer). In order to obtain the results 
# provided in the manuscript, 3 chains with more diffuse starting points were run for
# 2000 iterations, discarding the first 1000 iterations as burnin. This takes ~6 
# hours per model to run on the HPC. 
#
################################                         ################################   

# create folders 

dir.create("model_fitting/Results")
dir.create("model_fitting/Results/inc")



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


file_path_inc = paste0("Results","/inc")

# RUN MODEL VARINATS on HPC ----------------------------------------------------

# Symp model 1 -----------------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/1"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
    
  )
)

obj$task_get("f769b631c80f6e79b83ea58c8b9245ab")$log()
fit = obj$task_get("f769b631c80f6e79b83ea58c8b9245ab")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}

# Asymp model 1 ----------------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test/1"),
    model_path = "model/est_test_asymp.stan", 
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )
)


obj$task_get("dac7484fdeb89976f37695fe04b24e13")$log()
fit = obj$task_get("dac7484fdeb89976f37695fe04b24e13")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}


# Asymp2 model 1 ---------------------------------------------------------------

 obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test_2/1"),
    model_path = "model/est_test_asymp2.stan", 
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )
 )


obj$task_get("62ca5578233449270df530f7ab3f7e34")$log()
fit = obj$task_get("62ca5578233449270df530f7ab3f7e34")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}



# Asymp / symp model 1----------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_and_symp_test/1"),
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )
)


obj$task_get("ecba5548a79a8660f9e128f349017510")$log()
fit = obj$task_get("")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}





# RUN MODEL VARIANTS LOCALLY ---------------------------------------------------

# Symp model 1 -----------------------------------------------------------------

symp =  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/1"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
    
  )


# Asymp model 1 ----------------------------------------------------------------

asymp =  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test/1"),
    model_path = "model/est_test_asymp.stan", 
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )

# Asymp2 model 1 ---------------------------------------------------------------

asymp2 =  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_test_2/1"),
    model_path = "model/est_test_asymp2.stan", 
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )

# Asymp / symp model 1----------------------------------------------------------


asymp_symp =  run_model_fitting(
    file_path= paste0(file_path_inc,"/asymp_and_symp_test/1"),
    model_path = "model/est_test_asymp_and_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )





# RUN SENSITIVITY ANALYSES  ON SYMP MODEL ON HPC -------------------------------

# M2 - (0.643) -----------------------------------------------------------------
obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/2"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.643
  )
)

obj$task_get("55f33eca06b6181ca8940e2791496c51")$log()
fit = obj$task_get("55f33eca06b6181ca8940e2791496c51")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}

# M3 - (0.875) ----------------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/3"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.875
  )
)

obj$task_get("ba4b73cf60a0bbc83272311195402ebc")$log()

fit = obj$task_get("")$result()
if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}


# M4 - (5.6 inc period) -------------------------------------------------------
obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/4"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99 ,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689 ,
    sigma = 1 / (5.6 - 1.31),
    gamma = 1 / (7.2 - 5.6)
  )
)


obj$task_get("f1aa00140d983be2263a54eb4508dc35")$log()
fit = obj$task_get("")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}



# M5 - (2.15 latent period) -------------------------------------------------------
obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/5"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99 ,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689 ,
    epsilon = 1/ 2.15, 
    sigma = 1 / (5.1 - 2.15),
    gamma = 1 / (7.2 - 5.1)
  )
)


obj$task_get("4ec136d0948f3f53b7ceb2d1e1c3f5a6")$log()
fit = obj$task_get("")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}




# M6 - (0.689, VE = 0.7) ----------------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/6"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689,
    VE = 0.7
    
  )
)


obj$task_get("8b0b24a6364cf9cf32aea664d44c984a")$log()

fit = obj$task_get("")$result()


if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}



# M7 - (0.689, VE = 0.8) ----------------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/7"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689,
    VE = 0.8
    
  )
)


obj$task_get("c1dde9e272627d2779d0b0684b279979")$log()
fit = obj$task_get("")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}


# M8 - (alpha) ----------------------------------------------------------------

obj$enqueue(
  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/8"),
    model_path = "model/est_test_symp_alpha.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k", "alpha"
    ),
    adapt_delta = 0.95,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )
)


obj$task_get("734157ea8bfba34bafee3db3f8600632")$log()

fit = obj$task_get("")$result()

if(any(rhat(fit) > 1.05, na.rm = T)){
  print("Rhat too high")
} else{
  print("Rhat good")}





# RUN SENSITIVITY ANALYSES  ON SYMP MODEL LOCALLY -------------------------------

# M2 - (0.643) -----------------------------------------------------------------

 M2= run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/2"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.643
  )


# M3 - (0.875) ----------------------------------------------------------------


M3=  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/3"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.875
  )


# M4 - (5.6 inc period) -------------------------------------------------------

M4=  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/4"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99 ,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689 ,
    sigma = 1 / (5.6 - 1.31),
    gamma = 1 / (7.2 - 5.6)
  )

# M5 - (2.15 latent period) -------------------------------------------------------

M5=  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/5"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99 ,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689 ,
    epsilon = 1/ 2.15, 
    sigma = 1 / (5.1 - 2.15),
    gamma = 1 / (7.2 - 5.1)
  )



# M6 - (0.689, VE = 0.7) ----------------------------------------------------------------


M6=  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/6"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689,
    VE = 0.7
    
  )


# M7 - (0.689, VE = 0.8) ----------------------------------------------------------------


M7=  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/7"),
    model_path = "model/est_test_symp.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k"
    ),
    adapt_delta = 0.99,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689,
    VE = 0.8
    
  )


# M8 - (alpha) ----------------------------------------------------------------

M8 =  run_model_fitting(
    file_path=  paste0(file_path_inc,"/symp_test/8"),
    model_path = "model/est_test_symp_alpha.stan",
    n_iter =2000,
    n_warmups = 1000,
    prev = F,
    scale_time_step = 2,
    pars =  c("lp__", 
              "beta[1]","beta[2]","beta[3]","beta[4]",
              "rho_it" , "rho_ven" ,
              "omega[1]", "omega[2]", # "omega[3]","omega[4]",
              "I0_it_M", "I0_it_A","I0_it_O","I0_it_Al",
              "I0_ven_M", "I0_ven_A","I0_ven_O","I0_ven_Al",
              "k", "alpha"
    ),
    adapt_delta = 0.95,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020",
    phi_Ag =  0.689
  )
