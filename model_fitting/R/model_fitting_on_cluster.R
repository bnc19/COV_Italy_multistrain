


model_fittting_on_cluster = function(file_path,
                                   model_path){
  
# Create file to save folder 
  
  folder = file_path
  
  if (file.exists(folder)) {
    
    cat("The folder already exists")
    
  } else {
    
    dir.create(folder)
    
  }

# Load packages 
  library(bayesplot)
  library(tidyverse)
  library(Hmisc)
  library(cowplot)
  
# Source function 
  source("R/run_stan_model.R")
  source("R/diagnose_stan_fit.R")
  source("R/plot_model_fit.R")
  
  
# Run model 
  fit = run_stan_model(modelPath=model_path)
  
# Diagnostics 
  diagnostics = diagnose_stan_fit(fit=fit,file_path = file_path)
  
# Plot 
  plot = plot_stan_fit(fit=fit,file_path =file_path)
  
  return(list(fit,diagnostics,plot))
}
