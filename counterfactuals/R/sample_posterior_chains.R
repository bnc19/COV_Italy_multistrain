# Function to sample from posterior chains -------------------------------------

sample_posterior_chains = function(
  posterior_chains,
  number_of_samples,
  seed = 1){
  
  set.seed(seed)
  posterior_index = sample(size = number_of_samples,
                           1:dim(posterior_chains)[1],
                           replace = T)
  
  bootdata = list()
  
  for (i in 1:number_of_samples) {
    bootdata[[i]] = posterior_chains[posterior_index[i],]
  }
  
  posterior_samples = bind_rows(bootdata)
  
  return(posterior_samples)
  
}



# Function to calculate mean  from posterior chains for DIC --------------------



summarise_posterior_chains_DIC = function( posterior_chains)  {
  
  
  summary = data.frame(
    mean  = apply(posterior_chains, 2,  mean),
    lower = apply(posterior_chains, 2, quantile, probs = 0.025),
    upper = apply(posterior_chains, 2, quantile, probs = 0.975)
  )
  return(summary[-1,])
  
}

# Function to calculate mean and 95% CrI ---------------------------------------
# from posterior chains for baseline model manuscript --------------------------



summarise_posterior_chains_symp = function( posterior_chains, alpha = FALSE, 
                                       sigma = 1 / (5.1 - 1.31),
                                       gamma = 1 / 2.1 ,
                                       mu = 0.59 ,
                                       phi_PCR = 0.920)  {
  
  source("R/calculate_R0.R")
  
  
  posterior_chains_beta = posterior_chains %>%  
    select(beta_M,beta_A,beta_O,beta_Al)
  posterior_chains_rho = posterior_chains %>%  
    select(rho_v)
  
  
  if(alpha==T){
    
    
    posterior_chains_alpha = posterior_chains %>%  
      select(alpha)
    
  posterior_chains = posterior_chains %>%  as.data.frame() %>% 
    mutate( rho_v =  rho_v*100,
            rho_i =  rho_i*100, 
            omega_i1=  omega_i1*100,  
            omega_i2=  omega_i2*100,
            alpha = (1 - alpha)* 100 )
  

  
  R0 = calculate_R0_alpha (
    sigma = sigma,
    gamma = gamma ,
    mu = mu ,
    phi_PCR = phi_PCR,
    posterior_chains_beta = posterior_chains_beta,
    posterior_chains_rho = posterior_chains_rho,
    posterior_chains_alpha = posterior_chains_alpha
  )
  
  } else {
    posterior_chains = posterior_chains %>%  as.data.frame() %>% 
      mutate( rho_v =  rho_v*100,
              rho_i =  rho_i*100, 
              omega_i1=  omega_i1*100,  
              omega_i2=  omega_i2*100)
    
    R0 = calculate_R0 (
      sigma = sigma,
      gamma = gamma ,
      mu = mu ,
      phi_PCR = phi_PCR,
      posterior_chains_beta = posterior_chains_beta,
      posterior_chains_rho = posterior_chains_rho)
    
  }
  
  summary = data.frame(
    mean  = apply(posterior_chains, 2,  mean),
    lower = apply(posterior_chains, 2, quantile, probs = 0.025),
    upper = apply(posterior_chains, 2, quantile, probs = 0.975)
  )
  
  
  summary = round(summary, 2)
  
  sum = summary  %>% 
    mutate(parameters = rownames(summary)) %>%  
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(Mean = ifelse(
      grepl("rho|omega|alpha", parameters), paste0(mean, "%"), mean), 
      ) %>% 
    mutate(out = paste(Mean, CrI)) %>%  
    filter(parameters!= "rho_i") %>%  
    select(!Mean)
  
  
 
  
  out=rbind(R0,sum[-1,])
  
  
  return(out)
  
}




# Function to calculate mean and 95% CrI ---------------------------------------
# from posterior chains for baseline model manuscript --------------------------



summarise_posterior_chains_not_symp = function( posterior_chains, alpha = FALSE, 
                                            sigma = 1 / (5.1 - 1.31),
                                            gamma = 1 / 2.1 ,
                                            mu = 0.59 ,
                                            phi_PCR = 0.920,
                                            asymp=F,
                                            asymp2=F,
                                            asymp_symp=F)  {
  
  source("R/calculate_R0.R")
  
  
  posterior_chains_beta = posterior_chains %>%  
    select(beta_M,beta_A,beta_O,beta_Al)
  posterior_chains_rho = posterior_chains %>%  
    select(rho_v)
  
  
  posterior_chains = posterior_chains %>%  as.data.frame() %>% 
    mutate( rho_v =  rho_v*100,
            rho_i =  rho_i*100, 
            omega_i1=  omega_i1*100,  
            omega_i2=  omega_i2*100)
  
  if(asymp==T){
    R0 = calculate_R0_asymp (
      sigma = sigma,
      gamma = gamma ,
      mu = mu ,
      phi_PCR = phi_PCR,
      posterior_chains_beta = posterior_chains_beta,
      posterior_chains_rho = posterior_chains_rho
    )
    
  } else if (asymp2 == T){
  
    R0 = calculate_R0_asymp2 (
      sigma = sigma,
      gamma = gamma ,
      mu = mu ,
      phi_PCR = phi_PCR,
      posterior_chains_beta = posterior_chains_beta,
      posterior_chains_rho = posterior_chains_rho)
    
  } else if ( asymp_symp == T){
    R0 = calculate_R0_asymp_symp (
      sigma = sigma,
      gamma = gamma ,
      mu = mu ,
      phi_PCR = phi_PCR,
      posterior_chains_beta = posterior_chains_beta,
      posterior_chains_rho = posterior_chains_rho)
    } 


  
  summary = data.frame(
    mean  = apply(posterior_chains, 2,  mean),
    lower = apply(posterior_chains, 2, quantile, probs = 0.025),
    upper = apply(posterior_chains, 2, quantile, probs = 0.975)
  )
  
  
  summary = round(summary, 2)
  
  sum = summary  %>% 
    mutate(parameters = rownames(summary)) %>%  
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(Mean = ifelse(
      grepl("rho|omega|alpha", parameters), paste0(mean, "%"), mean), 
    ) %>% 
    mutate(out = paste(Mean, CrI)) %>%  
    filter(parameters!= "rho_i") %>%  
    select(!Mean)
  
  
  
  
  out=rbind(R0,sum[-1,])
  
  
  return(out)
  
}

