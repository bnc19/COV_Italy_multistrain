 
# func calc R0 for main baseline model ----------------------------------------- 
calculate_R0 =   function(
    sigma = 1 / (5.1 - 1.31),
    gamma = 1 / 2.1 ,
    mu = 0.59 ,
    phi_PCR = 0.920,
    posterior_chains_beta,
    posterior_chains_rho
  ){
  
  R0 = matrix(0, nrow = length(posterior_chains_rho$rho_v), ncol = 4)
  R0_sum =  matrix(0, nrow = 4, ncol = 3)
  delta =  as.numeric(posterior_chains_rho$rho_v) * phi_PCR
  
  for(i in 1:4){
    R0[,i] = as.numeric(posterior_chains_beta[,i]) / sigma + ((1-delta) * mu * as.numeric(posterior_chains_beta[,i])) / gamma + ((1-mu) * as.numeric(posterior_chains_beta[,i])) / gamma
    
    R0_sum[i, ] = c(
      mean = mean(R0[,i]),
      lower = quantile(R0[ ,i], probs = 0.025),
      upper = quantile(R0[,i], probs = 0.975)
    )
  
  }
  
  R0_sum = round(R0_sum,2)
  
  out = R0_sum  %>% 
    as.data.frame() %>%  
    mutate(parameters = c("R0_M", "R0_A", "R0_O", "R0_Al")) %>%  
    rename(mean = V1, lower= V2, upper = V3) %>% 
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(out = paste(mean, CrI))   
 
    
     return(out)
}


# Func to calc R0 for baseline model including alpha ---------------------------


calculate_R0_alpha =   function(
  sigma = 1 / (5.1 - 1.31),
  gamma = 1 / 2.1 ,
  mu = 0.59 ,
  phi_PCR = 0.920,
  posterior_chains_beta,
  posterior_chains_rho,
  posterior_chains_alpha 
){
  
  
  R0 = matrix(0, nrow = length(posterior_chains_rho$rho_v), ncol = 4)
  R0_sum =  matrix(0, nrow = 4, ncol = 3)
  delta =  as.numeric(posterior_chains_rho$rho_v) * phi_PCR
  
  alpha = 1 - as.numeric(posterior_chains_alpha$alpha)
  
  for(i in 1:4){
    R0[,i] = as.numeric(posterior_chains_beta[,i]) / sigma + ((1-delta) * mu * alpha * as.numeric(posterior_chains_beta[,i])) / gamma + ((1-mu) * as.numeric(posterior_chains_beta[,i])) / gamma
    
    R0_sum[i, ] = c(
      mean(R0[,i]),
      quantile(R0[ ,i], probs = 0.025),
      quantile(R0[,i], probs = 0.975)
    )
    
    
    
  }
  
  R0_sum = round(R0_sum,2)

  out = R0_sum  %>% 
    as.data.frame() %>%  
    rename(mean = V1, lower= V2, upper = V3) %>%
    mutate(parameters = c("R0_M", "R0_A", "R0_O", "R0_Al")) %>%  
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(out = paste(mean, CrI))

  
  return(out)
}


# Func calculate R0 asymp model ------------------------------------------------

calculate_R0_asymp =   function(
  sigma = 1 / (5.1 - 1.31),
  gamma = 1 / 2.1 ,
  mu = 0.59 ,
  phi_PCR = 0.920,
  posterior_chains_beta,
  posterior_chains_rho
){
  
  R0 = matrix(0, nrow = length(posterior_chains_rho$rho_v), ncol = 4)
  R0_sum =  matrix(0, nrow = 4, ncol = 3)
  delta =  as.numeric(posterior_chains_rho$rho_v) * phi_PCR
  
  for(i in 1:4){
    R0[,i] = as.numeric(posterior_chains_beta[,i]) / sigma + ((1-delta) * (1-mu) * as.numeric(posterior_chains_beta[,i])) / gamma
    
    R0_sum[i, ] = c(
      mean = mean(R0[,i]),
      lower = quantile(R0[ ,i], probs = 0.025),
      upper = quantile(R0[,i], probs = 0.975)
    )
    
  }
  
  R0_sum = round(R0_sum,2)
  
  out = R0_sum  %>% 
    as.data.frame() %>%  
    mutate(parameters = c("R0_M", "R0_A", "R0_O", "R0_Al")) %>%  
    rename(mean = V1, lower= V2, upper = V3) %>% 
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(out = paste(mean, CrI))   
  
  
  return(out)
}



# Func calculate R0 asymp2 model ------------------------------------------------

calculate_R0_asymp2 =   function(
  sigma = 1 / (5.1 - 1.31),
  gamma = 1 / 2.1 ,
  mu = 0.59 ,
  phi_PCR = 0.920,
  posterior_chains_beta,
  posterior_chains_rho
){
  
  R0 = matrix(0, nrow = length(posterior_chains_rho$rho_v), ncol = 4)
  R0_sum =  matrix(0, nrow = 4, ncol = 3)
  delta_symp =   phi_PCR
  delta_asymp =  as.numeric(posterior_chains_rho$rho_v) * phi_PCR
  
  for(i in 1:4){
    R0[,i] = as.numeric(posterior_chains_beta[,i]) / sigma + ((1-delta_asymp) * (1-mu) * as.numeric(posterior_chains_beta[,i])) / gamma +  ((1-delta_symp) * mu * as.numeric(posterior_chains_beta[,i])) / gamma 
    
    R0_sum[i, ] = c(
      mean = mean(R0[,i]),
      lower = quantile(R0[ ,i], probs = 0.025),
      upper = quantile(R0[,i], probs = 0.975)
    )
    
  }
  
  R0_sum = round(R0_sum,2)
  
  out = R0_sum  %>% 
    as.data.frame() %>%  
    mutate(parameters = c("R0_M", "R0_A", "R0_O", "R0_Al")) %>%  
    rename(mean = V1, lower= V2, upper = V3) %>% 
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(out = paste(mean, CrI))   
  
  
  return(out)
}



# func calc R0 for asymp / symp model ------------------------------------------ 
calculate_R0_asymp_symp =   function(
  sigma = 1 / (5.1 - 1.31),
  gamma = 1 / 2.1 ,
  mu = 0.59 ,
  phi_PCR = 0.920,
  posterior_chains_beta,
  posterior_chains_rho
){
  
  R0 = matrix(0, nrow = length(posterior_chains_rho$rho_v), ncol = 4)
  R0_sum =  matrix(0, nrow = 4, ncol = 3)
  delta =  as.numeric(posterior_chains_rho$rho_v) * phi_PCR
  
  for(i in 1:4){
    R0[,i] = as.numeric(posterior_chains_beta[,i]) / sigma + ((1-delta) * as.numeric(posterior_chains_beta[,i])) / gamma 
    
    R0_sum[i, ] = c(
      mean = mean(R0[,i]),
      lower = quantile(R0[ ,i], probs = 0.025),
      upper = quantile(R0[,i], probs = 0.975)
    )
    
  }
  
  R0_sum = round(R0_sum,2)
  
  out = R0_sum  %>% 
    as.data.frame() %>%  
    mutate(parameters = c("R0_M", "R0_A", "R0_O", "R0_Al")) %>%  
    rename(mean = V1, lower= V2, upper = V3) %>% 
    mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
    mutate(out = paste(mean, CrI))   
  
  
  return(out)
}

