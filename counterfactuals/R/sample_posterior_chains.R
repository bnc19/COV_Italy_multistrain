# Function to sample from posterior chains -------------------------------------

sample_posterior_chains = function(
  posterior_chains,
  number_of_samples,
  seed = 1){
  
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


# Function to calculate mean and 95% CrI from posterior chains -----------------



summarise_posterior_chains = function( posterior_chains)  {

  posterior_chains = posterior_chains %>%  as.data.frame() %>% 
   mutate( rho_v =  rho_v*100,
           rho_i =  rho_i*100, 
           omega_i1=  omega_i1*100,  
           omega_i2=  omega_i2*100)
  
  summary = data.frame(
  mean  = apply(posterior_chains, 2,  mean),
  lower = apply(posterior_chains, 2, quantile, probs = 0.025),
  upper = apply(posterior_chains, 2, quantile, probs = 0.975)
)


summary = round(summary, 1)
 
Out = summary  %>% 
  mutate(parameters = rownames(summary)) %>%  
  mutate(CrI =  paste0("(", lower, "-", upper,")")) %>% 
  mutate(Mean = ifelse(
    grepl("rho|omega", parameters), paste0(mean, "%"), mean)) %>% 
  mutate(out = paste(Mean, CrI)) %>%  
  filter(parameters!= "rho_i")



return(Out[-1,])

}

