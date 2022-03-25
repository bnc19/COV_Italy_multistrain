#################################################################
# Sample posterior chains 
#################################################################


#################################################################
# Set up 
#################################################################

# rm(list = ls())


number_of_samples = 100

# posterior chains Veneto
posterior_chains = read.csv("Veneto/posteriorChains.csv")[-1]


# posterior chains Italy
posterior_chains_italy = read.csv("Italy/posteriorChains.csv")[-1]

#################################################################
# Sample from posterior chains 
#################################################################
set.seed(23122023)

# sample for Italy and Veneto
posterior_index = sample(size = number_of_samples,
                         1:dim(posterior_chains)[1],
                         replace = T)
set.seed(23122022)

posterior_index_italy = sample(size = number_of_samples,
                               1:dim(posterior_chains_italy)[1],
                               replace = T)

bootdata = list()
for (i in 1:number_of_samples) {
  bootdata[[i]] = posterior_chains[posterior_index[i],]
}


bootdata_italy = list()
for (i in 1:number_of_samples) {
  bootdata_italy[[i]] = posterior_chains_italy[posterior_index_italy[i],]
}


posterior_samples = bind_rows(bootdata)
posterior_samples_italy = bind_rows(bootdata_italy)



#################################################################
# Calculate mean and 95% CrI posterior estimates 
#################################################################

Veneto_post_est  = cbind(
  apply (posterior_samples, 2,  mean),
  apply(posterior_samples, 2, quantile, probs = 0.025),
  apply(posterior_samples, 2, quantile, probs = 0.975)
)
Italy_post_est = cbind(
  apply (posterior_samples_italy, 2,  mean),
  apply(posterior_samples_italy, 2, quantile, probs = 0.025),
  apply(posterior_samples_italy, 2, quantile, probs = 0.975)
)



#################################################################
# Save
#################################################################

write.csv(posterior_samples, "Veneto/100_posterior_samples_veneto.csv")
write.csv(posterior_samples_italy, "Italy/100_posterior_samples_italy.csv")



write.csv(Veneto_post_est, "Veneto/Veneto_post_est.csv")
write.csv(Italy_post_est, "Italy/Italy_post_est.csv")
