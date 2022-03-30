
################    Sample posterior chains     ################   


################    Set up ################   

library(tidyverse)

# rm(list = ls())


number_of_samples = 100

################    Import data ################   


# posterior chains Veneto
posterior_chains = read.csv("counterfactuals/Veneto/posteriorChains.csv")[-1]


# posterior chains Italy
posterior_chains_italy = read.csv("counterfactuals/Italy/posteriorChains.csv")[-1]




# posterior chains Veneto SA
posterior_chains_sens = read.csv("counterfactuals/Veneto/posteriorChains_sens.csv")[-1]


# posterior chains Italy SA
posterior_chains_italy_sens = read.csv("counterfactuals/Italy/posteriorChains_sens.csv")[-1]

################    Sample from posterior chains ################   
set.seed(23122023)

# sample for Italy and Veneto
posterior_index = sample(size = number_of_samples,
                         1:dim(posterior_chains)[1],
                         replace = T)

set.seed(23122023)

# sample for Italy and Veneto
posterior_index_sens = sample(size = number_of_samples,
                              1:dim(posterior_chains_sens)[1],
                              replace = T)

set.seed(23122022)

posterior_index_italy = sample(size = number_of_samples,
                               1:dim(posterior_chains_italy)[1],
                               replace = T)

set.seed(23122022)

posterior_index_italy_sens = sample(size = number_of_samples,
                               1:dim(posterior_chains_italy_sens)[1],
                               replace = T)



bootdata = bootdata_italy = bootdata_sens =  bootdata_italy_sens= list()


for (i in 1:number_of_samples) {
  bootdata[[i]] = posterior_chains[posterior_index[i],]
  bootdata_sens[[i]] = posterior_chains_sens[posterior_index_sens[i],]
  bootdata_italy[[i]] = posterior_chains_italy[posterior_index_italy[i],]
  bootdata_italy_sens[[i]] = posterior_chains_italy_sens[posterior_index_italy_sens[i],]
  
}

posterior_samples = bind_rows(bootdata)
posterior_samples_italy = bind_rows(bootdata_italy)
posterior_samples_sens = bind_rows(bootdata_sens)
posterior_samples_italy_sens= bind_rows(bootdata_italy_sens)


################    Calculate mean and 95% CrI posterior estimates ################   

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


Veneto_post_est_sens  = cbind(
  apply (posterior_samples_sens, 2,  mean),
  apply(posterior_samples_sens, 2, quantile, probs = 0.025),
  apply(posterior_samples_sens, 2, quantile, probs = 0.975)
)
Italy_post_est_sens = cbind(
  apply (posterior_samples_italy_sens, 2,  mean),
  apply(posterior_samples_italy_sens, 2, quantile, probs = 0.025),
  apply(posterior_samples_italy_sens, 2, quantile, probs = 0.975)
)


################ save ################    

write.csv(posterior_samples, "counterfactuals/Veneto/100_posterior_samples_veneto.csv")
write.csv(posterior_samples_italy, "counterfactuals/Italy/100_posterior_samples_italy.csv")
write.csv(posterior_samples_sens, "counterfactuals/Veneto/100_posterior_samples_veneto_sens.csv")
write.csv(posterior_samples_italy_sens, "counterfactuals/Italy/100_posterior_samples_italy_sens.csv")



write.csv(Veneto_post_est, "counterfactuals/Results/Veneto_post_est.csv")
write.csv(Italy_post_est, "counterfactuals/Results/Italy_post_est.csv")
write.csv(Veneto_post_est_sens, "counterfactuals/Results/Veneto_post_est_sens.csv")
write.csv(Italy_post_est_sens, "counterfactuals/Results/Italy_post_est_sens.csv")

