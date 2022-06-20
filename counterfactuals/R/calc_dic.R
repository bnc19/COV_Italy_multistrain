# DIC functions

calc_dic = function (
  posterior_chains,
  log_lik,
  fixed_model_path 
){

source("R/sample_posterior_chains.R")


post_chains_sum = summarise_posterior_chains(posterior_chains = posterior_chains)

model_posts = replicate_rstan_fixed(model_path = fixed_model_path,
                                    posterior_sample_row = data.frame(t(post_chains_sum))[1,])


# Calculate model deviance 

E_log_lik =  model_posts$sum_LL # mean post value 

pDIC = 2 * (E_log_lik - mean(rowSums(log_lik)))  # sum across the data points to get total log likelihood then calculate the mean for all iterations 

dic = -2 * E_log_lik  + 2 * pDIC

return(dic)
}
