# DIC functions

# Function to calculate the log-lik

# input:
# - log_lik: an N x I matrix of the pointwise log likelihood extracted during model fitting, 
#            N = no. data points, I = no. iterations 


calculate_log_like = function(
  log_lik
){
  
  # calculate total model log likelihood 
  
  sum_log_lik = rowSums(log_lik) # total LL by iteration 
  
   out  = data.frame(
   mean = mean(sum_log_lik), # mean LL 
   lower = quantile(sum_log_lik, probs = 0.025),
   upper = quantile(sum_log_lik, probs = 0.975))

  
  
  return(round(out,1))
}
