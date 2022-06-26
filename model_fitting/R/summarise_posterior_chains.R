

summarise_posterior_chains = function( posterior_chains)  {
  summary = data.frame(
    mean  = apply(posterior_chains, 2,  mean),
    lower = apply(posterior_chains, 2, quantile, probs = 0.025),
    upper = apply(posterior_chains, 2, quantile, probs = 0.975)
  )
  
  return(summary[-1 ,])
  
}

