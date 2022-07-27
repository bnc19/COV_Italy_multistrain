

save_post_chains = function(fit,file_path){
  
  
  # extract data from posterior
  fit_posts =  rstan::extract(fit)
  
  
  beta_M =  fit_posts$beta[,1]
  beta_A =  fit_posts$beta[,2]
  beta_O =  fit_posts$beta[,3]
  beta_Al =  fit_posts$beta[,4]
  
  rho_v =  fit_posts$rho_ven
  rho_i =  fit_posts$rho_it

  omega_i1 =  fit_posts$omega[,1]
  omega_i2 =  fit_posts$omega[,2]
  
  # omega_v1 =  fit_posts$omega[,3]
  # omega_v2 =  fit_posts$omega[,4]
  
  I0_i_M =  fit_posts$I0_it_M
  I0_i_A =  fit_posts$I0_it_A
  I0_i_O =  fit_posts$I0_it_O
  I0_i_Al =  fit_posts$I0_it_Al
  
  
  I0_v_M =  fit_posts$I0_ven_M
  I0_v_A =  fit_posts$I0_ven_A
  I0_v_O =  fit_posts$I0_ven_O
  I0_v_Al =  fit_posts$I0_ven_Al
  
  
  k = fit_posts$k
  
  alpha = fit_posts$alpha
  
  if(is.null(alpha)){
  posterior_chains = data.frame(
    beta_M,
    beta_A,
    beta_O,
    beta_Al,
    rho_v,
    rho_i,
    # omega_v1,
    # omega_v2,
    omega_i1,
    omega_i2,
    I0_i_M,
    I0_i_A,
    I0_i_O,
    I0_i_Al,
    I0_v_M,
    I0_v_A,
    I0_v_O,
    I0_v_Al,
    k)} else{
      posterior_chains = data.frame(
        beta_M,
        beta_A,
        beta_O,
        beta_Al,
        rho_v,
        rho_i,
        # omega_v1,
        # omega_v2,
        omega_i1,
        omega_i2,
        I0_i_M,
        I0_i_A,
        I0_i_O,
        I0_i_Al,
        I0_v_M,
        I0_v_A,
        I0_v_O,
        I0_v_Al,
        k, 
        alpha)
      
    }
  
  write.csv(posterior_chains, 
            file =  paste0(file_path,"/posterior_chains.csv"))
  
}
