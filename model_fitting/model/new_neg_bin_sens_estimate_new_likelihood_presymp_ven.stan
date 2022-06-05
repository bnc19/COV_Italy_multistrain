data {
  
  int<lower=1> n_var; 
  int<lower = 1> n_data_ven[n_var];      // number of months observed 
  
  int index_M_ven[n_data_ven[1]];
  int index_A_ven[n_data_ven[2]];
  int index_O_ven[n_data_ven[3]];
  int index_Al_ven[n_data_ven[4]];
  
  int<lower = 1> n_months_ven; 
  int<lower = 1> n_days_ven;    // number of months fitted, accounting for initial seed
  int<lower = 1> n_pop_ven;       // population
  int<lower = 1> n_recov_ven;     // seropositive 

  int<lower = 1> n_ts_ven;

  int y_M_ven[n_data_ven[1]];            // positive M  
  int y_A_ven[n_data_ven[2]];            // positive A  
  int y_O_ven[n_data_ven[3]];            // positive O
  int y_Al_ven[n_data_ven[4]];           // positive Alpha 

  real gamma;
  real sigma;
  real mu;
  real epsilon;
  
  real phi_PCR;
  real phi_Ag;

  int PCR_daily_ven[n_ts_ven];
  int Ag_daily_ven[n_ts_ven];
  real vac_ven[n_ts_ven];
  
  int month_index_ven[n_months_ven+1]; // vectors to index months from days and of the number of days in each month 

  int time_switch1_ven;
  int time_switch2_ven;
  
  int time_seed_M_ven;
  int time_seed_alpha_ven;


}

transformed data {
  int S0_ven = n_pop_ven-n_recov_ven;
}

parameters {
  real<lower = 0>             beta[n_var];  // transmission parameter 
  real<lower = 0>             I0[n_var];    // seed
  real<lower = 0, upper = 1>  rho;      // probability of reporting 
  real<lower = 0, upper = 1>  omega[2]; // reduction in transmission
  real<lower = 0>             k;        // overdispersion parameter 

}   

transformed parameters{
  real S_ven[(n_ts_ven +1)];
  real E_ven[(n_ts_ven +1), n_var];
  real PS_ven[(n_ts_ven+1), n_var];
  real IS_ven[(n_ts_ven+1), n_var];
  real IA_ven[(n_ts_ven+1), n_var];
  real Q_ven[(n_ts_ven +1)];
  real R_ven[(n_ts_ven +1)];
  
  real FOI_ven[n_ts_ven, n_var]; 
  real incidence_ven[n_ts_ven,n_var];
  
  // 
  real pPCR_daily_ven[n_ts_ven];     // daily probability of PCR test
  real p_daily_ven[n_ts_ven];        // daily concordant incidence 
  real delta_ven[n_ts_ven, n_var];


  real beta2_ven[n_ts_ven,n_var]; 
  
// Initial values 

  S_ven[1] = S0_ven - I0[2] - I0[3] ;

 for(i in 1:n_var) E_ven[1,i] = 0 ;
 for(i in 1:n_var) PS_ven[1,i] = 0 ;
 for(i in 1:n_var) IA_ven[1,i] = 0 ;
  
  IS_ven[1,1] = 0 ;
  IS_ven[1,2] = I0[2] ;
  IS_ven[1,3] = I0[3] ;
  IS_ven[1,4] = 0 ;
  
  Q_ven[1] = 0 ;
  R_ven[1] = 0 ; 
  
  // model interventions 
  
    
 for(t in 1:n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] =  beta[i];
 } 
 
 for(t in time_switch1_ven :  (time_switch2_ven-1)) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[1]) *beta[i];
 }  
 
 for(t in  time_switch2_ven :  n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[2]) *beta[i];
 }  
 

//--- Simulate --- // 
 for (t in 1:(n_ts_ven)){
   
   IS_ven[time_seed_M_ven+1,1] = I0[1];
   IS_ven[time_seed_alpha_ven+1,4] = I0[4];

   p_daily_ven[t] = (sigma * sum(E_ven[t,2:4])) / S0_ven;
   pPCR_daily_ven[t] = (PCR_daily_ven[t] - Ag_daily_ven[t] * p_daily_ven[t]) / (PCR_daily_ven[t] + Ag_daily_ven[t] * (1-p_daily_ven[t]));
   delta_ven[t,1] = pPCR_daily_ven[t] * rho * phi_PCR;
   for(i in 2:n_var) delta_ven[t,i] = rho * (phi_PCR * pPCR_daily_ven[t] + phi_Ag * (1 - pPCR_daily_ven[t]));


   

   for(i in 1:n_var) FOI_ven[t,i] = beta2_ven[t,i] * (PS_ven[t,i] + IA_ven[t,i] + IS_ven[t,i]) / S0_ven;
   
   S_ven[t+1] = S_ven[t] - S_ven[t]*sum(FOI_ven[t,]) - vac_ven[t]*S_ven[t];
   for(i in 1:n_var) E_ven[t+1,i] = E_ven[t,i] + S_ven[t] * FOI_ven[t,i] - epsilon * E_ven[t,i];
   
   for(i in 1:n_var) PS_ven[t+1,i] = PS_ven[t,i] + epsilon * E_ven[t,i] - sigma * PS_ven[t,i];
   
   for(i in 1:n_var) IA_ven[t+1,i] = IA_ven[t,i] + (1-mu) * sigma * PS_ven[t,i] - gamma * IA_ven[t,i];
   for(i in 1:n_var) IS_ven[t+1,i] = IS_ven[t,i] + (1-delta_ven[t,i]) * mu * sigma * PS_ven[t,i] - gamma * IS_ven[t,i];
   
   Q_ven[t+1] = Q_ven[t] + (delta_ven[t,1] * mu * sigma * PS_ven[t,1]) + (delta_ven[t,2] * mu * sigma * sum(PS_ven[t,2:4])) - gamma * Q_ven[t];
   
   R_ven[t+1] = R_ven[t] + gamma * Q_ven[t] + gamma * sum(IS_ven[t,])  + gamma * sum(IA_ven[t,]) + vac_ven[t] * S_ven[t];
   
   
   for(i in 1:n_var) incidence_ven[t,i] = (delta_ven[t,i]*mu*sigma*PS_ven[t,i]) / S0_ven * 100000  ;
 }
  
}



model {

  real monthly_incidence_ven[n_months_ven, n_var];  // monthly  variant incidence

  // poisson rate paramater 
  
  real lambda_M_ven[n_data_ven[1]];
  real lambda_A_ven[n_data_ven[2]]; 
  real lambda_O_ven[n_data_ven[3]];
  real lambda_Al_ven[n_data_ven[4]]; 
  

  // used for for loop
  
  int index;
  int ind ; 
  
  // index by month for daily average of variant incidences and month average concordant incidence, in month i 
   
 index = 1; 
 
 for (m in 1:n_months_ven){ 
   ind = index + 1;
   for (i in 1:n_var){
     monthly_incidence_ven[m,i] =  mean( incidence_ven[month_index_ven[index]:(month_index_ven[ind]-1),i] ) + 0.000001;    }
   index = index + 1;
 }
 
 
  
  // variant specific mean parameters for likelihood  



  lambda_M_ven = monthly_incidence_ven[index_M_ven, 1] ; 
  
  lambda_A_ven =  monthly_incidence_ven[index_A_ven, 2] ;
  
  lambda_O_ven =  monthly_incidence_ven[index_O_ven, 3];
  
  lambda_Al_ven =  monthly_incidence_ven[index_Al_ven, 4] ;
  
  
  // likelihood 
  
  target += neg_binomial_2_lpmf(y_M_ven | lambda_M_ven, 1/k);

  target += neg_binomial_2_lpmf(y_A_ven | lambda_A_ven, 1/k);

  target += neg_binomial_2_lpmf(y_O_ven | lambda_O_ven, 1/k);

  target += neg_binomial_2_lpmf(y_Al_ven | lambda_Al_ven, 1/k);




// priors

  beta  ~ lognormal(0.5,0.8);
  rho   ~ beta(2,2);
  omega ~ beta(2,2); 
  I0    ~ normal(1,100);
  k ~ exponential(0.01);
}



generated quantities {
  real R0_ven[n_var];
  for(i in 1:n_var)  R0_ven[i] = (beta[i]* (1-rho) )/ gamma;

}