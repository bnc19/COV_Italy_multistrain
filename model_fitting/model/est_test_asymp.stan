data {
  
  // joint data 
  
  int<lower=1> n_var; 
  
  real gamma;
  real sigma;
  real mu;
  real epsilon;
  
  real phi_PCR;
  real phi_Ag;
  
  // italy data 

  int<lower = 1> n_data_it[n_var];      // number of months observed 
  
  int index_M_it[n_data_it[1]];
  int index_A_it[n_data_it[2]];
  int index_O_it[n_data_it[3]];
  int index_Al_it[n_data_it[4]];
  
  int<lower = 1> n_months_it; 
  int<lower = 1> n_days_it;    // number of months fitted, accounting for initial seed
  int<lower = 1> n_pop_it;       // population
  int<lower = 1> n_recov_it;     // seropositive 

  int<lower = 1> n_ts_it;

  int y_M_it[n_data_it[1]];            // positive M  
  int y_A_it[n_data_it[2]];            // positive A  
  int y_O_it[n_data_it[3]];            // positive O
  int y_Al_it[n_data_it[4]];           // positive Alpha 

  int PCR_daily_it[n_ts_it];
  int Ag_daily_it[n_ts_it];
  real vac_it[n_ts_it];
  
  int month_index_it[n_months_it+1]; // vectors to index months from days and of the number of days in each month 

  int time_switch1_it;
  int time_switch2_it;
  
  int time_seed_M_it;
  int time_seed_alpha_it;

 // Veneto data 
 
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
  int S0_it = n_pop_it-n_recov_it;
  int S0_ven = n_pop_ven-n_recov_ven;
}

parameters {
  real<lower = 0>             beta[n_var];  // transmission parameter 
  real<lower = 0>             I0_ven[n_var];    // seed
  real<lower = 0, upper = 1>  rho_ven;      // probability of reporting 
  real<lower = 0>             I0_it[n_var];    // seed
  real<lower = 0, upper = 1>  rho_it;      // probability of reporting 
  real<lower = 0, upper = 1>  omega[4]; // reduction in transmission
  real<lower = 0>             k;        // overdispersion parameter 

}   

transformed parameters{
  
// Define Italy param ----------------------------------------------------------
  
  real S_it[(n_ts_it +1)];
  real E_it[(n_ts_it +1), n_var];
  real PS_it[(n_ts_it+1), n_var];
  real IA_it[(n_ts_it+1), n_var];
  real Q_it[(n_ts_it +1)];
  real R_it[(n_ts_it +1)];
  
  real FOI_it[n_ts_it, n_var]; 
  real incidence_it[n_ts_it,n_var];
  
  // 
  real pPCR_daily_it[n_ts_it];     // daily probability of PCR test
  real p_daily_it[n_ts_it];        // daily concordant incidence 
  real delta_it[n_ts_it, n_var];


  real beta2_it[n_ts_it,n_var]; 
  
// Define Veneto param ---------------------------------------------------------
  
  real S_ven[(n_ts_ven +1)];
  real E_ven[(n_ts_ven +1), n_var];
  real PS_ven[(n_ts_ven+1), n_var];
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
  
// Initial values  Italy 

  S_it[1] = S0_it - I0_it[2] - I0_it[3] ;

 for(i in 1:n_var) E_it[1,i] = 0 ;
 for(i in 1:n_var) PS_it[1,i] = 0 ;

  IA_it[1,1] = 0 ;
  IA_it[1,2] = I0_it[2] ;
  IA_it[1,3] = I0_it[3] ;
  IA_it[1,4] = 0 ;
  
  Q_it[1] = 0 ;
  R_it[1] = 0 ; 
  
// Initial values Veneto -------------------------------------------------------

  S_ven[1] = S0_ven - I0_ven[2] - I0_ven[3] ;

 for(i in 1:n_var) E_ven[1,i] = 0 ;
 for(i in 1:n_var) PS_ven[1,i] = 0 ;

  IA_ven[1,1] = 0 ;
  IA_ven[1,2] = I0_ven[2] ;
  IA_ven[1,3] = I0_ven[3] ;
  IA_ven[1,4] = 0 ;
  
  Q_ven[1] = 0 ;
  R_ven[1] = 0 ; 
  
  
  // model interventions Italy -------------------------------------------------
    
 for(t in 1:n_ts_it) {
   for(i in 1: n_var)  beta2_it[t,i] =  beta[i];
 } 
 
 for(t in time_switch1_it :  (time_switch2_it-1)) {
   for(i in 1: n_var)  beta2_it[t,i] = (1-omega[1]) *beta[i];
 }  
 
 for(t in  time_switch2_it :  n_ts_it) {
   for(i in 1: n_var)  beta2_it[t,i] = (1-omega[2]) *beta[i];
 }  
 

// model interventions Veneto --------------------------------------------------
  
    
 for(t in 1:n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] =  beta[i];
 } 
 
 for(t in time_switch1_ven :  (time_switch2_ven-1)) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[3]) *beta[i];
 }  
 
 for(t in  time_switch2_ven :  n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[4]) *beta[i];
 }  
 

// Simulate Italy --------------------------------------------------------------
 for (t in 1:(n_ts_it)){
   
   IA_it[time_seed_M_it+1,1] = I0_it[1];
   IA_it[time_seed_alpha_it+1,4] = I0_it[4];

   p_daily_it[t] = (sigma * sum(E_it[t,2:4])) / S0_it;
   
   pPCR_daily_it[t] = 
   (PCR_daily_it[t] - Ag_daily_it[t] * p_daily_it[t]) / 
   (PCR_daily_it[t] + Ag_daily_it[t] * (1-p_daily_it[t]));
   
   delta_it[t,1] = pPCR_daily_it[t] * rho_it * phi_PCR;
   for(i in 2:n_var) delta_it[t,i] = 
   rho_it * (phi_PCR * pPCR_daily_it[t] + phi_Ag * (1 - pPCR_daily_it[t]));

   for(i in 1:n_var) FOI_it[t,i] =
   beta2_it[t,i] * (PS_it[t,i] + IA_it[t,i]) / S0_it;
   
   S_it[t+1] = 
   S_it[t] - S_it[t]*sum(FOI_it[t,]) - vac_it[t]*S_it[t];
   
   for(i in 1:n_var) E_it[t+1,i] = 
   E_it[t,i] + S_it[t] * FOI_it[t,i] - epsilon * E_it[t,i];
   
   for(i in 1:n_var) PS_it[t+1,i] = 
   PS_it[t,i] + epsilon * E_it[t,i] - sigma * PS_it[t,i];
   
   for(i in 1:n_var) IA_it[t+1,i] = 
   IA_it[t,i] + (1-mu) * (1-delta_it[t,i]) * sigma * PS_it[t,i] - gamma * IA_it[t,i];

   
   Q_it[t+1] = 
   Q_it[t] + (delta_it[t,1] * (1-mu) * sigma * PS_it[t,1]) 
   + (delta_it[t,2] * (1-mu) * sigma * sum(PS_it[t,2:4])) 
   + mu * sigma * sum(PS_it[t,]) - gamma * Q_it[t];
   
   R_it[t+1] = 
   R_it[t] + gamma * Q_it[t] + gamma * sum(IA_it[t,]) + vac_it[t] * S_it[t];
   
   
   for(i in 1:n_var) incidence_it[t,i] = 
  ( (delta_it[t,i]*(1-mu)*sigma*PS_it[t,i]) + (mu * sigma * PS_it[t,i]) ) / S0_it * 100000  ;
 }
 
 
 
// Simulate Veneto -------------------------------------------------------------
 for (t in 1:(n_ts_ven)){
   
   IA_ven[time_seed_M_ven+1,1] = I0_ven[1];
   IA_ven[time_seed_alpha_ven+1,4] = I0_ven[4];

   p_daily_ven[t] = (sigma * sum(E_ven[t,2:4])) / S0_ven;
   pPCR_daily_ven[t] = 
   (PCR_daily_ven[t] - Ag_daily_ven[t] * p_daily_ven[t]) / 
   (PCR_daily_ven[t] + Ag_daily_ven[t] * (1-p_daily_ven[t]));
   
   delta_ven[t,1] = pPCR_daily_ven[t] * rho_ven * phi_PCR;
   for(i in 2:n_var) delta_ven[t,i] = 
   rho_ven * (phi_PCR * pPCR_daily_ven[t] + phi_Ag * (1 - pPCR_daily_ven[t]));

   for(i in 1:n_var) FOI_ven[t,i] =
   beta2_ven[t,i] * (PS_ven[t,i] + IA_ven[t,i] ) / S0_ven;
   
   S_ven[t+1] = 
   S_ven[t] - S_ven[t]*sum(FOI_ven[t,]) - vac_ven[t]*S_ven[t];
   
   for(i in 1:n_var) E_ven[t+1,i] = 
   E_ven[t,i] + S_ven[t] * FOI_ven[t,i] - epsilon * E_ven[t,i];
   
   for(i in 1:n_var) PS_ven[t+1,i] = 
   PS_ven[t,i] + epsilon * E_ven[t,i] - sigma * PS_ven[t,i];
   
   for(i in 1:n_var) IA_ven[t+1,i] = 
   IA_ven[t,i] +(1-delta_ven[t,i]) * (1-mu) * sigma * PS_ven[t,i] - gamma * IA_ven[t,i];
   
   Q_ven[t+1] = 
   Q_ven[t] + (delta_ven[t,1] * (1-mu) * sigma * PS_ven[t,1]) 
   + (delta_ven[t,2] * (1-mu) * sigma * sum(PS_ven[t,2:4]))
   + mu * sigma * sum(PS_ven[t,]) - gamma * Q_ven[t];
   
   R_ven[t+1] = 
   R_ven[t] + gamma * Q_ven[t] + gamma * sum(IA_ven[t,]) + vac_ven[t] * S_ven[t];
   
   
   for(i in 1:n_var) incidence_ven[t,i] =
  ( (delta_ven[t,i]*(1-mu) *sigma*PS_ven[t,i]) + (mu * sigma *PS_ven[t,i]) ) / S0_ven * 100000  ;
 }
  
  
}



model {

// Define Italy model param ----------------------------------------------------


  real monthly_incidence_it[n_months_it, n_var];  // monthly variant incidence

  // poisson rate paramater 
  
  real lambda_M_it[n_data_it[1]];
  real lambda_A_it[n_data_it[2]]; 
  real lambda_O_it[n_data_it[3]];
  real lambda_Al_it[n_data_it[4]]; 
  

  // used for for loop
  
  int index_it;
  int ind_it ; 
  
  
// Define veneto model param ---------------------------------------------------
  real monthly_incidence_ven[n_months_ven, n_var];  // monthly variant incidence

  // poisson rate paramater 
  
  real lambda_M_ven[n_data_ven[1]];
  real lambda_A_ven[n_data_ven[2]]; 
  real lambda_O_ven[n_data_ven[3]];
  real lambda_Al_ven[n_data_ven[4]]; 
  

  // used for for loop
  
  int index_ven;
  int ind_ven ; 
 

// Italy  model ----------------------------------------------------------------

  // index by month for daily average of variant incidences and month average 
  // concordant incidence, in month i 
   
 index_it = 1; 
 
 for (m in 1:n_months_it){ 
   ind_it = index_it + 1;
   for (i in 1:n_var){
     monthly_incidence_it[m,i] =  
     mean( incidence_it[month_index_it[index_it]:(month_index_it[ind_it]-1),i] ) + 0.000001;}
  
   index_it = index_it + 1;
 }
 
   // variant specific mean parameters for likelihood  
  lambda_M_it = monthly_incidence_it[index_M_it, 1] ; 
  
  lambda_A_it =  monthly_incidence_it[index_A_it, 2] ;
  
  lambda_O_it =  monthly_incidence_it[index_O_it, 3];
  
  lambda_Al_it =  monthly_incidence_it[index_Al_it, 4] ;
  
  // likelihood 
  
  target += neg_binomial_2_lpmf(y_M_it | lambda_M_it, 1/k);

  target += neg_binomial_2_lpmf(y_A_it | lambda_A_it, 1/k);

  target += neg_binomial_2_lpmf(y_O_it | lambda_O_it, 1/k);

  target += neg_binomial_2_lpmf(y_Al_it | lambda_Al_it, 1/k);

// Veneto model ----------------------------------------------------------------

 // index by month for daily average of variant incidences and month average 
 // concordant incidence, in month i 
   
 index_ven= 1; 
 
 for (m in 1:n_months_ven){ 
   ind_ven = index_ven + 1;
   for (i in 1:n_var){
     monthly_incidence_ven[m,i] =  
     mean( incidence_ven[month_index_ven[index_ven]:(month_index_ven[ind_ven]-1),i] ) + 0.000001;    }
   index_ven = index_ven + 1;
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


// priors ----------------------------------------------------------------------

  beta    ~ lognormal(0.5,0.8);
  rho_it  ~ beta(1,1);
  rho_ven ~ beta(1,1);
  omega   ~ beta(1,1); 
  I0_it   ~ normal(1,200);
  I0_ven  ~ normal(1,200);
  k       ~ exponential(0.01);
}



generated quantities {
  real R0_it[n_var];
  real R0_ven[n_var];
  
 for(i in 1:n_var)  R0_it[i] = (beta[i]* (1-rho_it) )/ gamma;
 for(i in 1:n_var)  R0_ven[i] = (beta[i]* (1-rho_ven) )/ gamma;
}