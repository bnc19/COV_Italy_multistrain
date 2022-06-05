data {
  
  int<lower=1> n_locations;     // Italy and Veneto 
  int<lower=1> n_var; 
  int<lower = 1> n_data[n_locations, n_var];      // number of months observed 
  
  
  int index_M_ven[n_data[1,1]];
  int index_A_ven[n_data[1,2]];
  int index_O_ven[n_data[1,3]];
  int index_Al_ven[n_data[1,4]];
  
  int index_M_it[n_data[2,1]];
  int index_A_it[n_data[2,2]];
  int index_O_it[n_data[2,3]];
  int index_Al_it[n_data[2,4]];
  
  
  int<lower = 1> n_months[n_locations];  
  int<lower = 1> n_days[n_locations];    // number of months fitted, accounting for initial seed

  int<lower = 1> n_pop[n_locations];       // population
  int<lower = 1> n_recov[n_locations];     // seropositive 
  int<lower = 1> n_ts_ven;
  int<lower = 1> n_ts_it;
  
  int y_M_ven[n_data[1,1]];            // positive M  
  int y_A_ven[n_data[1,2]];            // positive A  
  int y_O_ven[n_data[1,3]];            // positive O
  int y_Al_ven[n_data[1,4]];          // positive Alpha 
  
  int y_M_it[n_data[2,1]];            // positive M  
  int y_A_it[n_data[2,2]];            // positive A  
  int y_O_it[n_data[2,3]];            // positive O
  int y_Al_it[n_data[2,4]];          // positive Alpha 

  real gamma;
  real sigma;
  real mu;
  real epsilon;
  
  real phi_PCR;
  real phi_Ag;

  int PCR_daily_ven[n_ts_ven];
  int Ag_daily_ven[n_ts_ven];
    
  int PCR_daily_it[n_ts_it];
  int Ag_daily_it[n_ts_it];
  
  
  real vac_ven[n_ts_ven]; 
  real vac_it[n_ts_it]; 
  
  int month_index_ven[n_months[1]+1]; // vectors to index months from days and of the number of days in each month 
  int month_index_it[n_months[2]+1]; // vectors to index months from days and of the number of days in each month 


  int time_switch1[n_locations];
  int time_switch2[n_locations];
  
  int time_seed_M[n_locations];
  int time_seed_alpha[n_locations];


}

transformed data {
  int S0[n_locations];
  for(i in 1:n_locations) S0[i]= n_pop[i]-n_recov[i];
}

parameters {
  real<lower = 0>             beta[n_var];  // transmission parameter 
  real<lower = 0>             I0[n_locations,n_var];    // seed
  real<lower = 0, upper = 1>  rho[n_locations];      // probability of reporting 
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
  
  // real S_it[(n_ts_it +1)];
  // real E_it[(n_ts_it +1), n_var];
  // real PS_it[(n_ts_it+1), n_var];
  // real IS_it[(n_ts_it+1), n_var];
  // real IA_it[(n_ts_it+1), n_var];
  // real Q_it[(n_ts_it +1)];
  // real R_it[(n_ts_it +1)];
  
  real FOI_ven[n_ts_ven, n_var]; 
  real incidence_ven[n_ts_ven,n_var];
  // 
  // real FOI_it[n_ts_it, n_var]; 
  // real incidence_it[n_ts_it,n_var];
  
  real pPCR_daily_ven[n_ts_ven];     // daily probability of PCR test
  real p_daily_ven[n_ts_ven];        // daily concordant incidence 
  real delta_ven[n_ts_ven, n_var];


  // real pPCR_daily_it[n_ts_it];     // daily probability of PCR test
  // real p_daily_it[n_ts_it];        // daily concordant incidence 
  // real delta_it[n_ts_it, n_var];


  real beta2_ven[n_ts_ven,n_var]; 
  // real beta2_it[n_ts_it,n_var]; 


// Initial values for Veneto

  S_ven[1] = S0[1] - I0[1,2] - I0[1,3] ;

 for(i in 1:n_var) E_ven[1,i] = 0 ;
 for(i in 1:n_var) PS_ven[1,i] = 0 ;
 for(i in 1:n_var) IA_ven[1,i] = 0 ;
  
  IS_ven[1,1] = 0 ;
  IS_ven[1,2] = I0[1,2] ;
  IS_ven[1,3] = I0[1,3] ;
  IS_ven[1,4] = 0 ;
  
  Q_ven[1] = 0 ;
  R_ven[1] = 0 ; 
  
  
 
// // Initial values for Italy
// 
//   S_it[1] = S0[1] - I0[2,2] - I0[2,3] ;
// 
//  for(i in 1:n_var) E_it[1,i] = 0 ;
//  for(i in 1:n_var) PS_it[1,i] = 0 ;
//  for(i in 1:n_var) IA_it[1,i] = 0 ;
//   
//   IS_it[1,1] = 0 ;
//   IS_it[1,2] = I0[2,2] ;
//   IS_it[1,3] = I0[2,3] ;
//   IS_it[1,4] = 0 ;
//   
//   Q_it[1] = 0 ;
//   R_it[1] = 0 ; 
//    
  
  
  // model interventions in Veneto
  
 for(t in 1:n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] =  beta[i];
 } 
 
 for(t in time_switch1[1] :  (time_switch2[1]-1)) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[1]) * beta[i];
 }  
 
 for(t in  time_switch2[1] :  n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[2]) * beta[i];
 }  
 
 
 
  
  // model interventions in Italy
  
 //    
 // for(t in 1:n_ts_it) {
 //   for(i in 1: n_var)  beta2_it[t,i] =  beta[i];
 // } 
 // 
 // for(t in time_switch1[2] :  (time_switch2[2]-1)) {
 //   for(i in 1: n_var)  beta2_it[t,i] = (1-omega[1]) *beta[i];
 // }  
 // 
 // for(t in  time_switch2[2] :  n_ts_it) {
 //   for(i in 1: n_var)  beta2_it[t,i] = (1-omega[2]) *beta[i];
 // }  
 // 
 

//--- Simulate in Veneto --- // 
 for (t in 1:(n_ts_ven)){
   
   IS_ven[time_seed_M[1]+1,1] = I0[1,1];
   IS_ven[time_seed_alpha[1]+1,4] = I0[1,4];

   p_daily_ven[t] = (epsilon * sum(E_ven[t,2:4])) / S0[1];
   pPCR_daily_ven[t] = (PCR_daily_ven[t] - Ag_daily_ven[t] * p_daily_ven[t]) / (PCR_daily_ven[t] + Ag_daily_ven[t] * (1-p_daily_ven[t]));
   delta_ven[t,1] = pPCR_daily_ven[t] * rho[1] * phi_PCR;
   for(i in 2:n_var) delta_ven[t,i] = rho[1] * (phi_PCR * pPCR_daily_ven[t] + phi_Ag * (1 - pPCR_daily_ven[t]));


   

   for(i in 1:n_var) FOI_ven[t,i] = beta2_ven[t,i] * (PS_ven[t,i] + IA_ven[t,i] + IS_ven[t,i]) / S0[1];
   
   S_ven[t+1] = S_ven[t] - S_ven[t]*sum(FOI_ven[t,]) - vac_ven[t]*S_ven[t];
  
   for(i in 1:n_var) E_ven[t+1,i] = E_ven[t,i] + S_ven[t] * FOI_ven[t,i] - epsilon * E_ven[t,i];
   
   for(i in 1:n_var) PS_ven[t+1,i] = PS_ven[t,i] + epsilon * E_ven[t,i] - sigma * PS_ven[t,i];
   
   for(i in 1:n_var) IA_ven[t+1,i] = IA_ven[t,i] + (1-mu) * sigma * PS_ven[t,i] - gamma * IA_ven[t,i];
   for(i in 1:n_var) IS_ven[t+1,i] = IS_ven[t,i] + (1-delta_ven[t,i]) * mu * sigma * PS_ven[t,i] - gamma * IS_ven[t,i];
   
   Q_ven[t+1] = Q_ven[t] + (delta_ven[t,1] * mu * sigma * PS_ven[t,1]) + (delta_ven[t,2] * mu * sigma * sum(PS_ven[t,2:4])) - gamma * Q_ven[t];
   
   R_ven[t+1] = R_ven[t] + gamma * Q_ven[t] + gamma * sum(IS_ven[t,])  + gamma * sum(IA_ven[t,]) + vac_ven[t] * S_ven[t];
   
   
   for(i in 1:n_var) incidence_ven[t,i] = (delta_ven[t,i]*mu*sigma*PS_ven[t,i]) / S0[1] * 100000  ;
 }
  
  
 //  //--- Simulate in Italy --- // 
 // for (t in 1:(n_ts_it)){
 //   
 //   IS_it[time_seed_M[2]+1,1] = I0[2,1];
 //   IS_it[time_seed_alpha[2]+1,4] = I0[2,4];
 // 
 //   p_daily_it[t] = (sigma * sum(E_it[t,2:4])) / S0[2];
 //   pPCR_daily_it[t] = (PCR_daily_it[t] - Ag_daily_it[t] * p_daily_it[t]) / (PCR_daily_it[t] + Ag_daily_it[t] * (1-p_daily_it[t]));
 //   delta_it[t,1] = pPCR_daily_it[t] * rho[2] * phi_PCR;
 //   for(i in 2:n_var) delta_it[t,i] = rho[2] * (phi_PCR * pPCR_daily_it[t] + phi_Ag * (1 - pPCR_daily_it[t]));
 // 
 // 
 //   
 // 
 //   for(i in 1:n_var) FOI_it[t,i] = beta2_it[t,i] * (PS_it[t,i] + IA_it[t,i] + IS_it[t,i]) / S0[2];
 //   
 //   S_it[t+1] = S_it[t] - S_it[t]*sum(FOI_it[t,]) - vac_it[t]*S_it[t];
 //   for(i in 1:n_var) E_it[t+1,i] = E_it[t,i] + S_it[t] * FOI_it[t,i] - epsilon * E_it[t,i];
 //   
 //   for(i in 1:n_var) PS_it[t+1,i] = PS_it[t,i] + epsilon * E_it[t,i] - sigma * PS_it[t,i];
 //   
 //   for(i in 1:n_var) IA_it[t+1,i] = IA_it[t,i] + (1-mu) * sigma * PS_it[t,i] - gamma * IA_it[t,i];
 //   for(i in 1:n_var) IS_it[t+1,i] = IS_it[t,i] + (1-delta_it[t,i]) * mu * sigma * PS_it[t,i] - gamma * IS_it[t,i];
 //   
 //   Q_it[t+1] = Q_it[t] + (delta_it[t,1] * mu * sigma * PS_it[t,1]) + (delta_it[t,2] * mu * sigma * sum(PS_it[t,2:4])) - gamma * Q_it[t];
 //   
 //   R_it[t+1] = R_it[t] + gamma * Q_it[t] + gamma * sum(IS_it[t,])  + gamma * sum(IA_it[t,]) + vac_it[t] * S_it[t];
 //   
 //   
 //   for(i in 1:n_var) incidence_it[t,i] = (delta_it[t,i]*mu*sigma*PS_it[t,i]) / S0[2] * 100000  ;
 // }
 //  
  
}



model {
  real monthly_incidence_ven[n_months[1], n_var];  // monthly  variant incidence
  // real monthly_incidence_it[n_months[2], n_var];  // monthly  variant incidence

  // poisson rate paramater 
  
  real lambda_M_ven[n_data[1,1]];
  real lambda_A_ven[n_data[1,2]]; 
  real lambda_O_ven[n_data[1,3]];
  real lambda_Al_ven[n_data[1,4]]; 
  
  // 
  // real lambda_M_it[n_data[2,1]];
  // real lambda_A_it[n_data[2,2]]; 
  // real lambda_O_it[n_data[2,3]];
  // real lambda_Al_it[n_data[2,4]]; 
  // 

  // used for for loop
  
  int index;
  int ind ; 
  
  // index by month for daily average of variant incidences and month average concordant incidence, in month i 
   
 index = 1; 
 
 for (m in 1:n_months[1]){ 
   ind = index + 1;
   for (i in 1:n_var){
     monthly_incidence_ven[m,i] =  mean( incidence_ven[month_index_ven[index]:(month_index_ven[ind]-1),i] );    }
   index = index + 1;
 }
 
 // 
 //     
 // index = 1; 
 // 
 // for (m in 1:n_months[2]){ 
 //   ind = index + 1;
 //   for (i in 1:n_var){
 //     monthly_incidence_it[m,i] =  mean( incidence_it[month_index_it[index]:(month_index_it[ind]-1),i] );    }
 //   index = index + 1;
 // }
 
  // variant specific mean parameters for likelihood  Veneto



  lambda_M_ven = monthly_incidence_ven[index_M_ven, 1]; 
  
  lambda_A_ven =  monthly_incidence_ven[index_A_ven, 2]; 
  
  lambda_O_ven =  monthly_incidence_ven[index_O_ven, 3];
  
  lambda_Al_ven =  monthly_incidence_ven[index_Al_ven, 4]; 
  
  
    // variant specific mean parameters for likelihood  Italy


// 
//   lambda_M_it = monthly_incidence_it[index_M_it, 1]; 
//   
//   lambda_A_it =  monthly_incidence_it[index_A_it, 2]; 
//   
//   lambda_O_it =  monthly_incidence_it[index_O_it, 3];
//   
//   lambda_Al_it =  monthly_incidence_it[index_Al_it, 4]; 
  
  // likelihood veneto 
  
  target += neg_binomial_2_lpmf(y_M_ven | lambda_M_ven, 1/k);

  target += neg_binomial_2_lpmf(y_A_ven | lambda_A_ven, 1/k);

  target += neg_binomial_2_lpmf(y_O_ven | lambda_O_ven, 1/k);

  target += neg_binomial_2_lpmf(y_Al_ven | lambda_Al_ven, 1/k);


  // // likelihood Italy 
  // 
  // target += neg_binomial_2_lpmf(y_M_it | lambda_M_it, 1/k);
  // 
  // target += neg_binomial_2_lpmf(y_A_it | lambda_A_it, 1/k);
  // 
  // target += neg_binomial_2_lpmf(y_O_it | lambda_O_it, 1/k);
  // 
  // target += neg_binomial_2_lpmf(y_Al_it | lambda_Al_it, 1/k);

// priors

  beta  ~ lognormal(0.5,0.8);
  rho   ~ beta(2,5);
  omega ~ beta(1,1); 
  for(i in 1:n_locations) I0[i]    ~ normal(1,200);
  k ~ exponential(0.01);
}



generated quantities {
  real R0[n_locations, n_var];

for(i in 1:n_locations)  
  for(v in 1:n_var)
   R0[i,v] = (beta[v]* (1-rho[i]) )/ gamma;

}