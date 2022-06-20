data {
  
  // joint data 
  
  int<lower=1> n_var; 
  int scale_time_step;

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
  
  // parameters
  
  real<lower = 0>             beta[n_var];  // transmission parameter 
  real<lower = 0>             I0_ven[n_var];    // seed
  real<lower = 0, upper = 1>  rho_ven;      // probability of reporting 
  real<lower = 0>             I0_it[n_var];    // seed
  real<lower = 0, upper = 1>  rho_it;      // probability of reporting 
  real<lower = 0, upper = 1>  omega[4]; // reduction in transmission
  real<lower = 0>             k;        // overdispersion parameter 
  real<lower = gamma>         tau;

}

transformed data {
  int S0_it = n_pop_it-n_recov_it;
  int S0_ven = n_pop_ven-n_recov_ven;
  real gamma2 = gamma / scale_time_step; 
  real sigma2 = sigma / scale_time_step; 
  real epsilon2 = epsilon / scale_time_step; 
}

 
model {
}


generated quantities {
  
    
 real beta2[n_var] ;  
 real tau2 = tau / scale_time_step;
  
// Define Italy param ----------------------------------------------------------
  
  real S_it[(n_ts_it +1)];
  real E_it[(n_ts_it +1), n_var];
  real PS_it[(n_ts_it+1), n_var];
  real IA_it[(n_ts_it+1), n_var];
  real IS_it[(n_ts_it+1), n_var];
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
  real IS_ven[(n_ts_ven+1), n_var];
  real Q_ven[(n_ts_ven +1)];
  real R_ven[(n_ts_ven +1)];
  
  real FOI_ven[n_ts_ven, n_var]; 
  real incidence_ven[n_ts_ven,n_var];
  
  // 
  real pPCR_daily_ven[n_ts_ven];     // daily probability of PCR test
  real p_daily_ven[n_ts_ven];        // daily concordant incidence 
  real delta_ven[n_ts_ven, n_var];


  real beta2_ven[n_ts_ven,n_var]; 
  
   vector[sum(n_data_it) + sum(n_data_ven)] log_lik;
   
   real sum_LL ;
  
// Define Italy model param ----------------------------------------------------
  real daily_incidence_it[n_days_it, n_var];  // monthly variant incidence
  real monthly_incidence_it[n_months_it, n_var];  // monthly variant incidence

  // poisson rate paramater 
  
  real lambda_M_it[n_data_it[1]];
  real lambda_A_it[n_data_it[2]]; 
  real lambda_O_it[n_data_it[3]];
  real lambda_Al_it[n_data_it[4]]; 
  

  // used for for loop
  
  int index_it;
  int ind_it ; 
  int index_i;
  int ind_i ;
  
  
// Define veneto model param ---------------------------------------------------
 
  real daily_incidence_ven[n_days_ven, n_var];  // monthly variant incidence
  real monthly_incidence_ven[n_months_ven, n_var];  // monthly variant incidence

  // poisson rate paramater 
  
  real lambda_M_ven[n_data_ven[1]];
  real lambda_A_ven[n_data_ven[2]]; 
  real lambda_O_ven[n_data_ven[3]];
  real lambda_Al_ven[n_data_ven[4]]; 
  

  // used for for loop
  
  int index_ven;
  int ind_ven ; 
  int index_v;
  int ind_v ;
  
   // Scale beta 

  for(i in 1:n_var) beta2[i] = beta[i] / scale_time_step;
  
// Initial values  Italy 

  S_it[1] = S0_it - I0_it[2] - I0_it[3] ;

 for(i in 1:n_var) E_it[1,i] = 0 ;
 for(i in 1:n_var) PS_it[1,i] = 0 ;
 for(i in 1:n_var) IS_it[1,i] = 0 ;

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
 for(i in 1:n_var) IS_ven[1,i] = 0 ; 
  
  IA_ven[1,1] = 0 ;
  IA_ven[1,2] = I0_ven[2] ;
  IA_ven[1,3] = I0_ven[3] ;
  IA_ven[1,4] = 0 ;
  
  Q_ven[1] = 0 ;
  R_ven[1] = 0 ; 
  
  
  // model interventions Italy -------------------------------------------------
    
 for(t in 1:n_ts_it) {
   for(i in 1: n_var)  beta2_it[t,i] =  beta2[i];
 } 
 
 for(t in time_switch1_it :  (time_switch2_it-1)) {
   for(i in 1: n_var)  beta2_it[t,i] = (1-omega[1]) *beta2[i];
 }  
 
 for(t in  time_switch2_it :  n_ts_it) {
   for(i in 1: n_var)  beta2_it[t,i] = (1-omega[2]) *beta2[i];
 }  
 

// model interventions Veneto --------------------------------------------------
  
    
 for(t in 1:n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] =  beta2[i];
 } 
 
 for(t in time_switch1_ven :  (time_switch2_ven-1)) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[3]) *beta2[i];
 }  
 
 for(t in  time_switch2_ven :  n_ts_ven) {
   for(i in 1: n_var)  beta2_ven[t,i] = (1-omega[4]) *beta2[i];
 }  
 

// Simulate Italy --------------------------------------------------------------

  for (t in 1:(n_ts_it)){
   
   IA_it[time_seed_M_it+1,1] = I0_it[1];
   IA_it[time_seed_alpha_it+1,4] = I0_it[4];

   p_daily_it[t] = (epsilon2 * sum(E_it[t,2:4])) / S0_it;
   
   pPCR_daily_it[t] = 
   (PCR_daily_it[t] - Ag_daily_it[t] * p_daily_it[t]) / 
   (PCR_daily_it[t] + Ag_daily_it[t] * (1-p_daily_it[t]));
   
   delta_it[t,1] = pPCR_daily_it[t] * rho_it * phi_PCR;
   for(i in 2:n_var) delta_it[t,i] = 
   rho_it * (phi_PCR * pPCR_daily_it[t] + phi_Ag * (1 - pPCR_daily_it[t]));

   for(i in 1:n_var) FOI_it[t,i] =
   beta2_it[t,i] * (PS_it[t,i] + IA_it[t,i] + IS_it[t,i]) / S0_it;
   
   S_it[t+1] = 
   S_it[t] - S_it[t]*sum(FOI_it[t,]) - vac_it[t]*S_it[t];
   
   for(i in 1:n_var) E_it[t+1,i] = 
   E_it[t,i] + S_it[t] * FOI_it[t,i] - epsilon2 * E_it[t,i];
   
   for(i in 1:n_var) PS_it[t+1,i] = 
   PS_it[t,i] + epsilon2 * E_it[t,i] - sigma2 * PS_it[t,i];
   
   for(i in 1:n_var) IA_it[t+1,i] = 
   IA_it[t,i] + (1-mu) * sigma2 * PS_it[t,i] 
   - gamma2 * IA_it[t,i];

    for(i in 1:n_var) IS_it[t+1,i] = 
   IS_it[t,i] + mu * sigma2 * PS_it[t,i] 
   - gamma2  * (1-delta_it[t,i])  * IS_it[t,i] 
   - tau2 * delta_it[t,i] * IS_it[t,i]  ;
   
   Q_it[t+1] = 
   Q_it[t] + (delta_it[t,1] * tau2 * IS_it[t,1]) 
   + (delta_it[t,2] * tau2 * sum(IS_it[t,2:4])) 
   - gamma2 * Q_it[t];
   
   R_it[t+1] = 
   R_it[t] + gamma2 * Q_it[t] +  gamma2 *  sum(IA_it[t,])
   + gamma2  * (1-delta_it[t,1])  * IS_it[t,1]
   + gamma2  * (1-delta_it[t,2])  * sum(IS_it[t,2:4])  
   + vac_it[t] * S_it[t];
   
   
   for(i in 1:n_var) incidence_it[t,i] = 
  ( (delta_it[t,i] * tau2 * IS_it[t,i]) ) / S0_it * 100000  ;
 }
 
 
  
// Simulate Veneto -------------------------------------------------------------
  for (t in 1:(n_ts_ven)){
   
   IA_ven[time_seed_M_ven+1,1] = I0_ven[1];
   IA_ven[time_seed_alpha_ven+1,4] = I0_ven[4];

   p_daily_ven[t] = (epsilon2 * sum(E_ven[t,2:4])) / S0_ven;
   
   pPCR_daily_ven[t] = 
   (PCR_daily_ven[t] - Ag_daily_ven[t] * p_daily_ven[t]) / 
   (PCR_daily_ven[t] + Ag_daily_ven[t] * (1-p_daily_ven[t]));
   
   delta_ven[t,1] = pPCR_daily_ven[t] * rho_ven * phi_PCR;
   for(i in 2:n_var) delta_ven[t,i] = 
   rho_ven * (phi_PCR * pPCR_daily_ven[t] + phi_Ag * (1 - pPCR_daily_ven[t]));

   for(i in 1:n_var) FOI_ven[t,i] =
   beta2_ven[t,i] * (PS_ven[t,i] + IA_ven[t,i] + IS_ven[t,i]) / S0_ven;
   
   S_ven[t+1] = 
   S_ven[t] - S_ven[t]*sum(FOI_ven[t,]) - vac_ven[t]*S_ven[t];
   
   for(i in 1:n_var) E_ven[t+1,i] = 
   E_ven[t,i] + S_ven[t] * FOI_ven[t,i] - epsilon2 * E_ven[t,i];
   
   for(i in 1:n_var) PS_ven[t+1,i] = 
   PS_ven[t,i] + epsilon2 * E_ven[t,i] - sigma2 * PS_ven[t,i];
   
   for(i in 1:n_var) IA_ven[t+1,i] = 
   IA_ven[t,i] + (1-mu) * sigma2 * PS_ven[t,i] 
   - gamma2 * IA_ven[t,i];

    for(i in 1:n_var) IS_ven[t+1,i] = 
   IS_ven[t,i] + mu * sigma2 * PS_ven[t,i] 
   - gamma2  * (1-delta_ven[t,i])  * IS_ven[t,i] 
   - tau2 * delta_ven[t,i] * IS_ven[t,i]  ;
   
   Q_ven[t+1] = 
   Q_ven[t] + (delta_ven[t,1] * tau2 * IS_ven[t,1]) 
   + (delta_ven[t,2] * tau2 * sum(IS_ven[t,2:4])) 
   - gamma2 * Q_ven[t];
   
   R_ven[t+1] = 
   R_ven[t] + gamma2 * Q_ven[t] +  gamma2 *  sum(IA_ven[t,])
   + gamma2  * (1-delta_ven[t,1])  * IS_ven[t,1]
   + gamma2  * (1-delta_ven[t,2])  * sum(IS_ven[t,2:4])  
   + vac_ven[t] * S_ven[t];
   
   
   
   for(i in 1:n_var) incidence_ven[t,i] = 
  ((delta_ven[t,i] * tau2 * IS_ven[t,i]) ) / S0_ven * 100000  ;
 }
 
 
// Italy  model ----------------------------------------------------------------

  index_i = 1;

  for (t in 1:n_days_it){
  ind_i = index_i + (scale_time_step-1);

  for(i in 1:n_var) daily_incidence_it[t,i] =  sum(incidence_it[index_i:ind_i,i]);

  index_i = index_i + scale_time_step;
}

  // index by month for daily average of variant incidences and month average 
  // concordant incidence, in month i 
   
 index_it = 1; 
 
 for (m in 1:n_months_it){ 
   ind_it = index_it + 1;
   for (i in 1:n_var){
     monthly_incidence_it[m,i] =  
     mean( daily_incidence_it[month_index_it[index_it]:(month_index_it[ind_it]-1),i] ) + 0.000001;}
  
   index_it = index_it + 1;
 }
 
   // variant specific mean parameters for likelihood  
  lambda_M_it = monthly_incidence_it[index_M_it, 1] ; 
  
  lambda_A_it =  monthly_incidence_it[index_A_it, 2] ;
  
  lambda_O_it =  monthly_incidence_it[index_O_it, 3];
  
  lambda_Al_it =  monthly_incidence_it[index_Al_it, 4] ;


// Veneto model ----------------------------------------------------------------

 // index by month for daily average of variant incidences and month average 
 // concordant incidence, in month i 
 
   index_v = 1;

  for (t in 1:n_days_ven){
  ind_v = index_v + (scale_time_step-1);

  for(i in 1:n_var) daily_incidence_ven[t,i] =  sum(incidence_ven[index_v:ind_v,i]);

  index_v = index_v + scale_time_step;
}

   
 index_ven= 1; 
 
 for (m in 1:n_months_ven){ 
   ind_ven = index_ven + 1;
   for (i in 1:n_var){
     monthly_incidence_ven[m,i] =  
     mean( daily_incidence_ven[month_index_ven[index_ven]:(month_index_ven[ind_ven]-1),i] ) + 0.000001;    }
   index_ven = index_ven + 1;
 }
 
  // variant specific mean parameters for likelihood  

  lambda_M_ven = monthly_incidence_ven[index_M_ven, 1] ; 
  
  lambda_A_ven =  monthly_incidence_ven[index_A_ven, 2] ;
  
  lambda_O_ven =  monthly_incidence_ven[index_O_ven, 3];
  
  lambda_Al_ven =  monthly_incidence_ven[index_Al_ven, 4] ;
  
  
  // Log likelihood 


 for (i in 1:n_data_it[1]){
    log_lik[i] =10 * neg_binomial_2_lpmf(y_M_it [i] | lambda_M_it [i], k);
  }
  
 for (i in (n_data_it[1] + 1): ( sum(n_data_it[1:2]))){
    log_lik[i] = neg_binomial_2_lpmf(y_A_it [(i-n_data_it[1])] | lambda_A_it[(i-n_data_it[1])], k);
  }

 for (i in (sum(n_data_it[1:2])+1): sum(n_data_it[1:3]) ){
    log_lik[i] = neg_binomial_2_lpmf(y_O_it [(i- sum(n_data_it[1:2]))] | lambda_O_it[(i- sum(n_data_it[1:2]))], k);
  }

 for (i in (sum(n_data_it[1:3])+1): (sum(n_data_it))){
    log_lik[i] = neg_binomial_2_lpmf(y_Al_it [(i- sum(n_data_it[1:3]))] |lambda_Al_it[(i-sum(n_data_it[1:3]))], k);
  }

 for (i in ( sum(n_data_it) +1):(sum(n_data_it) + n_data_ven[1] )){
    log_lik[i] = 10 * neg_binomial_2_lpmf(y_M_ven [(i-sum(n_data_it))] | lambda_M_ven[(i-sum(n_data_it))], k);
  }

 for (i in ( sum(n_data_it) + n_data_ven[1]+1):(sum(n_data_it) + n_data_ven[1] + n_data_ven[2] )){
    log_lik[i] = neg_binomial_2_lpmf(y_A_ven [(i-sum(n_data_it) -  n_data_ven[1])] | 
    lambda_A_ven[(i-sum(n_data_it)- n_data_ven[1])], k);
  }
  
 for (i in ( sum(n_data_it) + sum(n_data_ven[1:2])+1):(sum(n_data_it) + sum(n_data_ven[1:3]) )){
    log_lik[i] = neg_binomial_2_lpmf(y_O_ven [(i-sum(n_data_it) - sum(n_data_ven[1:2]))] | 
    lambda_O_ven[(i-sum(n_data_it)- sum(n_data_ven[1:2]))], k);
  }  
  
   for (i in ( sum(n_data_it) + sum(n_data_ven[1:3])+1):(sum(n_data_it) + sum(n_data_ven) )){
    log_lik[i] = neg_binomial_2_lpmf(y_Al_ven [(i-sum(n_data_it) -  sum(n_data_ven[1:3]))] | 
    lambda_Al_ven[(i-sum(n_data_it)- sum(n_data_ven[1:3]))], k);
  }  

sum_LL = sum(log_lik); 

}