data {
  
  int<lower = 1> n_data_A;      // number of months observed 
  int<lower = 1> n_data_M;
  int<lower = 1> n_data_O;
  int<lower = 1> n_data_Al;
  
  
  int index_M[n_data_M];
  int index_O[n_data_O];
  int index_A[n_data_A];
  int index_Al[n_data_Al];
  
  int<lower = 1> n_months; 
  int<lower = 1> n_days;    // number of months fitted, accounting for initial seed
  int<lower = 1> n_pop;       // population
  int<lower = 1> n_recov;     // seropositive 
  int<lower = 1> n_variants;
  int<lower = 1> n_ts;

  
  int y_M[n_data_M];            // positive M  
  int y_A[n_data_A];            // positive A  
  int y_O[n_data_O];            // positive O
  int y_Al[n_data_Al];          // positive Alpha 
  

  real gamma;
  real sigma;
  real phi_PCR;
  real phi_Ag;

  int PCR_daily[n_ts];
  int Ag_daily[n_ts];
    
  real vac[n_ts]; 
  
  int month_index[n_months+1]; // vectors to index months from days and of the number of days in each month 

  int time_switch1;
  int time_switch2;
  
  int time_seed_M;
  int time_seed_alpha;


}

transformed data {
  int S0 = n_pop-n_recov;
}

parameters {
  real<lower = 0>             beta[n_variants];  // transmission parameter 
  real<lower = 0>             I0[n_variants];    // seed
  real<lower = 0, upper = 1>  rho;      // probability of reporting 
  real<lower = 0, upper = 1>  omega[2]; // reduction in transmission
  real<lower = 0>             k;        // overdispersion parameter 
}   

transformed parameters{
  
  real S[(n_ts+1)];
  real E[(n_ts+1), n_variants];
  real I[(n_ts+1), n_variants];
  real Q[(n_ts+1)];
  real R[(n_ts+1)];
  
  
  real FOI[n_ts, n_variants]; 
  real incidence[n_ts,n_variants];
  
  real pPCR_daily[n_ts];     // daily probability of PCR test
  real p_daily[n_ts];        // daily concordant incidence 
  real delta[n_ts, n_variants];

  real beta2[n_ts,n_variants]; 

// Initial values 

  S[1] = S0 - I0[2] - I0[3] ;

 for(i in 1:n_variants) E[1,i] = 0 ;

  I[1,1] = 0 ;
  I[1,2] = I0[2] ;
  I[1,3] = I0[3] ;
  I[1,4] = 0 ;
  
  Q[1] = 0 ;
  R[1] = 0 ; 
  
  // model interventions 
  
    
 for(t in 1:n_ts) {
   for(i in 1: n_variants)  beta2[t,i] =  beta[i];
 } 
 
 for(t in time_switch1 :  (time_switch2-1)) {
   for(i in 1: n_variants)  beta2[t,i] = (1-omega[1]) *beta[i];
 }  
 
 for(t in  time_switch2 :  n_ts) {
   for(i in 1: n_variants)  beta2[t,i] = (1-omega[2]) *beta[i];
 }  
 

//--- Simulate --- // 
 for (t in 1:(n_ts)){
   
   I[time_seed_M+1,1] = I0[1];
   I[time_seed_alpha+1,4] = I0[4];

   p_daily[t] = (sigma * sum(E[t,2:4])) / S0;
   pPCR_daily[t] = (PCR_daily[t] - Ag_daily[t] * p_daily[t]) / (PCR_daily[t] + Ag_daily[t] * (1-p_daily[t]));
   delta[t,1] = pPCR_daily[t] * rho * phi_PCR;
   for(i in 2:n_variants) delta[t,i] = rho * (phi_PCR * pPCR_daily[t] + phi_Ag * (1 - pPCR_daily[t]));


   

   for(i in 1:n_variants) FOI[t,i] = beta2[t,i] * I[t,i] / S0;
   
   S[t+1] = S[t] - S[t]*sum(FOI[t,]) - vac[t]*S[t];
   for(i in 1:n_variants) E[t+1,i] = E[t,i] + S[t] * FOI[t,i] - sigma * E[t,i];
   for(i in 1:n_variants) I[t+1,i] = I[t,i] + (1-delta[t,i]) * sigma * E[t,i] - gamma * I[t,i];
   Q[t+1] = Q[t] + (delta[t,1] * sigma * E[t,1]) + (delta[t,2] * sigma * sum(E[t,2:4])) - gamma * Q[t];
   R[t+1] = R[t] + gamma * Q[t] + gamma * sum(I[t,]) + vac[t] * S[t];
   
   
   for(i in 1:n_variants) incidence[t,i] = (delta[t,i]*sigma*E[t,i]) / S0 * 100000  ;
 }
  
}



model {

  real monthly_incidence[n_months, n_variants];  // monthly  variant incidence

  // poisson rate paramater 
  
  real lambda_M[n_data_M];
  real lambda_A[n_data_A]; 
  real lambda_O[n_data_O];
  real lambda_Al[n_data_Al]; 
  

  // used for for loop
  
  int index;
  int ind ; 
  
  // index by month for daily average of variant incidences and month average concordant incidence, in month i 
   
 index = 1; 
 
 for (m in 1:n_months){ 
   ind = index + 1;
   for (i in 1:n_variants){
     monthly_incidence[m,i] =  mean( incidence[month_index[index]:(month_index[ind]-1),i] );    }
   index = index + 1;
 }
 
 
  
  // variant specific mean parameters for likelihood  



  lambda_M = monthly_incidence[index_M, 1]; 
  
  lambda_A =  monthly_incidence[index_A, 2]; 
  
  lambda_O =  monthly_incidence[index_O, 3];
  
  lambda_Al =  monthly_incidence[index_Al, 4]; 
  
  
  // likelihood 
  
  target += neg_binomial_2_lpmf(y_M | lambda_M, 1/k);

  target += neg_binomial_2_lpmf(y_A | lambda_A, 1/k);

  target += neg_binomial_2_lpmf(y_O | lambda_O, 1/k);

  target += neg_binomial_2_lpmf(y_Al | lambda_Al, 1/k);




// priors

  beta  ~ lognormal(0.5,0.8);
  rho   ~ beta(2,2);
  omega ~ beta(2,2); 
  I0    ~ normal(1,200);
k ~ exponential(0.01);
}



generated quantities {
  real R0[n_variants];
  
  real monthly_incidence[n_months, n_variants];  // monthly  variant incidence

// used for for loop
  
  int index;
  int ind ; 
  
  // index by month for daily average of variant incidences and month average concordant incidence, in month i 
   
 index = 1; 
 
 for (m in 1:n_months){ 
   ind = index + 1;
   for (i in 1:n_variants){
     monthly_incidence[m,i] =  mean( incidence[month_index[index]:(month_index[ind]-1),i] );    }
   index = index + 1;
 }
 
 
 for(i in 1:n_variants)  R0[i] = (beta[i]* (1-rho) )/ gamma;

}