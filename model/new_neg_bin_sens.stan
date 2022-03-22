functions {
  
  
  real[] SIR(real t,  // time
  real[] y,           // system state 
  real[] theta,       // parameters 
  real[] x_r,         // real data
  int [] x_i) {        // integer data
  
  
  // define integer data 
  
  int n_months = x_i[1] ;           // number of months to run the model 
  int vec[n_months+1] = x_i[2:(n_months+2)];    // position of the first of each month 
  int PCR_v[n_months] = x_i[(n_months+3):(2*n_months+2)];  // # PCR
  int Ag_v[n_months] =  x_i[(2*n_months+3):(3*n_months+2)]; // # Ag 
  int S0 =     x_i[(3*n_months+3)];    // S0
  int Rec =   x_i[(3*n_months+4)];    // recovered at day 0 
  int Days =  x_i[(3*n_months+5)];    // number of days to run the model 
  int time_switch1 = x_i[(3*n_months+6)] ; // day interventions implemented  
  int time_switch2 = x_i[(3*n_months+7)] ; // day interventions implemented  
  int seed_alpha = x_i[(3*n_months+8)] ; // date to seed alpha 
  int seed_M = x_i[(3*n_months+9)] ;    // datte to seed M varaint 

  
  
   // define rates  
  
  real vaccine[n_months] = x_r[1:n_months];  // monthly average of the daily rate of vaccination 
  real sigma = x_r[n_months + 1];    // progression rate 
  real gamma = x_r[n_months + 2];    // recovery rate
  real phi_PCR = x_r[n_months+3]; // PCR test sens
  real phi_Ag = x_r[n_months+4];  // Ag test sens 
  real rho = theta[5];    // testing rate
  
  real omega;   // intervention 
  real delta_C; // rate of detection of concordant variants
  real delta_D; // rate of detection of discortant variants
  

  
  // define states 
  
  real dS_dt ; 
  real dEM_dt ;
  real dIM_dt ;
  real dEA_dt ;
  real dIA_dt ;
  real dEO_dt ;
  real dIO_dt ;
  real dEAl_dt ;
  real dIAl_dt ;
  real dQ_dt ;
  real dR_dt ;                                          
  
  // define transmission parameters
  
  real beta_M;  
  real beta_A;
  real beta_O;
  real beta_Al;
  
  // define initial seeds
  
  real I0_M;  
  real I0_A;
  real I0_O;
  real I0_Al; 
  
  // define state variables (y)
  
   
  real S;
  real E_M;
  real I_M;
  real E_A;
  real I_A;
  real E_O;
  real I_O;
  real E_Al;
  real I_Al;
  real Q;
  real R;
  
  real N;
  
  // define time dependent variables 

  real vac ; 
  int  PCR;  
  int  Ag;
  real p;
  real pPCR;
  
  
  // define varaibles to index over for loop 

  int index ;
  int ind ;

 
    
  // initial values (varaint specific seeds)
  
  I0_M =  theta[8]; 
  I0_A =  theta[9];
  I0_O =  theta[10];
  I0_Al = theta[11]; 
  
  
  //  set initial values: 
  
  
  if (t < seed_M){
  S = y[1] ;
  } else if (t >= seed_M && t < seed_alpha){
  S = y[1] - I0_M;
  } else if (t >= seed_alpha){
    S = y[1] - I0_Al;
  }
  
  E_M = y[2] ;
  
    if (t < seed_M){
    I_M = y[3] ; 
  } else if (t >= seed_M) {    //seed M when t hits correct date 
  I_M = y[3] + I0_M;       // n.b. it says when t >= seed_m but have checked it only seeds once at time seed_M
  }
  
  E_A = y[4] ;
  I_A = y[5] ;
  E_O = y[6] ;
  I_O = y[7] ;
  E_Al = y[8] ;
  
  if (t < seed_alpha){
    I_Al = y[9] ; 
  } else if (t >= seed_alpha) {    
  I_Al = y[9] + I0_Al;
  }
  
  Q = y[10] ;
  R = y[11] ; 
  
  
  
  
  
  
// testing and vaccine is time dependent 
// vec is a vector which index by the first of each month for n_months + 1
// this allows to sequence over time and extract the right N_PCR or N_Ag or vaccine rate for the rate month 
// this for loop works when I test it outside of R, I thought it was working within Rstan but maybe its not


index = 1;

  PCR = 1;
  Ag = 0;
  vac = 0;

      
   for (x in 1:n_months){  
        
   ind = index + 1;
    
    if (t >= vec[index] && t<= (vec[ind]-1)) {
      PCR = PCR_v[x];
      Ag = Ag_v[x]; 
      vac = vaccine[x];
    } else{
      
      PCR = PCR;
      Ag = Ag;
      vac = vac;
    }
    
        index = index + 1;
        
      }
   

 
  
  // model interventions
  
  if (t < time_switch1 ) { 
    omega = 0; 
  } else if  (t >= time_switch1 && t < time_switch2) {
    omega = theta[6];
  } else if (t >= time_switch2) {
    omega = theta[7];
  }
  
  
  
  beta_M = (1- omega) * theta[1];
  beta_A = (1- omega) * theta[2];
  beta_O = (1- omega) * theta[3];
  beta_Al = (1- omega) * theta[4];
  
  
   N = S0 ;   // I tried to define N as the sum of all state variables, so that it increases when we add the late seeds of Al and M varaints but rstan won't work 
   
 
//  calculate p+ and pPCR. I assume this is redone at each step using the current values but can't be 100% certain. 
  
    p = (sigma * E_A + sigma * E_O + sigma * E_Al) / N ;
            
    pPCR = (PCR - Ag * p) / (PCR + Ag * ( 1 - p)) ;
   
    delta_C = rho * ( (phi_PCR * pPCR ) + (phi_Ag * (1 - pPCR)) ) ; 
    
    delta_D = rho * phi_PCR * pPCR ; 

   // ODE 
  
  dS_dt = - (beta_M * I_M / N  * S) - (beta_A * I_A / N * S)- (beta_O * I_O / N * S) - (beta_Al * I_Al / N * S) - vac * S; 
  
  dEM_dt = beta_M * I_M / N  * S -  sigma * E_M ;
  
  dIM_dt = (1 -  delta_D  ) * sigma * E_M - gamma * I_M;
  
  dEA_dt = beta_A * I_A / N * S -  sigma * E_A;
  
  dIA_dt = (1 - delta_C) * sigma * E_A - gamma * I_A ;
  
  dEO_dt = beta_O * I_O / N * S -  sigma * E_O;
  
  dIO_dt = (1 - delta_C) * sigma * E_O - gamma * I_O ;
  
  dEAl_dt = beta_Al * I_Al / N * S -  sigma * E_Al;
  
  dIAl_dt = (1 - delta_C) * sigma * E_Al - gamma * I_Al ;
  
  dQ_dt = delta_C * sigma * (E_A + E_Al + E_O) +   delta_D   * sigma  * E_M  - gamma * Q;
  
  dR_dt = gamma * (I_M  +  I_A + I_O + Q + I_Al) + vac * S;  
  
  
  return {dS_dt, dEM_dt, dIM_dt, dEA_dt, dIA_dt, dEO_dt,dIO_dt, dEAl_dt, dIAl_dt, dQ_dt, dR_dt};
  }
  }


data {
  
  int<lower = 1> n_data_A;      // number of months observed 
  int<lower = 1> n_data_M;
  int<lower = 1> n_data_O;
  int<lower = 1> n_data_Al;
  
  int<lower = 1> n_months; 
  int<lower = 1> n_days;    // number of months fitted, accounting for initial seed
  int<lower = 1> n_difeq;     // number of differential equations 
  int<lower = 1> n_pop;       // population
  int<lower = 1> n_recov;     // seropositive 

  real n_seq_M[n_data_M]; // number of samples sequenced 
  real n_seq_A[n_data_A]; 
  real n_seq_O[n_data_O]; 
  real n_seq_Al[n_data_Al]; 
  
  real n_reported_M[n_data_M];  // total reported cases
  real n_reported_A[n_data_A];
  real n_reported_O[n_data_O];
  real n_reported_Al[n_data_Al];

  
  int y_M[n_data_M];            // positive M  
  int y_A[n_data_A];            // positive A  
  int y_O[n_data_O];            // positive O
  int y_Al[n_data_Al];          // positive Alpha 
  
  int index_M[n_data_M];
  int index_O[n_data_O];
  int index_A[n_data_A];
  int index_Al[n_data_Al];

  real gamma;
  real sigma;
  real phi_PCR;
  real phi_Ag;

  real t0;                    // initial time point 
  real ts[n_days];          // number of time steps
  
  int PCR_daily[n_days];
  int Ag_daily[n_days];
  
  int PCR[n_months];
  int Ag[n_months];
  
  real x_r_data[n_months + 4]; // all data type real to feed into SEIR function 
  int x_i_data[3*n_months+9];  // all data type integer to feed into SEIR function 

 
  int month_index[n_months+1]; // vectors to index months from days and of the number of days in each month 
  
  real seed_mean ;  // prior mean and SD for variant seeding (so can change betwen Veneto and Italy)
  real seed_sd ; 
  real rho_a;
  real rho_b;
  
  real   rel_tol ;
  real    abs_tol ;
  real   max_num_steps;
  
}

transformed data {
  real x_r[n_months + 4] = x_r_data ;   // data for SEIR model 
  int x_i[3*n_months+9] =  x_i_data ; 
}

parameters {
  real<lower = 0>             beta[4];  // transmission parameter 
  real<lower = 0>             I0[4];    // seed
  real<lower = 0, upper = 1>  rho;      // probability of reporting 
  real<lower = 0, upper = 1>  omega[2]; // reduction in transmission
  real<lower = 0>             k;        // overdispersion parameter 
}   

transformed parameters{
  
  real y_hat[n_days, n_difeq]; // solution from the ODE solver 
  real pPCR_daily[n_days];     // daily probability of PCR test
  real p_daily[n_days];        // daily concordant incidence 
  real delta_C_daily[n_days];

    // parameters for rk45
  real theta[11] = {beta[1], beta[2], beta[3], beta[4], rho, omega[1], omega[2], I0[1], I0[2], I0[3], I0[4]} ;
  
  real init[11]  = {n_pop  - I0[2] - I0[3]  - n_recov ,  0, 0  ,   0, I0[2] ,  0, I0[3],  0,  0,  0,  n_recov };
  
 // ODE solver 
  y_hat = integrate_ode_rk45(SIR, init, t0, ts, theta, x_r, x_i); 
  
  
  
   // calculate the daily condordant varaint incidence (proportion)
  for (i in 1:n_days){
    p_daily[i] = ( sigma * y_hat[i,4]  +  sigma *y_hat[i,6] +  sigma *y_hat[i,8]) / (n_pop- n_recov);

  // calculate daily pPCR - this is definitely correct as I can check the output 
    pPCR_daily[i] = (PCR_daily[i] - Ag_daily[i] * p_daily[i]) / (PCR_daily[i] + Ag_daily[i] * (1-p_daily[i])); 
    
    delta_C_daily[i] = rho * ( (phi_PCR * pPCR_daily[i] ) + (phi_Ag * (1 - pPCR_daily[i])) ) ; 
  
 
 }  
  

  
 
  
}

model {
  
  // monthly concordant incidence, variant incidence, probability of a PCR test 
  
  real pm[n_months];                // monthly concordant incidence
  real y_hat_monthly[n_months, 4];  // monthly  variant incidence
  real pPCR[n_months];              // monthly probability of a PCR test 
  
  real delta_C[n_months];
  
// extract months to fit 
  
  real y_hat_M[n_data_M];
  real y_hat_A[n_data_A];
  real y_hat_O[n_data_O];
  real y_hat_Al[n_data_Al];
  
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

for (i in 1:n_months){ 
  ind = index + 1;
  pm[i] = mean(p_daily[month_index[index]:month_index[ind]-1]);
  for (x in 1:4){
    y_hat_monthly[i,x] =  mean( y_hat[month_index[index]:(month_index[ind]-1), (2 * x) ] ); 
  }
  index = index + 1;
}


 
  y_hat_M = y_hat_monthly[index_M, 1]; 
  
  y_hat_A =  y_hat_monthly[index_A, 2]; 
  
  y_hat_O =  y_hat_monthly[index_O, 3];
  
  y_hat_Al =  y_hat_monthly[index_Al, 4]; 
  
    
     
    
  // calculate pPCR in month i 
  
  for (i in 1:n_months){
    pPCR[i] = (PCR[i] - Ag[i] * pm[i]) / (PCR[i] + Ag[i] * (1-pm[i])); 
    
    delta_C[i] = rho * ( (phi_PCR * pPCR[i] ) + (phi_Ag * (1 - pPCR[i])) ) ; 
  }
  
    
    
    

  
  // variant specific mean parameters for likelihood  

 for (i in 1:n_data_M){  
    lambda_M[i] = (rho  * pPCR[i+(n_months-n_data_M)] * phi_PCR  * sigma * y_hat_M[i]) / n_reported_M[i] * n_seq_M[i] + 0.0001;
  }
  
  for (i in 1:n_data_A){  
    lambda_A[i] = (delta_C[i + (n_months - n_data_A)] *  sigma * y_hat_A[i]) / n_reported_A[(i)] * n_seq_A[(i)]+ 0.0001;
  }
  
  for (i in 1:n_data_O){  
    lambda_O[i] = (delta_C[i + (n_months - n_data_O)] *  sigma * y_hat_O[i]) / n_reported_O[(i)] * n_seq_O[(i)]+ 0.0001;
  }
  
  for (i in 1:n_data_Al){  
    lambda_Al[i] = (delta_C[i + (n_months - n_data_Al)] * sigma * y_hat_Al[i]) / n_reported_Al[(i)] * n_seq_Al[(i)]+ 0.0001;
  }
  
  
  // likelihood 
  
  target += neg_binomial_2_lpmf(y_M | lambda_M, 1/ k);

  target += neg_binomial_2_lpmf(y_A | lambda_A,1/ k);

  target += neg_binomial_2_lpmf(y_O | lambda_O, 1/k);

  target += neg_binomial_2_lpmf(y_Al | lambda_Al,1/k);




// priors

  beta[1:4]  ~ normal(2,3);
  rho        ~ beta(rho_a,rho_b);
  omega[1:2] ~ beta(2,2); // uniform
  I0[1:4]    ~ normal(seed_mean,seed_sd);
  k ~ exponential(0.01);


// 
// beta[1:4]  ~ lognormal(2,1);
// rho        ~ beta(rho_a,rho_b);
// omega[1:2] ~ beta(2,2); // uniform
// I0[1:4]    ~ normal(seed_mean,seed_sd);
// k ~ exponential(0.01);

}



generated quantities {


// these values are estimated from posterior after discarding burn in 

  real R_0[4];
  real reported_incidence[n_days,4];
  real true_incidence [n_days,4];
  real total_incidence[n_days];
  real susceptible[n_days];
  real recovered[n_days];






 
  // monthly concordant incidence, variant incidence, probability of a PCR test 
  
  real pm[n_months];                // monthly concordant incidence
  real y_hat_monthly[n_months, 4];  // monthly  variant incidence
  real pPCR[n_months];              // monthly probability of a PCR test 
  
  real delta_C[n_months];
  
// extract months to fit 
  
  real y_hat_M[n_data_M];
  real y_hat_A[n_data_A];
  real y_hat_O[n_data_O];
  real y_hat_Al[n_data_Al];
  
  // poisson rate paramater 
  
  real lambda_M[n_data_M];
  real lambda_A[n_data_A]; 
  real lambda_O[n_data_O];
  real lambda_Al[n_data_Al]; 
  

  // used for for loop
  
  int index;
  int ind ; 
  
  // index by month for daily average of variant incidences and month average concordant incidence, in month i 
   
   
   // calculate LL 


  vector[(n_data_M + n_data_A + n_data_O + n_data_Al)] log_lik;
   
   
     // recalculate lambda as not recognised from model block 
   
index = 1; 

for (i in 1:n_months){ 
  ind = index + 1;
  pm[i] = mean(p_daily[month_index[index]:month_index[ind]-1]);
  for (x in 1:4){
    y_hat_monthly[i,x] =  mean( y_hat[month_index[index]:(month_index[ind]-1), (2 * x) ] ); 
  }
  index = index + 1;
}


 
  y_hat_M = y_hat_monthly[index_M, 1]; 
  
  y_hat_A =  y_hat_monthly[index_A, 2]; 
  
  y_hat_O =  y_hat_monthly[index_O, 3];
  
  y_hat_Al =  y_hat_monthly[index_Al, 4]; 
  
    
     
    
  // calculate pPCR in month i 
  
  for (i in 1:n_months){
    pPCR[i] = (PCR[i] - Ag[i] * pm[i]) / (PCR[i] + Ag[i] * (1-pm[i])); 
    
    delta_C[i] = rho * ( (phi_PCR * pPCR[i] ) + (phi_Ag * (1 - pPCR[i])) ) ; 
  }
  
    
    
    

  
  // variant specific mean parameters for likelihood  

 for (i in 1:n_data_M){  
    lambda_M[i] = (rho  * pPCR[i+(n_months-n_data_M)] * phi_PCR  * sigma * y_hat_M[i]) / n_reported_M[i] * n_seq_M[i]+ 0.0001;
  }
  
  for (i in 1:n_data_A){  
    lambda_A[i] = (delta_C[i + (n_months - n_data_A)] *  sigma * y_hat_A[i]) / n_reported_A[(i)] * n_seq_A[(i)]+ 0.0001;
  }
  
  for (i in 1:n_data_O){  
    lambda_O[i] = (delta_C[i + (n_months - n_data_O)] *  sigma * y_hat_O[i]) / n_reported_O[(i)] * n_seq_O[(i)]+ 0.0001;
  }
  
  for (i in 1:n_data_Al){  
    lambda_Al[i] = (delta_C[i + (n_months - n_data_Al)] * sigma * y_hat_Al[i]) / n_reported_Al[(i)] * n_seq_Al[(i)]+ 0.0001;
  }
  
  
  
  // log likelihood
  
  
  
  // LL
  
  
 for (i in 1:n_data_M){
    log_lik[i] = neg_binomial_2_lpmf(y_M[i] | lambda_M[i], k);
  }


  for(i in (n_data_M+1):(n_data_M+n_data_A)){
    log_lik[i] = neg_binomial_2_lpmf(y_A[(i- n_data_M)] | lambda_A[(i- n_data_M)], k);
  }


  for(i in (n_data_M+ n_data_A +1):(n_data_M+n_data_A + n_data_O)){
    log_lik[i] = neg_binomial_2_lpmf(y_O [(i - n_data_M- n_data_A)] | lambda_O[(i-n_data_M-n_data_A)], k);
  }
  
  
  for(i in (n_data_M+ n_data_A + n_data_O +1):(n_data_M+n_data_A + n_data_O + n_data_Al)){
    log_lik[i] = neg_binomial_2_lpmf(y_Al [(i - n_data_M - n_data_A - n_data_O)] | lambda_Al[(i-n_data_M - n_data_A - n_data_O)], k);
  }
  
  
  // generated quantities 
  
  for (i in 1: n_days){
    
    for (x in 1:4){
    
    if (x == 1){  
      reported_incidence[i,x] = y_hat[i,(2*x)] * sigma * rho * pPCR_daily[i] * phi_PCR; // M variant reported incidence 
      } else {
        
      reported_incidence[i,x] = y_hat[i,(2*x)] * sigma * delta_C_daily[i] ;  // A, O, Al variant reported incidence 
      
      }
      
      true_incidence[i,x] = y_hat[i,(2*x)] * sigma ;  // true (unseen) incidence of all variants  
      
    } 
    
    
    total_incidence[i] = true_incidence[i, 1] + true_incidence[i,2]  + true_incidence[i,3] + true_incidence[i,4];  // total incidence, not variant specific
    susceptible[i] = y_hat[i,1];
    recovered[i]  = y_hat[i,11];
  }


for (i in 1:4){  
  R_0[i] =  beta[i] * (1 - rho) / gamma ;
}



}