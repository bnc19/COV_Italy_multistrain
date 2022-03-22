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
  real rho = x_r[n_months + 3];    // testing rate
 
  real phi_PCR = x_r[n_months+14]; // PCR test sens
  real phi_Ag = x_r[n_months+15];  // Ag test sens
  
  real omega;   // intervention 
  real delta_C; // rate of detection of concordant variant 
  real delta_D; // rate of detection of discordant variant 
  
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
  
  I0_M =  x_r[n_months + 10];
  I0_A =  x_r[n_months + 11];
  I0_O =  x_r[n_months + 12];
  I0_Al = x_r[n_months + 13];
  
  
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
  } else if (t >= seed_alpha) {    //seed alpha October 1st 
  I_Al = y[9] + I0_Al;
  }
  
  Q = y[10] ;
  R = y[11] ; 
  
  
  
  
  
  
// testing and vaccine is time dependent 
// vec is a vector which index by the first of each month for n_months + 1
// this allows to sequence over time and extract the right N_PCR or N_Ag or vaccine rate for the rate month 
// this for loop works when I test it outside of R, I thought it was working within Rstan but maybe its not


index = 1;

  vac = 0;

      
   for (x in 1:n_months){  
        
   ind = index + 1;
    
    if (t >= vec[index] && t<= (vec[ind]-1)) {
      vac = vaccine[x];
    } else{
      vac = vac;
    }
    
        index = index + 1;
        
      }
   

 
  
  // model interventions
  
  if (t < time_switch1 ) { 
    omega = 0; 
  } else if  (t >= time_switch1 && t < time_switch2) {
    omega = x_r[n_months + 8];
  } else if (t >= time_switch2) {
    omega = x_r[n_months + 9];
  }
  
  
  
  beta_M = (1- omega) * x_r[n_months + 4];
  beta_A = (1- omega) * x_r[n_months + 5];
  beta_O = (1- omega) * x_r[n_months + 6];
  beta_Al = (1- omega) * x_r[n_months + 7];
  
  
   N = S0 ;   
   
 
//  calculate p+ and pPCR
  
    p = (sigma * E_A + sigma * E_O + sigma * E_Al) / N ;
    
    pPCR = 0;
   
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
  int<lower = 1> n_months; 
  int<lower = 1> n_days;    // number of months fitted, accounting for initial seed
  int<lower = 1> n_difeq;     // number of differential equations 
  int<lower = 1> n_pop;       // population
  int<lower = 1> n_recov;     // seropositive 

  real gamma;
  real sigma;
  real t0;                    // initial time point 
  real ts[n_days];          // number of time steps
  real beta [4];
  real rho ; 
  real seed[4];
  real phi_PCR;
  real phi_Ag;

  
  int PCR_daily[n_days];
  int Ag_daily[n_days];
  
  
  real x_r_data[n_months + 15]; // all data type real to feed into SEIR function 
  int x_i_data[3*n_months+9];  // all data type integer to feed into SEIR function 

}

transformed data {
  real x_r[n_months + 15] = x_r_data ;   // data for SEIR model 
  int x_i[3*n_months+9] =  x_i_data ; 
}



transformed parameters{
  // parameters for rk45
  real theta[0] ;
}



  



generated quantities {




  real R_0[n_days, 4];
  real reported_incidence[n_days,4];
  real true_incidence [n_days,4];
  real total_incidence[n_days];
  real total_reported_incidence[n_days];
  // real susceptible[n_days];
  // real recovered[n_days];


  real y_hat[n_days, n_difeq]; // solution from the ODE solver 
  real pPCR_daily[n_days];     // daily probability of PCR test
  real delta_C_daily[n_days];



  real init[11]  = {n_pop  - seed[2] - seed[3]  - n_recov ,  0, 0  ,   0, seed[2] ,  0, seed[3],  0,  0,  0,  n_recov };
 // ODE solver 
  y_hat = integrate_ode_rk45(SIR, init, t0, ts, theta, x_r, x_i); 
  
  
 
  
   
 // calculate the daily condordant varaint incidence (proportion)
  for (i in 1:n_days){
 
  pPCR_daily[i] = 0.0;
  delta_C_daily[i] = rho * ( (phi_PCR * pPCR_daily[i] ) + (phi_Ag * (1 - pPCR_daily[i])) ) ; 
  
  }
  
  
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
    total_reported_incidence[i] = reported_incidence[i, 1] + reported_incidence[i,2]  + reported_incidence[i,3] + reported_incidence[i,4];  // total incidence, not variant specific

  // susceptible[i] = y_hat[i,1];
  //   recovered[i]  = y_hat[i,11];
  }




for (i in 1: n_days){
    
    for (x in 1:4){
    
    if (x == 1){  
       R_0[i,x] = (beta[x] * (1 - (rho * pPCR_daily[i]  * phi_PCR))) / gamma  ;  
      } else {
        
      R_0[i,x] = (beta[x] * (1 - delta_C_daily[i])) / gamma    ;   
      
      }
      
    }
}

}
