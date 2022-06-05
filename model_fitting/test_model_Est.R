
rm(list = ls())

setwd("Q:/COV_Italy_multistrain/model_fitting")

A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq 
M_data = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq 
O_data = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq 
Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq 
n_seq = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq

n_difeq = 11
n_pop = 4847026
n_recov = 93401
scale_time_step = 1

modelPath = "model/new_neg_bin_sens_estimate.stan"
Location = "Veneto"
n_chains =3
n_warmups =500
n_iter = 1000
n_thin = 1

pars = c("lp__", "beta[1]", "beta[2]","beta[3]","beta[4]", "rho"  , "omega[1]"  , "omega[2]"
         , "I0[1]", "I0[2]", "I0[3]", "I0[4]")

filePath = "Results/Veneto" 
ini_1_SEIR = function(){
  list(  beta = replicate(4,runif(1,1,3)),
         I0 = replicate(4, runif(1, 1,100)),
         omega = replicate(2,runif(1,0.2,0.8)),
         rho = runif(1,0.2,0.8),
         k = runif(1,0.01,2)
  )}
average_daily_reported_incidence = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily
daily_reported_incidence = read.csv("data/dailyReportedIncidence_veneto.csv")$new_case 
daily_PCR = round(read.csv("data/Veneto_daily_test_data.csv")$pcr_daily)
daily_Ag= round(read.csv("data/Veneto_daily_test_data.csv")$antigen_daily)

average_daily_vaccination =  read.csv("data/data_vac_veneto_day.csv")$prop_vac
index_M = 8:14
index_A = c(5,7:14)
index_O = c(5,7:9)
index_Al = 9:14 
start_date = "01-07-2020"
end_date = "31-05-2021"
time_intervention = c("15-11-2020" , "15-03-2021")
time_seed_alpha = "01-11-2020"
time_seed_M = "01-10-2020"
sigma = 1 / 5.1
gamma = 1 / 2.1 
phi_PCR = 0.920
phi_Ag =  0.643 # change antigen sensitivity 
prior_seed_mean = 1
prior_seed_sd = 100
################################################################################



##########  If index not given, no missing data   ########## 

# index for GISAID data

if (is.null(index_A)) {
  index_M = 1:length(M_data)
}

if (is.null(index_A)) {
  index_A = 1:length(A_data)
}

if (is.null(index_O)) {
  index_O = 1:length(O_data)
}

if (is.null(index_Al)) {
  index_Al = 1:length(Al_data)
}


##########  data wrangling to work with dates   ########## 



start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
end_date_d = as.Date.character(end_date, format = "%d-%m-%Y")

daily_date = seq.Date(from = start_date_d, to = end_date_d ,  by = "days")

index_switch1 =  which(daily_date == as.Date.character(time_intervention[1], format = "%d-%m-%Y")) * scale_time_step
index_switch2 =  which(daily_date == as.Date.character(time_intervention[2], format = "%d-%m-%Y")) * scale_time_step
index_seed_alpha =  which(daily_date == as.Date.character(time_seed_alpha, format = "%d-%m-%Y")) * scale_time_step
index_seed_M =  which(daily_date == as.Date.character(time_seed_M, format = "%d-%m-%Y")) * scale_time_step


n_days = length(daily_date)
n_months = round(n_days / 30)

index_1st_month = c(which(as.numeric(format(daily_date, "%d")) == 1), (tail(n_days) + 1)) # index by 1st of each month and the final day of the last month



# format daily data to length of seed plus data fitting 

daily_PCR_i = daily_PCR[(length(daily_PCR) - n_days + 1):length(daily_PCR)] / scale_time_step
daily_PCR_i2 = rep(daily_PCR_i, each = scale_time_step)

daily_Ag_i = daily_Ag[(length(daily_Ag) - n_days + 1):length(daily_Ag)]  / scale_time_step
daily_Ag_i2 = rep(daily_Ag_i, each = scale_time_step)


## format vaccination 

time_vac =  which(daily_date == as.Date.character("26-12-2020", format = "%d-%m-%Y"))

average_daily_vaccination_i = rep(0, (time_vac+length(average_daily_vaccination)))

average_daily_vaccination_i[(time_vac +1):(time_vac  + length(average_daily_vaccination))] = average_daily_vaccination / scale_time_step

average_daily_vaccination_i2 = rep(average_daily_vaccination_i[1:n_days], each = scale_time_step) 



# # format monthly data to length of seed plus data fitting 

average_daily_reported_incidence_i = average_daily_reported_incidence[(length(average_daily_reported_incidence)-n_months+1): length(average_daily_reported_incidence)]


n_seq_i = n_seq[(length(n_seq) - n_months + 1):length(n_seq)]


### having extracted GISAID data from first time of fitting (first occurrence of A variant), this allows to index variant specific data, removing NA values 

index_M_i = index_M - (length(M_data)- n_months)
index_A_i = index_A - (length(A_data)- n_months)
index_O_i = index_O - (length(O_data)- n_months)
index_Al_i  = index_Al - (length(Al_data)- n_months)

################################################################################

 # data
 
 n_variants = 4
 S0 = n_pop-n_recov

 time_switch1 = index_switch1
 time_switch2 =index_switch2
 
 time_seed_M = index_seed_M
 time_seed_alpha = index_seed_alpha
 
 month_index = index_1st_month
 
 seed_mean = prior_seed_mean 
 seed_sd = prior_seed_sd

 n_ts = n_days * scale_time_step

 
 Ag_daily  = round(daily_Ag_i2)
 PCR_daily = round(daily_PCR_i2)
 vac = average_daily_vaccination_i2
 
 
 n_seq_M = n_seq[index_M]
 n_seq_A = n_seq[index_A]
 n_seq_O = n_seq[index_O]
 n_seq_Al = n_seq[index_Al]
 
 n_reported_M = average_daily_reported_incidence[index_M]
 n_reported_A = average_daily_reported_incidence[index_A]
 n_reported_O = average_daily_reported_incidence[index_O]
 n_reported_Al = average_daily_reported_incidence[index_Al]
 
 y_M = round(M_data[index_M])
 y_A = round(A_data[index_A])
 y_O = round(O_data[index_O])
 y_Al = round(Al_data[index_Al])
 
 
 n_data_A  =  length(index_A)
 n_data_M =  length(index_M)
 n_data_O =  length(index_O)
 n_data_Al =  length(index_Al)
## empty states 

 S = Q = R = rep(0, (n_ts+1))
 I = E = matrix(0, nrow = (n_ts+1), ncol = n_variants)
 FOI = incidence = delta= matrix(0, nrow = n_ts, ncol = n_variants)
 pPCR_daily = p_daily = rep(0, n_ts)
 
 
 beta2 = matrix(0, nrow = n_ts, ncol =n_variants)
 
 
 
 # parameters 
 posterior = read.csv("100_posterior_samples_veneto.csv")
 
 beta = as.numeric(posterior[3,c(7,6,8,9)])
 
 I0 =  as.numeric(posterior[3,14:17])
 
 rho = as.numeric(posterior[3,10])
 omega = as.numeric(posterior[3,11:12])

# Initial values 
 
 S[1] = S0 - I0[2] - I0[3] ;
 
 for(i in 1:n_variants) E[1,i] = 0 ;
 
 I[1,1] = 0 ;
 I[1,2] = I0[2] ;
 I[1,3] = I0[3] ;
 I[1,4] = 0 ;
 
 Q[1] = 0 ;
 R[1] = 0 ; 
 

 
 # model interventions 
 
 for(t in 1:n_ts) {
   for(i in 1: n_variants)  beta2[t,i] =  (beta[i] / scale_time_step);
 } 
 
 for(t in time_switch1 :  (time_switch2-1)) {
   for(i in 1: n_variants)  beta2[t,i] = (1-omega[1]) * (beta[i]/ scale_time_step);
 }  
 
 for(t in  time_switch2 :  n_ts) {
   for(i in 1: n_variants)  beta2[t,i] = (1-omega[2]) * (beta[i]/ scale_time_step);
 }  
 

 #--- Simulate ---#
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
   
   
   for(i in 1:n_variants) incidence[t,i] = (sigma*E[t,i] * delta[t,i]);
 }
 
 
 

 monthly_incidence  =matrix(0, nrow = n_months, ncol = n_variants)

 ## aggregate to months  
 index = 1; 
 
 for (m in 1:n_months){ 
   ind = index + 1;
   for (i in 1:n_variants){
     monthly_incidence[m,i] =  mean( incidence[month_index[index]:(month_index[ind]-1),i] ); 
   }
   index = index + 1;
 }
 
 
 
 # variant specific mean parameters for likelihood  
 lambda_M =lambda_A = lambda_O = lambda_Al = c()
 
 
 for (i in 1:n_data_M){  
   lambda_M[i] = n_reported_M[i] * (y_M[i] /  n_seq_M[i]);
 }
 
 for (i in 1:n_data_A){  
   lambda_A[i] =  n_reported_A[(i)] * (y_A[i] / n_seq_A[(i)]);
 }
 
 for (i in 1:n_data_O){  
   lambda_O[i] = n_reported_O[(i)] * (y_O[i] / n_seq_O[(i)]);
 }
 
 for (i in 1:n_data_Al){  
   lambda_Al[i] = n_reported_Al[(i)] *(y_Al[i] /  n_seq_Al[(i)]);
 }
 
 
 
 plot(monthly_incidence[,1],)

 plot(monthly_incidence[,2])

 plot(monthly_incidence[,3])

 plot(monthly_incidence[,4])

  
 