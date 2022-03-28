########################################################################################################
# Function to format variant-specific genomic data to run multivariant models assuming fixed parameters.
########################################################################################################

Format_variant_data = function(index_M ,
                                       index_A ,
                                       index_O ,
                                       index_Al  ,
                                       start_date,
                                       end_date ,
                                       A_data ,
                                       M_data ,
                                       O_data ,
                                       Al_data, 
                                       average_daily_reported_incidence,
                                       n_seq) {
  
  
  
  ############ function ################
  
  
  library(tidyverse)

  # data wrangling to work with dates
  
  start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
  end_date_d = as.Date.character(end_date, format = "%d-%m-%Y")  
  
  daily_date = seq.Date(from = start_date_d, to = end_date_d ,  by = "days")
  
  
  n_days = length(daily_date)
  n_months = round(n_days / 30)
  
  
  ### having extracted GISAID data from first time of fitting (first occurrence of A variant), this allows to index variant specific data, removing NA values 
  
  index_M_i = index_M - (length(M_data)- n_months)
  index_A_i = index_A - (length(A_data)- n_months)
  index_O_i = index_O - (length(O_data)- n_months)
  index_Al_i  = index_Al - (length(Al_data)- n_months)
  
  
  
  data_list = list(M_data[index_M], A_data[index_A], O_data[index_O], Al_data[index_Al])  # data
  index_list = list(index_M_i, index_A_i, index_O_i, index_Al_i)  # index data by month
  

  

 
  # account for seeding and aggregating from day to month
  
  monthly_date = daily_date[ which(as.numeric(format(daily_date, "%d")) == 1) ]

  
  average_daily_reported_incidence_i = average_daily_reported_incidence[(length(average_daily_reported_incidence)-n_months+1): length(average_daily_reported_incidence)]
  
  
  n_seq_i = n_seq[(length(n_seq) - n_months + 1):length(n_seq)]
  
  

  
  dataConf = c()
  
  # calculate binomial confidence intervals
  
  for (i in 1:length(data_list)) {
    dataConf[[i]] = data.frame(
      Date = monthly_date[index_list[[i]]],
      binconf(data_list[[i]], n_seq_i[index_list[[i]]]) * average_daily_reported_incidence_i[index_list[[i]]]
    )
    
  }
  
  out = left_join(data.frame(Date = daily_date), dataConf[[1]])
  
  for (i in 2:length(dataConf)) {
    out = left_join(out, dataConf[[i]], by = "Date")
  }
  
  
 
  
  colnames(out) = c(
    "Date",
    "M_Variant_M",
    "M_Variant_L",
    "M_Variant_U",
    "A_Variant_M",
    "A_Variant_L",
    "A_Variant_U",
    "O_Variant_M",
    "O_Variant_L",
    "O_Variant_U",
    "Al_Variant_M",
    "Al_Variant_L",
    "Al_Variant_U"
  )
  

  
  return(out)
  
}
