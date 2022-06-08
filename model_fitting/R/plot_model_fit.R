
plot_stan_fit = function(stan_fit,
                         file_path,
                         varaints = c("M234I-A376T", "A220V", "Other", "Alpha"),
                         # initial conditions
                         n_pop_it = 59257566 - 4847026,
                         n_recov_it = 1482377 - 93401,
                         n_pop_veneto= 4847026,
                         n_recov_veneto= 93401,
                         # index data Italy 
                         index_M_it = 6:14,
                         index_A_it = 3:14,
                         index_O_it =  3:9,
                         index_Al_it = 9:14,
                         # index data Veneto 
                         index_M_veneto = 8:14,
                         index_A_veneto = c(5,7:14),
                         index_O_veneto = c(5,7:9),
                         index_Al_veneto = 9:14,
                         # Dates Italy 
                         start_date_it =  "01-05-2020", 
                         end_date_it = "31-05-2021",
                         # Dates Veneto 
                         start_date_veneto=  "01-07-2020",
                         end_date_veneto= "31-05-2021"
){

  # Import Italy data ------------------------------------------------------------
  A_data_it = read.csv("data/Dataset_italy_A_v5.csv")$Freq 
  M_data_it  = read.csv("data/Dataset_italy_M_v5.csv")$Freq 
  O_data_it  = read.csv("data/Dataset_italy_O_v1.csv")$Freq 
  Al_data_it  = read.csv("data/Dataset_italy_Alpha_v1.csv")$Freq 
  n_seq_it  = read.csv("data/Dataset_italy_A_v5.csv")$TotSeq
  avg_daily_rep_inc_it  = read.csv("data/Dataset_italy_A_v5.csv")$new_reported_cases_daily
 
  # Import Veneto data -----------------------------------------------------------
  A_data_veneto = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq 
  M_data_veneto  = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq 
  O_data_veneto  = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq 
  Al_data_veneto  = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq 
  n_seq_veneto  = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq
  avg_daily_rep_inc_veneto  = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily
 
#  Italy------------------------------------------------------------------------
  start_date_d_it = as.Date.character(start_date_it, format = "%d-%m-%Y")
  end_date_d_it = as.Date.character(end_date_it, format = "%d-%m-%Y")
  
  all_dates_it = seq.Date(from = start_date_d_it, to = end_date_d_it, by = "days")
  
  n_days_it = length(all_dates_it)
  n_months_it = round(n_days_it / 30)
  
  n_seq_i_it = n_seq_it[(length(n_seq_it) - n_months_it + 1):length(n_seq_it)]
  avg_daily_rep_inc_i_it = 
    avg_daily_rep_inc_it[(length(avg_daily_rep_inc_it)-n_months_it+1): 
                           length(avg_daily_rep_inc_it)]
  
  # having extracted GISAID data from first date of fitting 
  # index variant specific data, removing NA values 
  
  index_M_i_it  = index_M_it - (length(M_data_it)- n_months_it)
  index_A_i_it  = index_A_it - (length(A_data_it)- n_months_it)
  index_O_i_it  = index_O_it - (length(O_data_it)- n_months_it)
  index_Al_i_it = index_Al_it - (length(Al_data_it)- n_months_it)
  
# Dates Veneto -----------------------------------------------------------------
  
  start_date_d_veneto = as.Date.character(start_date_veneto, format = "%d-%m-%Y")
  end_date_d_veneto = as.Date.character(end_date_veneto, format = "%d-%m-%Y")
  
  all_dates_veneto = 
    seq.Date(from = start_date_d_veneto, to = end_date_d_veneto,  by = "days")
 
  n_days_veneto = length(all_dates_veneto)
  n_months_veneto = round(n_days_veneto/ 30)
  
  
  n_seq_i_veneto = 
    n_seq_veneto[(length(n_seq_veneto) - n_months_veneto+ 1):length(n_seq_veneto)]
  
  avg_daily_rep_inc_i_veneto = 
    avg_daily_rep_inc_veneto[(length(avg_daily_rep_inc_veneto)-n_months_veneto+1): 
                               length(avg_daily_rep_inc_veneto)]
  
  # having extracted GISAID data from first date of fitting 
  # index variant specific data, removing NA values 
  
  index_M_i_veneto  = index_M_veneto- (length(M_data_veneto)- n_months_veneto)
  index_A_i_veneto  = index_A_veneto- (length(A_data_veneto)- n_months_veneto)
  index_O_i_veneto  = index_O_veneto- (length(O_data_veneto)- n_months_veneto)
  index_Al_i_veneto = index_Al_veneto- (length(Al_data_veneto)- n_months_veneto)
  
  
# extract posterior estimates --------------------------------------------------
SEIR__fit_ext = rstan::extract(fit)


# Plot model Italy -------------------------------------------------------------

dataConf_it = dataDay_it  = list()
# account for seeding and aggregating from day to month

monthly_date_it  = all_dates_it [ which(as.numeric(format(all_dates_it , "%d")) == 1) ]

dataList_it  = list(M_data_it[index_M_it ], A_data_it[index_A_it ], 
                    O_data_it[index_O_it ], Al_data_it[index_Al_it ])  # data

indexList_it  = list(index_M_i_it , index_A_i_it , index_O_i_it , index_Al_i_it )  # index data by month

# calculate binomial confidence intervals

for (i in 1:length(dataList_it )) {
  dataConf_it [[i]] = data.frame(
    Date = monthly_date_it [indexList_it [[i]]],
    binconf(dataList_it [[i]], n_seq_i_it[indexList_it [[i]]]) * 
      (avg_daily_rep_inc_i_it [indexList_it [[i]]] / (n_pop_it  -n_recov_it ) * 100000) ,
    variant = varaints[i]
  )
  
}

dataConfDf_it  = bind_rows(dataConf_it )
dataDay_it  = left_join(data.frame(Date = all_dates_it ),dataConfDf_it )


Italy_fit  = SEIR__fit_ext$incidence_it %>% as.data.frame.table() %>%
  rename(ni = iterations, time = Var2, variant = Var3, value = Freq) %>%
  dplyr::mutate(ni = as.numeric(ni),
                time = as.numeric(time)) %>%
  group_by(time,variant) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>% 
  ungroup( ) %>%  
  mutate(Date = rep(all_dates_it , each = 4)) %>% 
  left_join(dataDay_it , by = "Date") %>%
  mutate(Variant = ifelse(variant.x == "A", "M234I-A376T",
                          ifelse(variant.x == "B", "A220V",
                                 ifelse(variant.x=="C", "Other", "Alpha")))) %>% 
  ggplot(aes(x = Date , y = mean)) +
  geom_line(aes(color = Variant)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst ,color = variant.y)) +
  geom_errorbar(aes(ymin=Lower, ymax = Upper, color = variant.y)) +
  labs(y = paste0("")) +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "none") +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Italy")+
  ylim(c(0,150))





# Plot model Veneto -------------------------------------------------------------------
dataConf_veneto = dataDay_veneto  = list()

monthly_date_veneto = # account for seeding and aggregating from day to month
  all_dates_veneto [ which(as.numeric(format(all_dates_veneto, "%d")) == 1) ]

dataList_veneto = # data list 
  list(M_data_veneto[index_M_veneto], A_data_veneto[index_A_veneto], 
       O_data_veneto[index_O_veneto], Al_data_veneto[index_Al_veneto])  

indexList_veneto = # index data by month
  list(index_M_i_veneto , index_A_i_veneto, index_O_i_veneto, index_Al_i_veneto )  


# calculate binomial confidence intervals

for (i in 1:length(dataList_veneto )) {
  dataConf_veneto [[i]] = data.frame(
    Date = monthly_date_veneto [indexList_veneto [[i]]],
    binconf(dataList_veneto [[i]], n_seq_i_veneto[indexList_veneto [[i]]]) * 
      (avg_daily_rep_inc_i_veneto[indexList_veneto [[i]]] / (n_pop_veneto  -n_recov_veneto ) * 100000) ,
    variant = varaints[i]
  )
  
}

dataConfDf_veneto  = bind_rows(dataConf_veneto )
dataDay_veneto  = left_join(data.frame(Date = all_dates_veneto ),dataConfDf_veneto )


veneto_fit  = SEIR__fit_ext$incidence_ven %>% as.data.frame.table() %>%
  rename(ni = iterations, time = Var2, variant = Var3, value = Freq) %>%
  dplyr::mutate(ni = as.numeric(ni),
                time = as.numeric(time)) %>%
  group_by(time,variant) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>% 
  ungroup( ) %>%  
  mutate(Date = rep(all_dates_veneto, each = 4)) %>% 
  left_join(dataDay_veneto, by = "Date") %>%
  mutate(Variant = ifelse(variant.x == "A", "M234I-A376T",
                          ifelse(variant.x == "B", "A220V",
                                 ifelse(variant.x=="C", "Other", "Alpha")))) %>% 
  ggplot(aes(x = Date , y = mean)) +
  geom_line(aes(color = Variant)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst ,color = variant.y)) +
  geom_errorbar(aes(ymin=Lower, ymax = Upper, color = variant.y))+
  labs(y = paste0("Reported Incidence per 100,000 population")) +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c(0.12,0.86)) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Veneto") + 
  ylim(c(0,150))

#  Save Veneto and Italy fits --------------------------------------------------

rep_inc_plot = plot_grid(veneto_fit, Italy_fit)


ggsave(
  rep_inc_plot,
  file =  paste0(file_path,"reportedIncidencePlotp.jpg"),
  height = 20,
  width = 50,
  unit = "cm",
  dpi = 720
)

return(rep_inc_plot)
}