

# calculated binomial confint for variant specific data ------------------------
# plot the model reported incidence over time against binomial conf ------------
# returns a ggplot -------------------------------------------------------------

plot_model_fit = function(posts_df,
                          location,
                          start_date,
                          end_date,
                          variants = c(  "M234I-A376T", "A220V", "Other", "Alpha")
                          ) {
  
  
  
  
  library(Hmisc)
  
  if (location == "Rest of Italy") {
    # Import Italy data ----------------------------------------------------------
    A_data = read.csv("data/Dataset_italy_A_v5.csv")$Freq_new
    M_data  = read.csv("data/Dataset_italy_M_v5.csv")$Freq_new
    O_data  = read.csv("data/Dataset_italy_O_v1.csv")$Freq_new
    Al_data  = read.csv("data/Dataset_italy_Alpha_v1.csv")$Freq_new
    n_seq  = read.csv("data/Dataset_italy_A_v5.csv")$TotSeq_new
    avg_daily_rep_inc  = read.csv("data/Dataset_italy_A_v5.csv")$new_reported_cases_daily_new
    
    n_pop = 59257566 - 4847026
    n_recov  = 1482377 - 93401
    index_M  = 6:14
    index_A  = 3:14
    index_O  =  3:9
    index_Al  = 9:14
    
 
      leg = c("none")
       y = " "

    
    
  } else if (location == "Veneto"){
    # Import Veneto data ---------------------------------------------------------
  A_data = read.csv("data/Dataset_Veneto_A_v5.csv")$Freq
  M_data  = read.csv("data/Dataset_Veneto_M_v5.csv")$Freq
  O_data  = read.csv("data/Dataset_Veneto_O_v1.csv")$Freq
  Al_data = read.csv("data/Dataset_Veneto_Alpha_v1.csv")$Freq
  n_seq  = read.csv("data/Dataset_Veneto_A_v5.csv")$TotSeq
  avg_daily_rep_inc = read.csv("data/Dataset_Veneto_A_v5.csv")$new_reported_cases_daily
  
  n_pop = 4847026
  n_recov= 93401
  index_M = 8:14
  index_A = c(5,7:14)
  index_O = c(5,7:9)
  index_Al = 9:14
  
  
  leg = c(0.18,0.78)
  
  y = "Reported incidence per 100,000 population"
  }
  # Dates ------------------------------------------------------------------------
  
  start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
  end_date_d = as.Date.character(end_date, format = "%d-%m-%Y")
  
  all_dates = seq.Date(from = start_date_d, to = end_date_d, by = "days")
  
  n_days = length(all_dates)
  n_months = round(n_days / 30)
  
  n_seq_i = n_seq [(length(n_seq) - n_months  + 1):length(n_seq)]
  avg_daily_rep_inc_i  =
    avg_daily_rep_inc [(length(avg_daily_rep_inc) - n_months + 1):length(avg_daily_rep_inc)]
  
  # having extracted GISAID data from first date of fitting
  # index variant specific data, removing NA values
  
  index_M_i   = index_M  - (length(M_data) - n_months)
  index_A_i   = index_A  - (length(A_data) - n_months)
  index_O_i   = index_O  - (length(O_data) - n_months)
  index_Al_i  = index_Al  - (length(Al_data) - n_months)
  
  
  # Plot model -------------------------------------------------------------------
  
  data_conf  = data_day   = list()
  # account for seeding and aggregating from day to month
  
  monthly_date   = all_dates  [which(as.numeric(format(all_dates  , "%d")) == 1)]
  
  data_list   = list(M_data [index_M], A_data [index_A],
                    O_data [index_O], Al_data [index_Al])  # data
  
  index_list   = list(index_M_i  , index_A_i  , index_O_i  , index_Al_i)  # index data by month
  
  # calculate binomial confidence intervals
  
  for (i in 1:length(data_list)) {
    data_conf  [[i]] = data.frame(
      Date = monthly_date  [index_list  [[i]]],
      binconf(data_list  [[i]], n_seq_i [index_list  [[i]]]) *
        (avg_daily_rep_inc_i  [index_list  [[i]]] / (n_pop   - n_recov) * 100000) ,
      variant = variants[i]
    )
    
  }
  
  data_conf_df   = bind_rows(data_conf)
  data_day   = left_join(data.frame(Date = all_dates), data_conf_df)
  
  posts_df$Date = as.Date.character(posts_df$Date, format = "%Y-%m-%d")
  
  
  data_plot = data_day %>%
    left_join(posts_df  , by = "Date") %>%
    mutate(Variant = ifelse(
      variant.y == "M",
      "M234I-A376T",
      ifelse(
        variant.y == "A",
        "A220V",
        ifelse(variant.y == "O", "Other",
        ifelse(variant.y == "Al","Alpha", "X"
      ))
    ))) %>%
    mutate(Data = ifelse(!is.na (PointEst), "Data", NA )) %>%  
    mutate(Model = ifelse(is.na (PointEst), "Model", NA )) 
  
  
  plot =  ggplot(data_plot, aes(x = Date , y = mean_rep)) +
    geom_line(aes(color = Variant, linetype = Model)) +
    geom_ribbon(aes(
      ymin = lower_rep,
      ymax = upper_rep,
      fill = Variant
    ), alpha = 0.4) +
    geom_point(aes(y = PointEst , color = variant.x, shape = Data),size=1) +
    geom_errorbar(aes(
      ymin = Lower,
      ymax = Upper,
      color = variant.x
    ),  width = 2.5) +
    labs(y = paste0(y)) +
    theme_bw() +
    ggtitle(paste(location)) +
    ylim(c(0, 105)) + 
    scale_shape_manual(values = c('Data' = 16)) +
    scale_linetype_manual(values = c("Model" = "solid")) +
    theme(
      text = element_text(size = 7),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      legend.position = leg,
      legend.title=element_text(size=5), 
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0, "mm"),
      legend.text=element_text(size=5) ,
      legend.spacing.y = unit(0, "mm"), 
      axis.text.x = element_text(size=5),
      legend.key.size = unit(.5,"line")) +
    scale_x_date(date_labels = "%b/%Y", breaks = "3 months")
  
  
  
  #  return plot -----------------------------------------------------------------
  
  
  
  return(list(plot, data_plot))
}
