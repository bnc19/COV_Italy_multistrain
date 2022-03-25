
################# Plot transmission dynamics of all variants and counterfactuals (Figure 5) ################   

################    Set up ################   

# setwd("Q:/testing_SARS-CoV-2_variants")

# rm(list = ls())

library("wesanderson")
library(cowplot) 
library(tidyverse)

################  Import data and fit for Veneto ################  

n_pop = 4847026
n_recov = 93401

variant_data = read.csv("counterfactuals/Veneto/Veneto_variant_data.csv")[, -1]
variant_data$Date = as.Date.character(variant_data$Date, format = "%Y-%m-%d")

modelfit1_df_veneto = read.csv( "counterfactuals/Results/modelfit1_df_veneto.csv")[,-1]

modelfit1_df_veneto_inc = modelfit1_df_veneto %>%  
select(! c(M_R0_M,A_R0_M,O_R0_M,Al_R0_M,pPCR_M,M_R0_L,A_R0_L,O_R0_L,Al_R0_L,pPCR_L,M_R0_U,A_R0_U,O_R0_U,Al_R0_U,pPCR_U)) 



modelfit1_df_veneto_inc = modelfit1_df_veneto_inc / (n_pop - n_recov) * 100000


variant_data[, 2:13] = variant_data[, 2:13] / (n_pop - n_recov) * 100000

modelfit1_vd_veneto = cbind(variant_data, modelfit1_df_veneto_inc)


################  Plot fit for Veneto ################  

modelfit1_plot = ggplot(modelfit1_vd_veneto, aes(x = Date, y = M_Variant_M)) +

  geom_ribbon(aes(ymin = M_fit_L, ymax = M_fit_U),
              fill = "orange",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = A_fit_L, ymax = A_fit_U),
              fill = "sky blue",
              alpha = 0.3) +
  geom_ribbon(aes(ymin = O_fit_L, ymax = O_fit_U),
              fill = "pink",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = Al_fit_L, ymax = Al_fit_U),
              fill = "plum3",
              alpha = 0.5, 
              linetype = 2) +
  geom_line(aes(y = M_fit_M, color = "M234I-A376T model"), size = 1.5) +
  geom_line(aes(y = A_fit_M, color = "A220V model"), size = 1.5) +
  geom_line(aes(y = O_fit_M, color = "Other model"), size = 1.5) +
  geom_line(aes(y = Al_fit_M, color = "Alpha model"), size = 1.5) +

  geom_point(shape = 19, size = 3, (aes(y = A_Variant_M, color = "A220V data"))) +
  geom_errorbar(color = "blue", aes(ymin = A_Variant_L, ymax = A_Variant_U), size =1) +
 
  
  geom_point(shape = 19, size = 3, (aes(y = O_Variant_M, color = "Other data"))) +
  geom_errorbar(color = "red", aes(ymin = O_Variant_L, ymax = O_Variant_U), size =1) +

 
  geom_point(shape = 19, size = 3, (aes(y = Al_Variant_M, color = "Alpha data"))) +
  geom_errorbar(color = "mediumpurple4", aes(ymin = Al_Variant_L, ymax = Al_Variant_U), size =1) +

  geom_point(shape = 19, size = 3, (aes(color = "M234I-A376T data"))) +
  geom_errorbar(aes(ymin = M_Variant_L, ymax = M_Variant_U), size =1) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T data' = 'black',
      'M234I-A376T model' = 'orange',
      "A220V data" = "blue",
      'A220V model' = 'sky blue',
      'Other model' = 'pink',
      'Other data' = 'red',
      'Alpha model' = 'plum3',
      'Alpha data' = 'mediumpurple4'
    )
  ) +
  labs(x = " ",
       y = paste0("Reported incidence per 100,000 population")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    legend.position = c(0.2, 0.83),
    legend.key.height = unit(.5, 'cm'),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14),
    axis.text.x = element_text(face="bold")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  ggtitle("Veneto") +
  scale_y_continuous(limits = c(0, 143), breaks = seq(0, 140, 20))

modelfit1_plot



# modelfit1_plot_true = ggplot(modelfit1_vd_veneto, aes(x = Date, y = T_M_fit_M)) +
#   geom_ribbon(aes(ymin = T_M_fit_L, ymax = T_M_fit_U),
#               fill = "orange",
#               alpha = 0.5) +
#   geom_line(aes(y = T_M_fit_M, color = "M234I-A376T model"), size = 1) +
#   geom_ribbon(aes(ymin = T_A_fit_L, ymax = T_A_fit_U),
#               fill = "sky blue",
#               alpha = 0.5) +
#   geom_line(aes(y = T_A_fit_M, color = "A220V model"), size = 1) +
#   geom_ribbon(aes(ymin = T_O_fit_L, ymax = T_O_fit_U),
#               fill = "pink",
#               alpha = 0.5) +
#   geom_line(aes(y = T_O_fit_M, color = "Other model"), size = 1) +
#   
#   geom_ribbon(aes(ymin = T_Al_fit_L, ymax = T_Al_fit_U),
#               fill = "plum3",
#               alpha = 0.5) +
#   geom_line(aes(y = T_Al_fit_M, color = "Alpha model"), size = 1) +
#   scale_colour_manual(
#     name = '',
#     values = c(
#       'M234I-A376T' = 'black',
#       'M234I-A376T model' = 'orange',
#       "A220V" = "blue",
#       'A220V model' = 'sky blue',
#       'Other model' = 'pink',
#       'Other' = 'red',
#       'Alpha model' = 'plum3',
#       'Alpha' = 'mediumpurple4' )) +
#   labs(x = "Date",
#        y = paste0("True incidence  \n per 100,000 population")) +
#   theme_bw() + theme(
#     text = element_text(size = 14),
#     axis.title.y = element_text(angle = 90, vjust = 0.7),
#     plot.title = element_text(hjust = 0.4),
#     legend.position = "none"
#   ) +
#   scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
#   ggtitle(" ") +
#   scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 50))
# 
# modelfit1_plot_true


################  Import data and fit for Italy ################  

n_pop_italy = (59257566 - 4847026 )
n_recov_italy = (1482377 - 93401)

modelfit1_df_italy = read.csv("counterfactuals/Results/modelfit1_df_italy.csv")


variant_data_italy = read.csv("counterfactuals/Italy/Italy_variant_data.csv")[,-1]
variant_data_italy$Date = as.Date.character(variant_data_italy$Date, format = "%Y-%m-%d")


modelfit1_df_italy_inc = modelfit1_df_italy %>%  
  select(! c(M_R0_M,A_R0_M,O_R0_M,Al_R0_M,pPCR_M,M_R0_L,A_R0_L,O_R0_L,Al_R0_L,pPCR_L,M_R0_U,A_R0_U,O_R0_U,Al_R0_U,pPCR_U)) 



modelfit1_df_italy_inc = modelfit1_df_italy_inc / (n_pop_italy - n_recov_italy) * 100000


variant_data_italy[, 2:13] = variant_data_italy[, 2:13] / (n_pop_italy - n_recov_italy) * 100000


modelfit1_vd = cbind(variant_data_italy, modelfit1_df_italy_inc)

################  Plot fit for Italy ################  

modelfit1_plot_italy = ggplot(modelfit1_vd, aes(x = Date, y = M_Variant_M)) +
  geom_ribbon(aes(ymin = M_fit_L, ymax = M_fit_U),
              fill = "orange",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = A_fit_L, ymax = A_fit_U),
              fill = "sky blue",
              alpha = 0.3) +
  geom_ribbon(aes(ymin = O_fit_L, ymax = O_fit_U),
              fill = "pink",
              alpha = 0.5) +
  geom_ribbon(aes(ymin = Al_fit_L, ymax = Al_fit_U),
              fill = "plum3",
              alpha = 0.5, 
              linetype = 2) +
  geom_line(aes(y = M_fit_M, color = "M234I-A376T model"), size = 1.5) +
  geom_line(aes(y = A_fit_M, color = "A220V model"), size = 1.5) +
  geom_line(aes(y = O_fit_M, color = "Other model"), size = 1.5) +
  geom_line(aes(y = Al_fit_M, color = "Alpha model"), size = 1.5) +
  
  geom_point(shape = 19, size = 3, (aes(y = A_Variant_M, color = "A220V data"))) +
  geom_errorbar(color = "blue", aes(ymin = A_Variant_L, ymax = A_Variant_U), size =1) +
  
  
  geom_point(shape = 19, size = 3, (aes(y = O_Variant_M, color = "Other data"))) +
  geom_errorbar(color = "red", aes(ymin = O_Variant_L, ymax = O_Variant_U), size =1) +
  
  
  geom_point(shape = 19, size = 3, (aes(y = Al_Variant_M, color = "Alpha data"))) +
  geom_errorbar(color = "mediumpurple4", aes(ymin = Al_Variant_L, ymax = Al_Variant_U), size =1) +
  
  geom_point(shape = 19, size = 3, (aes(color = "M234I-A376T data"))) +
  geom_errorbar(aes(ymin = M_Variant_L, ymax = M_Variant_U), size =1) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T data' = 'black',
      'M234I-A376T model' = 'orange',
      "A220V data" = "blue",
      'A220V model' = 'sky blue',
      'Other model' = 'pink',
      'Other data' = 'red',
      'Alpha model' = 'plum3',
      'Alpha data' = 'mediumpurple4'
    )
  ) +
  labs(x = "", y = paste0(" ")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    axis.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face="bold")
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  ggtitle("Rest of Italy") +
  scale_y_continuous(limits = c(0, 143), breaks = seq(0, 140, 20))

modelfit1_plot_italy

 
# modelfit1_plot_true_italy = ggplot(modelfit1_vd, aes(x = Date, y = T_M_fit_M)) +
#   geom_ribbon(aes(ymin = T_M_fit_L, ymax = T_M_fit_U),
#               fill = "orange",
#               alpha = 0.5) +
#   geom_line(aes(y = T_M_fit_M, color = "M234I-A376T model"), size = 1) +
#   geom_ribbon(aes(ymin = T_A_fit_L, ymax = T_A_fit_U),
#               fill = "sky blue",
#               alpha = 0.5) +
#   geom_line(aes(y = T_A_fit_M, color = "A220V model"), size = 1) +
#   
#   geom_ribbon(aes(ymin = T_O_fit_L, ymax = T_O_fit_U),
#               fill = "pink",
#               alpha = 0.5) +
#   geom_line(aes(y = T_O_fit_M, color = "Other model"), size = 1) +
#   
#   geom_ribbon(aes(ymin = T_Al_fit_L, ymax = T_Al_fit_U),
#               fill = "lavender",
#               alpha = 0.5) +
#   geom_line(aes(y = T_Al_fit_M, color = "Alpha model"), size = 1) +
#   scale_colour_manual(
#     name = '',
#     values = c(
#       'M234I-A376T' = 'black',
#       'M234I-A376T model' = 'orange',
#       "A220V" = "blue",
#       'A220V model' = 'sky blue',
#       'Other model' = 'pink',
#       'Other' = 'red',
#       'Alpha model' = 'lavender',
#       'Alpha' = 'mediumpurple4'
#     )
#   ) +
#   labs(x = "Date", y = paste0(" ")) +
#   theme_bw() + theme(
#     text = element_text(size = 14),
#     axis.title.y = element_text(angle = 90, vjust = 0.7),
#     axis.text.y = element_blank(),
#     plot.title = element_text(hjust = 0.4),
#     legend.position = "none"
#   ) +
#   scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
#   ggtitle(" ") +
#   scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 50))
# 
# modelfit1_plot_true_italy





################  Import cumulative data Veneto ################  

file_names = dir("counterfactuals/Results") #where you have your files

list_cumdf = lapply(paste0("counterfactuals/Results/",file_names),read.csv, )
names(list_cumdf) = file_names


################  Plot proportion of incidence reported in Veneto and Italy ################  

 proportion_reported_plot = list_cumdf$cuminc1_veneto.csv %>%  
  bind_rows(list_cumdf$cuminc1_df_italy.csv) %>%  
  select(X, Rep_M,     Rep_A,     Rep_O ,   Rep_Al  , Rep_tot) %>% 
  rename(Stat = X) %>% 
  mutate(Country = factor(rep(c("Veneto", "Rest of Italy"), each = 3), levels = c("Veneto", "Rest of Italy"))) %>%  
  pivot_longer(cols= -c(Country,Stat) ) %>%  
  mutate(Variant = factor (rep ( c( "M234I-A376T", "A220V", "Other", "Alpha", "Total model"), 6)),
         value = value * 100) %>% 
  pivot_wider(id_cols = c(Country,Variant), names_from = Stat,values_from = value)  %>%  
  ggplot(aes(x = Country, y = mean)) +
  geom_point( shape = 19, size = 4,aes(color = rep(c("M234I-A376T", "A220V", "Other", "Alpha", "Total model"),2)), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = rep(c("M234I-A376T", "A220V", "Other", "Alpha", "Total model"),2)),
    position = position_dodge(width = 0.5), width =  0.4, size = 1) +
  labs(x = " ", y = paste0("Percentage of cumulative incidence reported")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c(0.12,0.92),
    axis.text.x = element_text(face="bold")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T' = 'orange',
      'A220V' = 'sky blue',
      'Other' = 'pink',
      'Alpha' = 'plum3',
      "Total model" = " darkgrey"),
    breaks = c('Total model'))

################  Plot cumulative incidence by test strategy ################  

 cum_inc_test_strat_plot = list_cumdf$cuminc1_veneto.csv %>%  
  bind_rows(list_cumdf$cuminc5_df_veneto.csv,list_cumdf$cuminc3_df_veneto.csv,list_cumdf$cuminc3_df_sens_veneto.csv,  list_cumdf$cuminc7_df_veneto.csv, ) %>%  
  select(X,  true_SEIR_M , true_SEIR_A,  true_SEIR_O, true_SEIR_Al, total_incidence) %>% 
  rename(Stat = X) %>% 
  mutate(model = factor(rep(c("DNCOV follows \nANCOV+ \n(baseline)","DNCOV follows \nANCOV-", 
                       "only ANCOV \nmass testing","only ANCOV \nmass testing \n 87.5% test sensitivity", "Italy testing"), each = 3), 
                       levels = c("DNCOV follows \nANCOV+ \n(baseline)","DNCOV follows \nANCOV-", "Italy testing",
                                  "only ANCOV \nmass testing","only ANCOV \nmass testing \n 87.5% test sensitivity") )) %>% 
  pivot_longer(cols= -c(model,Stat) ) %>%  
  mutate(Variant = factor (rep ( c( "M234I-A376T", "A220V", "Other", "Alpha", "Total"), 15) , 
                           levels =c( "M234I-A376T", "A220V", "Other", "Alpha", "Total")),
         value = value / (n_pop - n_recov) * 100) %>% 
  pivot_wider(id_cols = c(model,Variant), names_from = Stat,values_from = value)  %>%  
  ggplot(aes(x = model, y = mean)) +
  geom_point( shape = 19, size = 4,aes(color = rep(c("M234I-A376T", "A220V", "Other", "Alpha", "Total"),5)), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = rep(c("M234I-A376T", "A220V", "Other", "Alpha", "Total"),5)),
                position = position_dodge(width = 0.5), width =  0.4, size = 1) +
  labs(x = " ", y = paste0("Cumulative incidence (%)")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c("none"),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14) ,
    axis.text.x = element_text(face="bold"))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T' = 'orange',
      'A220V' = 'sky blue',
      'Other' = 'pink',
      'Alpha' = 'plum3',
      "Total" = "darkgrey"))


################  Plot cumulative incidence by R0m ################  

cum_inc_r0m_plot = list_cumdf$cuminc2_veneto.csv %>%  
  bind_rows(list_cumdf$cuminc4_df_veneto.csv,list_cumdf$cuminc6_df_veneto.csv, ) %>%
  filter(X == "T_M_fit_M"|X == "T_M_fit_L"| X =="T_M_fit_U"|X == "total_M"| X =="total_L"|X =="total_U") %>% 
  mutate(Stat = rep(rep(c("mean", "lower", "upper"), each = 2), 3)) %>% 
  pivot_longer(cols= -c(X,Stat)) %>%   
  mutate(model = factor(rep(c("DNCOV follows \nANCOV+ \n(baseline)",
                              "only ANCOV \nmass testing", "Italy testing"), each = 30), 
                        levels = c("DNCOV follows \nANCOV+ \n(baseline)","Italy testing",
                                   "only ANCOV \nmass testing") )) %>% 
  mutate(Variant = rep(rep ( c( "M234I-A376T", "Total"), each = 5), 9),
         value = value / (n_pop - n_recov) * 100,
         `R0m scaling factor` = rep(c("0.8", "1", "1.2", "1.4", "1.6"), 18)) %>% 
  pivot_wider(id_cols = c(model,Variant,  `R0m scaling factor`), names_from = Stat,values_from = value)  %>%  
  ggplot(aes(x = model, y = mean)) +
  geom_point( shape = 19, size = 4,aes(group = interaction(Variant,`R0m scaling factor`),
                                       color = `R0m scaling factor`, shape = Variant), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper,group = interaction(Variant,`R0m scaling factor`),
                    color = `R0m scaling factor`, linetype = Variant),
                position = position_dodge(width = 0.5), width = 0.4, size = 1) +
  labs(x = " ", y = paste0("Cumulative incidence (%)")) +
  theme_bw() + theme(
    text = element_text(size = 16),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c(0.05,0.64),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14),
    axis.text.x = element_text(face="bold")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(values = wes_palette(name = "Darjeeling1")) 




################# Save ################# 

incidence_fig = plot_grid(
  modelfit1_plot,
  modelfit1_plot_italy,
  proportion_reported_plot,
  ncol = 3,
  labels = c("a", "b", "c"),
  align = "h"
)

fig5 = plot_grid(
  incidence_fig,
  cum_inc_test_strat_plot,
  cum_inc_r0m_plot,
  ncol = 1,
  labels = c("", "d", "e"),
  rel_heights = c(1 /3, 1 / 4, 1 / 4),
  align = "v,h"
)


ggsave(
  plot = fig5,
  filename = "fig5.jpg",
  height = 45,
  width = 52,
  units = "cm",
  dpi = 300
)

