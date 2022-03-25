################ Plot R0 as a function of proportion of DNCOV tests (Figure S7) ################ 


################    Set up ################   

# setwd("Q:/testing_SARS-CoV-2_variants")

# rm(list = ls())

library(tidyverse)
library(ggpubr)
library(cowplot)

################  Import data ###############



variant_data = read.csv("counterfactuals/Veneto/Veneto_variant_data.csv")[, -1]
variant_data$Date = as.Date.character(variant_data$Date, format = "%Y-%m-%d")

modelfit1_df_veneto = read.csv( "counterfactuals/Results/modelfit1_df_veneto.csv")[,-1]


modelfit1_df_italy = read.csv("counterfactuals/Results/modelfit1_df_italy.csv")


variant_data_italy = read.csv("counterfactuals/Italy/Italy_variant_data.csv")[,-1]
variant_data_italy$Date = as.Date.character(variant_data_italy$Date, format = "%Y-%m-%d")


################  Plot R0 for Veneto ################  

modelfit1_df_veneto_R0 = modelfit1_df_veneto %>%  
  select( c(M_R0_M,A_R0_M,O_R0_M,Al_R0_M,pPCR_M,M_R0_L,A_R0_L,O_R0_L,Al_R0_L,pPCR_L,M_R0_U,A_R0_U,O_R0_U,Al_R0_U,pPCR_U)) %>%  
  mutate(Date = variant_data$Date)



xmin = min(modelfit1_df_veneto_R0$Date)
xmax = max(modelfit1_df_veneto_R0$Date)

R0_Veneto_plot = ggplot(modelfit1_df_veneto_R0, aes(x = Date, y = M_R0_M    )) +
  geom_ribbon(aes(ymin = M_R0_L, ymax = M_R0_U),
              fill = "orange",
              alpha = 0.5) +
  geom_line(aes(y = M_R0_M    , color = "M234I-A376T R0"), size = 1) +
  geom_ribbon(aes(ymin = A_R0_L, ymax = A_R0_U),
              fill = "sky blue",
              alpha = 0.5) +
  geom_line(aes(y = A_R0_M  , color = "A220V R0"), size = 1) +
  geom_ribbon(aes(ymin = O_R0_L, ymax = O_R0_U),
              fill = "pink",
              alpha = 0.5) +
  geom_line(aes(y = O_R0_M, color = "Other R0"), size = 1) +
  geom_ribbon(aes(ymin = Al_R0_L, ymax = Al_R0_U),
              fill = "plum3",
              alpha = 0.5) +
  geom_line(aes(y = Al_R0_M, color = "Alpha R0"), size = 1) +
  geom_ribbon(aes(ymin = pPCR_L, ymax = pPCR_U),
              fill = "black",
              alpha = 0.2) +
  geom_line(aes(y = pPCR_M, color = paste("p(DN)")), size = 1) +
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T R0' = 'orange',
      'A220V R0' = 'sky blue',
      'Other R0' = 'pink',
      'Alpha R0' = 'plum3',
      "p(DN)" = "black")
  ) +
  labs(x = " ",
       y = paste0(" ")) +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    legend.position = c(0.2, 0.82),
    legend.key.height = unit(.5, 'cm')
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Veneto") + 
  geom_segment(x = xmin, y = 1, xend = xmax, yend = 1, color = "red", linetype = 3, size = 1) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) 


R0_Veneto_plot



################  Plot R0 for Italy ################  



modelfit1_df_italy_R0 = modelfit1_df_italy %>%  
  select( c(M_R0_M,A_R0_M,O_R0_M,Al_R0_M,pPCR_M,M_R0_L,A_R0_L,O_R0_L,Al_R0_L,pPCR_L,M_R0_U,A_R0_U,O_R0_U,Al_R0_U,pPCR_U)) %>%  
  mutate(Date = variant_data_italy$Date)



xmin_italy = min(modelfit1_df_italy_R0$Date)
xmax_italy = max(modelfit1_df_italy_R0$Date)


R0_Italy_plot = ggplot(modelfit1_df_italy_R0, aes(x = Date, y = M_R0_M    )) +
  geom_ribbon(aes(ymin = M_R0_L, ymax = M_R0_U),
              fill = "orange",
              alpha = 0.5) +
  geom_line(aes(y = M_R0_M    , color = "M234I-A376T R0"), size = 1) +
  geom_ribbon(aes(ymin = A_R0_L, ymax = A_R0_U),
              fill = "sky blue",
              alpha = 0.5) +
  geom_line(aes(y = A_R0_M  , color = "A220V R0"), size = 1) +
  geom_ribbon(aes(ymin = O_R0_L, ymax = O_R0_U),
              fill = "pink",
              alpha = 0.5) +
  geom_line(aes(y = O_R0_M, color = "Other R0"), size = 1) +
  geom_ribbon(aes(ymin = Al_R0_L, ymax = Al_R0_U),
              fill = "plum3",
              alpha = 0.5) +
  geom_line(aes(y = Al_R0_M, color = "Alpha R0"), size = 1) +
  geom_ribbon(aes(ymin = pPCR_L, ymax = pPCR_U),
              fill = "black",
              alpha = 0.2) +
  geom_line(aes(y = pPCR_M, color = "p(DN)"), size = 1) + 
  scale_colour_manual(
    name = '',
    values = c(
      'M234I-A376T R0' = 'orange',
      'A220V R0' = 'sky blue',
      'Other R0' = 'pink',
      'Alpha R0' = 'plum3',
      "p(DN)" = "black")
  ) +
  labs(x = " ",
       y = paste0(" ")) +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    legend.position = c("none"),
    legend.key.height = unit(.5, 'cm'),
    axis.text.y = element_blank()
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Rest of Italy") + 
  geom_segment(x = xmin_italy, y = 1, xend = xmax_italy, yend = 1, color = "red", linetype = 3, size = 1) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) 

R0_Italy_plot



################# save ################# 



R0_fig = plot_grid(R0_Veneto_plot, R0_Italy_plot, ncol = 2, labels = c("a", "b"))
R0_fig_anot = annotate_figure(R0_fig, bottom = "Date", element_text(size = 22))

ggsave(
  plot = R0_fig_anot,
  filename = "figS7.jpg",
  height = 20,
  width = 40,
  units = "cm",
  dpi = 300
)



