
################ run analysis and plot genomic detection of concordant and         ################ 
################ discordant variants under alternative testing strategies  (fig 7) ################  


################   Set up  ################  


# rm(list = ls())


library(tidyverse)
library(cowplot)
library(ggpubr)


source("counterfactuals/R/plot_variant_detection_prob.R")


################  Produce plots ################

only_antigen  = calc_variant_detection_prob_AN_or_DN (
  sens_C = 0.643,
  sens_D = 0.0,
  Legend = T,
  y_axis_ticks = F,
  x_axis_ticks = F,
  y_axis_label =T,
  title = c("antigen only\nsequence antigen +",  "antigen only\nsequence antigen -", 
            "antigen only\nsequence 50% antigen -, 50% antigen +"))



only_molecular  = calc_variant_detection_prob_AN_or_DN (
  sens_C = 0.92,
  sens_D = 0.92,
  Legend = F,
  y_axis_label = T,
  y_axis_ticks = T,
  x_axis_ticks = F,
  title = c("molecular only\nsequence molecular +",  "molecular only\nsequence molecular -"))

follow_up50  =  calc_variant_detection_prob_AN_and_DN (
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.643,
  Ag_sens_D = 0.0,
  Ag_spec = 0.99,
  Legend = F,
  y_axis_ticks = F,
  x_axis_ticks = F,
  percent_follow_up = 0.5  ,
  title = c("molecular follows 50% antigen +\nsequence molecular +",  "molecular follows 50% antigen -\nsequence molecular +"))
  



follow_up100  =  calc_variant_detection_prob_AN_and_DN (
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.643,
  Ag_sens_D = 0.0,
  Ag_spec = 0.99,
  Legend = F,
  y_axis_ticks = F,
  y_axis_label =F,
  x_axis_ticks = F,
  percent_follow_up = .5,
  title = c("molecular follows 50% antigen +\nsequence molecular +",  "molecular follows 50% antigen -\nsequence molecular +"))

################  Save output ################

variant_detection_grid = plot_grid(only_antigen[[5]],only_antigen[[6]],  only_antigen[[7]],
                                   only_antigen[[2]],only_antigen[[3]],  only_antigen[[4]],
                                   only_molecular[[5]], follow_up100[[4]],follow_up100[[5]],
                                   only_molecular[[2]],follow_up100[[2]], follow_up100[[3]], 
                                   ncol =3,label_size = 16, align = "h",
                                   labels = c("a","b","c", "", "", "", "d","e","f"))

variant_detection_grid_anot = annotate_figure(variant_detection_grid,  
                                           bottom = text_grob("Specimens sequenced (%)", size = 32) )



ggsave(
  plot = variant_detection_grid_anot,
  filename = "fig7.jpg",
  height = 60,
  width = 60 * (3/4),
  units = "cm",
  dpi = 600
)

