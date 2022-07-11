# Diagnosis of concordant and discordant variants ------------------------------
# under alternative testing strategies (figure 6) ------------------------------ 

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")
# Set up ----------------------------------------------------------------------- 
# rm(list = ls())

library(tidyverse)
library(ggpubr)
library(cowplot)


source("R/plot_PPV_NPV_p_t+.R")


# define prevalence ------------------------------------------------------------

prevalence_C = seq(0,1, by = 0.005)
prevalence_D = seq(0,1, by = 0.005)

C_D_prev = expand.grid(prevalence_C, prevalence_D)

C_D_prev$Var1  = ifelse(C_D_prev$Var1 + C_D_prev$Var2 > 1, NA, C_D_prev$Var1  )

C_D_prev$Var2  = ifelse(C_D_prev$Var1 + C_D_prev$Var2 > 1, NA, C_D_prev$Var2  )

C_D_prev = drop_na(C_D_prev)


# Results ----------------------------------------------------------------------



test_spec_pos_ag = calc_test_met_pos_ag(C_D_prev = C_D_prev, 
              PCR_sens = 0.92 ,
              Ag_sens_C = 0.689,
              Ag_sens_D = 0.0,
              Ag_spec = 0.985,
              PCR_spec = 1,
              proportion_Ag = 1,
              x_axis_ticks = T,
              y_axis_ticks = T,
              percent_follow_up = 1,
              title = "molecular follows 100% antigen +")

test_spec_pos_ag_50_fu = calc_test_met_pos_ag(C_D_prev = C_D_prev, 
                                        PCR_sens = 0.92 ,
                                        Ag_sens_C = 0.689,
                                        Ag_sens_D = 0.0,
                                        Ag_spec = 0.985,
                                        PCR_spec = 1,
                                        proportion_Ag = 1,
                                        x_axis_ticks = T,
                                        y_axis_ticks = T,
                                        percent_follow_up = .5,
                                        title = "molecular follows 50% antigen +")



test_spec_neg_ag = calc_test_met_neg_ag(C_D_prev = C_D_prev, 
                                        PCR_sens = 0.92 ,
                                        Ag_sens_C = 0.689,
                                        Ag_sens_D = 0.0,
                                        Ag_spec = 0.985,
                                        PCR_spec = 1,
                                           proportion_Ag = 1,
                                        x_axis_ticks = T,
                                        y_axis_ticks = T,
                                           title = "molecular follows 100% antigen -",
                                           percent_follow_up = 1)



test_spec_neg_ag_50_fu = calc_test_met_neg_ag(C_D_prev = C_D_prev, 
                                              PCR_sens = 0.92 ,
                                              Ag_sens_C = 0.689,
                                              Ag_sens_D = 0.0,
                                              Ag_spec = 0.985,
                                              PCR_spec = 1,
                                        proportion_Ag = 1,
                                        x_axis_ticks = T,
                                        y_axis_ticks = T,
                                        title = "molecular follows 50% antigen -",
                                        percent_follow_up = 0.50)



test_spec_100 = calc_test_met(C_D_prev = C_D_prev, 
                              PCR_sens = 0.92 ,
                              Ag_sens_C = 0.689,
                              Ag_sens_D = 0.0,
                              Ag_spec = 0.985,
                              PCR_spec = 1,
              proportion_Ag = 1,
              x_axis_ticks = T,
              y_axis_ticks = T,
              Legend = T,
              title = "only antigen")


test_spec_0 = calc_test_met(C_D_prev = C_D_prev, 
                            PCR_sens = 0.92 ,
                            Ag_sens_C = 0.689,
                            Ag_sens_D = 0.0,
                            Ag_spec = 0.985,
                            PCR_spec = 1,
                              proportion_Ag = 0,
                            x_axis_ticks = T,
                            title = "only molecular")



test_spec_test_prop = plot_grid(test_spec_100[[2]]  ,test_spec_pos_ag_50_fu[[2]], test_spec_pos_ag[[2]], 
                               test_spec_neg_ag_50_fu[[2]],test_spec_neg_ag[[2]], test_spec_0[[2]], ncol =6,
                                labels = c("a", "b", "c", "d", "e", "f"), label_size = 16)


test_spec_test_prop_anot = annotate_figure(test_spec_test_prop, bottom = text_grob("Concordant variant prevalence", size = 32),  
                                           left = text_grob("Discordant variant prevalence", size = 32, rot = 90) )



ggsave(
  plot = test_spec_test_prop_anot,
  filename = "fig6.jpg",
  height = (80 * 4/6),
  width = 80,
  units = "cm",
  dpi = 600
)


