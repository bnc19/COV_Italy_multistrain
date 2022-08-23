################################################################################
# Script to run analysis and plot the performance of alternative antigen and   #
# molecular based testing strategies in diagnosing concordant and discordant   #
# variants (Figure 6).                                                         #
################################################################################


# Set up -----------------------------------------------------------------------

# setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")
# rm(list = ls())

# packages 
library(tidyverse)
library(ggpubr)
library(cowplot)

# functions and output
source("R/plot_PPV_NPV_p_t+.R")
dir.create("Figures")

# define prevalence ------------------------------------------------------------

prevalence_C = seq(0, 1, by = 0.005)
prevalence_D = seq(0, 1, by = 0.005)

C_D_prev = expand.grid(prevalence_C, prevalence_D)

C_D_prev$Var1  = ifelse(C_D_prev$Var1 + C_D_prev$Var2 > 1, NA, C_D_prev$Var1)

C_D_prev$Var2  = ifelse(C_D_prev$Var1 + C_D_prev$Var2 > 1, NA, C_D_prev$Var2)

C_D_prev = drop_na(C_D_prev)


# Results ----------------------------------------------------------------------



test_spec_pos_ag = calc_test_met_pos_ag(
  C_D_prev = C_D_prev,
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.689,
  Ag_sens_D = 0.0,
  Ag_spec = 0.985,
  PCR_spec = 1,
  proportion_Ag = 1,
  x_axis_ticks = T,
  y_axis_ticks = T,
  percent_follow_up = 1,
  title = "molecular follows 100% antigen +"
)

test_spec_pos_ag_50_fu = calc_test_met_pos_ag(
  C_D_prev = C_D_prev,
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.689,
  Ag_sens_D = 0.0,
  Ag_spec = 0.985,
  PCR_spec = 1,
  proportion_Ag = 1,
  x_axis_ticks = T,
  y_axis_ticks = T,
  percent_follow_up = .5,
  title = "molecular follows 50% antigen +"
)



test_spec_neg_ag = calc_test_met_neg_ag(
  C_D_prev = C_D_prev,
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.689,
  Ag_sens_D = 0.0,
  Ag_spec = 0.985,
  PCR_spec = 1,
  proportion_Ag = 1,
  x_axis_ticks = T,
  y_axis_ticks = T,
  title = "molecular follows 100% antigen -",
  percent_follow_up = 1
)



test_spec_neg_ag_50_fu = calc_test_met_neg_ag(
  C_D_prev = C_D_prev,
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.689,
  Ag_sens_D = 0.0,
  Ag_spec = 0.985,
  PCR_spec = 1,
  proportion_Ag = 1,
  x_axis_ticks = T,
  y_axis_ticks = T,
  title = "molecular follows 50% antigen -",
  percent_follow_up = 0.50
)



test_spec_100 = calc_test_met(
  C_D_prev = C_D_prev,
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.689,
  Ag_sens_D = 0.0,
  Ag_spec = 0.985,
  PCR_spec = 1,
  proportion_Ag = 1,
  x_axis_ticks = T,
  y_axis_ticks = T,
  Legend = T,
  title = "only antigen"
)

# prob testing pos if conc prev < 0.05
range(filter(test_spec_100[[1]], Var1  < 0.5)$PT_pos)
range(filter(test_spec_100[[1]], Var1  < 0.5)$PPV, na.rm = T)

range(filter(test_spec_100[[1]], Var1  < 0.5 &
               Var2 > 25)$NPV, na.rm = T)
range(filter(test_spec_100[[1]], Var1  < 0.5 &
               Var2 > 25)$PPV, na.rm = T)



test_spec_0 = calc_test_met(
  C_D_prev = C_D_prev,
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.689,
  Ag_sens_D = 0.0,
  Ag_spec = 0.985,
  PCR_spec = 1,
  proportion_Ag = 0,
  x_axis_ticks = T,
  title = "only molecular"
)



test_spec_test_prop = plot_grid(
  test_spec_100[[2]]  ,
  test_spec_pos_ag_50_fu[[2]],
  test_spec_pos_ag[[2]],
  test_spec_neg_ag_50_fu[[2]],
  test_spec_neg_ag[[2]],
  test_spec_0[[2]],
  ncol = 6,
  labels = c("a", "b", "c", "d", "e", "f"),
  label_size = 5
)


test_spec_test_prop_anot = annotate_figure(
  test_spec_test_prop,
  bottom = text_grob("Concordant variant prevalence", size = 7),
  left = text_grob(
    "Discordant variant prevalence",
    size = 7,
    rot = 90
  )
)



ggsave(
  plot = test_spec_test_prop_anot,
  filename = "Figures/Figure6.pdf",
  height = (22 * 4 / 6),
  width = 20.5,
  units = "cm",
  dpi = 800
)
