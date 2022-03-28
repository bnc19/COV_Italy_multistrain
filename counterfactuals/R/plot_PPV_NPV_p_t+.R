
############################################################################################
# Function to calculate PPV/NPV, p(T+) based on test sensitivity and
# and specificity and  alternative antigen and molecular based testing strategies
############################################################################################

calc_test_met = function(C_D_prev, 
                         PCR_sens = 0.92 ,
                         Ag_sens_C = 0.643,
                         Ag_sens_D = 0.0,
                         Ag_spec = 0.99,
                         PCR_spec = 0.995,
                         proportion_Ag = 1, 
                         Legend = F,
                         y_axis_ticks = F,
                         x_axis_ticks = F, 
                         title) {
  
  
  proportion_PCR = (1-proportion_Ag)
  
  # calculate TP / FP / TN / FN
  
  TP_ag = (Ag_sens_C * C_D_prev$Var1) + ( Ag_sens_D * C_D_prev$Var2)
  FP_ag = (1 - Ag_spec) * (1-(C_D_prev$Var1 +  C_D_prev$Var2))
  TN_ag = Ag_spec * (1-  (C_D_prev$Var1 +  C_D_prev$Var2))
  FN_ag = (1 - Ag_sens_C) * C_D_prev$Var1 + (1 - Ag_sens_D) *  C_D_prev$Var2 
  
  
  TP_PCR = (PCR_sens * C_D_prev$Var1) + ( PCR_sens *  C_D_prev$Var2)
  FP_PCR = (1 - PCR_spec) * (1-(C_D_prev$Var1 +  C_D_prev$Var2))
  TN_PCR = PCR_spec * (1-(C_D_prev$Var1 +  C_D_prev$Var2))
  FN_PCR = (1 - PCR_sens) * C_D_prev$Var1 + (1 - PCR_sens) *  C_D_prev$Var2 
  
  # Calculate PPV / NPV p(T+)
  
  PT_pos =  (proportion_Ag * (FP_ag + TP_ag)) + (proportion_PCR * (FP_PCR + TP_PCR))
  PT_neg = (proportion_Ag * (FN_ag + TN_ag)) + (proportion_PCR * (FN_PCR + TN_PCR))
  
  PPV = ((proportion_Ag * TP_ag) + (proportion_PCR * TP_PCR))  / PT_pos
  NPV = ((proportion_Ag * TN_ag) + (proportion_PCR * TN_PCR))  / PT_neg
  
  

  
  # output results 
  out = data.frame(C_D_prev, PPV, NPV, PT_pos, proportion_PCR,proportion_Ag) * 100
  

  
  if ( Legend == F) {
    leg = c("none")
  } else {
    leg = c(0.8,0.7)
  }
  
  if (x_axis_ticks == F){
    axis.text.x = element_blank()
  } else (
    axis.text.x = NULL
  )
  
  if (y_axis_ticks == F){
    axis.text.y = element_blank()
  } else (
    axis.text.y = NULL
  )
  
  
  
  PPV_plot = ggplot(out, aes(Var1 ,Var2 , fill= PPV)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name= "PPV") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10))  + 
    ggtitle(title)+ 
    geom_contour(aes(z = PPV), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = PPV), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = PPV), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  NPV_plot = ggplot(out, aes(Var1 ,Var2 , fill= NPV)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name= "NPV") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=10),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = NPV), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = NPV), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = NPV), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  
  PT_plot = ggplot(out, aes(Var1 ,Var2 , fill= PT_pos)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name = "p(T+)") +
    labs(x = " ", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  +  
    geom_contour(aes(z = PT_pos), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = PT_pos), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = PT_pos), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  pAg_plot = ggplot(out, aes(Var1 ,Var2 , fill= proportion_Ag)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name = "ANCOV %") +
    labs(x = " ", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = proportion_Ag), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = proportion_Ag), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = proportion_Ag), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  

  plot_out = plot_grid(PPV_plot, NPV_plot, PT_plot, pAg_plot, ncol = 1)
  
  return(list(out, plot_out))
  
}

#########################################################################
# Function to calculate PPV/NPV, p(T+) based on test sensitivity and
# and specificity and proportion of antigen vs. molecular tests conducted
# assuming molecular tests follows positive antigen test 
#########################################################################

calc_test_met_pos_ag = function(C_D_prev, 
                         PCR_sens = 0.92 ,
                         Ag_sens_C = 0.643,
                         Ag_sens_D = 0.0,
                         Ag_spec = 0.99,
                         PCR_spec = 0.995,
                         proportion_Ag = 1, 
                         Legend = F,
                         y_axis_ticks = F,
                         x_axis_ticks = F, 
                         title = "", 
                         percent_follow_up = 1) {
  

  # calculate TP / FP / TN / FN for antigen test 
  
  TP_ag = (Ag_sens_C * C_D_prev$Var1) + ( Ag_sens_D * C_D_prev$Var2)
  FP_ag = (1 - Ag_spec) * (1-(C_D_prev$Var1 +  C_D_prev$Var2))
  TN_ag = Ag_spec * (1-  (C_D_prev$Var1 +  C_D_prev$Var2))
  FN_ag = (1 - Ag_sens_C) * C_D_prev$Var1 + (1 - Ag_sens_D) *  C_D_prev$Var2 
  
  ## Probability of a positive antigen test = pPCR if positive Ag tests are followed up with PCR
  
  PT_pos_Ag = TP_ag +FP_ag  
  
  # follow up X% of negative Ag with PCR
  
  PCR_follow = percent_follow_up * PT_pos_Ag
  
  
  ## positive tests are positive PCR tests or Ag tests not followed up 
  
  PT_pos = (TP_ag * PCR_sens * percent_follow_up) + ((1- percent_follow_up) * TP_ag) + ( FP_ag * (1 - PCR_spec) * percent_follow_up) + (FP_ag * (1 - percent_follow_up))
  
  PT_neg = FN_ag + (TP_ag * (1 - PCR_sens) * percent_follow_up) + TN_ag + (FP_ag * PCR_spec * percent_follow_up)
  
  TP = TP_ag * PCR_sens * percent_follow_up + (1- percent_follow_up) * TP_ag
  
  TN = TN_ag + (FP_ag * PCR_spec * percent_follow_up)

  
  # rescale to add to 100%
  
  pPCR = PCR_follow / (1 + PCR_follow)
  pAg = 1 / (1 + PCR_follow)
  
  
  ### PPV and NPV 
  
  PPV = TP / PT_pos
  NPV = TN / PT_neg
  
  # output results 
  out = data.frame(C_D_prev, PPV, NPV, PT_pos, pPCR, pAg) * 100
  
  # plot results 
  
  # percent_Ag = proportion_Ag * 100
  # percent_PCR= proportion_PCR * 100
  # title = paste(paste0(percent_Ag,"%"), "ANCOV,", paste0(percent_PCR,"%"), "DNCOV")
  
  
  
  if ( Legend == F) {
    leg = c("none")
  } else {
    leg = c(0.8,0.7)
  }
  
  if (x_axis_ticks == F){
    axis.text.x = element_blank()
  } else (
    axis.text.x = NULL
  )
  
  if (y_axis_ticks == F){
    axis.text.y = element_blank()
  } else (
    axis.text.y = NULL
  )
  
  
  PPV_plot = ggplot(out, aes(Var1 ,Var2 , fill= PPV)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name= "PPV") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10))  +   
    ggtitle(title)+ 
    geom_contour(aes(z = PPV), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = PPV), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = PPV), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  NPV_plot = ggplot(out, aes(Var1 ,Var2 , fill= NPV)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name= "NPV") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = NPV), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = NPV), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = NPV), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  
  PT_plot = ggplot(out, aes(Var1 ,Var2 , fill= PT_pos)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name = "p(T+)") +
    labs(x = " ", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg)+
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = PT_pos), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = PT_pos), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = PT_pos), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  pAg_plot = ggplot(out, aes(Var1 ,Var2 , fill= pAg)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name = "ANCOV %") +
    labs(x = " ", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = pAg), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = pAg), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = pAg), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  

  plot_out = plot_grid(PPV_plot, NPV_plot, PT_plot, pAg_plot, ncol = 1)
  
  return(list(out, plot_out))
  
}



#########################################################################
# Function to calculate PPV/NPV, p(T+) based on test sensitivity and
# and specificity and proportion of antigen vs. molecular tests conducted
# assuming molecular tests follows X% negative antigen test 
#########################################################################

calc_test_met_neg_ag = function(C_D_prev, 
                                PCR_sens = 0.92 ,
                                Ag_sens_C = 0.645,
                                Ag_sens_D = 0.0,
                                Ag_spec = 0.9968,
                                PCR_spec = 1,
                                proportion_Ag = 1, 
                                Legend = F,
                                y_axis_ticks = F,
                                x_axis_ticks = F, 
                                title = "",
                                percent_follow_up = 1) {
  
  
  
  # calculate TP / FP / TN / FN for antigen test 
  
  TP_ag = (Ag_sens_C * C_D_prev$Var1) + ( Ag_sens_D * C_D_prev$Var2)
  FP_ag = (1 - Ag_spec) * (1-(C_D_prev$Var1 +  C_D_prev$Var2))
  TN_ag = Ag_spec * (1-  (C_D_prev$Var1 +  C_D_prev$Var2))
  FN_ag = (1 - Ag_sens_C) * C_D_prev$Var1 + (1 - Ag_sens_D) *  C_D_prev$Var2 
  
  ## Probability of a negative antigen test = pPCR if negative Ag tests are followed up with PCR
  
  PT_neg_Ag = TN_ag + FN_ag  
  
  # follow up X% of negative Ag with PCR
  
  PCR_follow = percent_follow_up * PT_neg_Ag
  
  
  ## positive tests are positive PCR tests
  
  
  PT_pos = TP_ag + (FN_ag * percent_follow_up * PCR_sens) + FP_ag + (TN_ag * (1 - PCR_spec) * percent_follow_up)
  
  
  PT_neg = (FN_ag * (1 - percent_follow_up) ) + (FN_ag * (1 - PCR_sens) * percent_follow_up) +  (TN_ag * (1 - percent_follow_up)) + (TN_ag * PCR_spec* percent_follow_up)
  
  TP = TP_ag +( FN_ag * PCR_sens * percent_follow_up)
  
  TN =  (TN_ag * (1 - percent_follow_up)) +  (TN_ag * PCR_spec * percent_follow_up)
  
  

   
  # rescale to add to 100%
  
  pPCR = PCR_follow / (1 + PCR_follow)
  pAg = 1 / (1 + PCR_follow)
  
  
  
  # Calculate PPV / NPV p(T+)
  
  PPV = TP / PT_pos
  NPV = TN / PT_neg
  
  
  # output results 
  out = data.frame(C_D_prev, PPV, NPV, PT_pos, pPCR, pAg)* 100

  
# plot results  
  
  if ( Legend == F) {
    leg = c("none")
  } else {
    leg = c(0.8,0.7)
  }
  
  if (x_axis_ticks == F){
    axis.text.x = element_blank()
  } else (
    axis.text.x = NULL
  )
  
  if (y_axis_ticks == F){
    axis.text.y = element_blank()
  } else (
    axis.text.y = NULL
  )
  
  
  
  PPV_plot = ggplot(out, aes(Var1 ,Var2 , fill= PPV)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name= "PPV") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14 ),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10))  +   
    ggtitle(title)+ 
    geom_contour(aes(z = PPV), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = PPV), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = PPV), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  NPV_plot = ggplot(out, aes(Var1 ,Var2 , fill= NPV)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name= "NPV") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = NPV), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = NPV), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = NPV), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  
  PT_plot = ggplot(out, aes(Var1 ,Var2 , fill= PT_pos)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name = "p(T+)") +
    labs(x = " ", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  + 
    geom_contour(aes(z = PT_pos), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = PT_pos), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = PT_pos), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  pAg_plot = ggplot(out, aes(Var1 ,Var2 , fill= pAg)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdPu", limits = c(0,100), name = "ANCOV %") +
    labs(x = " ", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5),
      legend.position = leg) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))  +  
    geom_contour(aes(z = pAg), 
                 breaks = 50, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = pAg), 
                 breaks = 25, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = pAg), 
                 breaks = 75, col = 'black', size = .7, linetype = 3)
  
  
  
  
  
  
  plot_out = plot_grid(PPV_plot, NPV_plot, PT_plot, pAg_plot, ncol = 1)
  return(list(out, plot_out))
  
}












