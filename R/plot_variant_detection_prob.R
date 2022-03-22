#########################################################################
# Function to calculate prob variant detection if PCR / Ag
#########################################################################

calc_variant_detection_prob_AN_or_DN = function( 
                         sens_C ,
                         sens_D ,
                         Legend = F,
                         y_axis_ticks = F,
                         x_axis_ticks = F,
                         y_axis_label = F,
                         title
                         ) {
  
  
prevelance = seq(0,1, by = 0.005)


seq_perc = seq(0,0.03, by = 0.0001)

prev_seq_perc = expand.grid(prevelance, seq_perc)



# sequence positive samples 


TP_C = (sens_C * prev_seq_perc$Var1) 
TP_D = (sens_D *prev_seq_perc$Var1)


# sequence negative samples 

FN_C = (1 - sens_C) * prev_seq_perc$Var1 
FN_D = (1 - sens_D) * prev_seq_perc$Var1


Prob_genome_detection_C_pos =  TP_C * prev_seq_perc$Var2 
Prob_genome_detection_D_pos =  TP_D * prev_seq_perc$Var2 
Prob_genome_detection_C_neg =  FN_C * prev_seq_perc$Var2  
Prob_genome_detection_D_neg =  FN_D * prev_seq_perc$Var2 

Prob_genome_detection_C_both =  (0.5 * FN_C * prev_seq_perc$Var2) + (0.5 *  TP_C * prev_seq_perc$Var2 )
Prob_genome_detection_D_both = (0.5 * FN_D * prev_seq_perc$Var2) + (0.5 *  TP_D * prev_seq_perc$Var2)


out = data.frame(prev_seq_perc,Prob_genome_detection_C_pos, Prob_genome_detection_D_pos, 
                 Prob_genome_detection_C_neg, Prob_genome_detection_D_neg, 
                 Prob_genome_detection_C_both,   Prob_genome_detection_D_both) * 100

if ( Legend == F) {
  leg = c("none")
} else {
  leg = c(0.25,0.62)
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

if (y_axis_label == F){
  y_label = c(" ", "")
} else (
  y_label = c("Concordant variant prevalence (%)", "Discordant varaint prevalence (%)")
)

seq_con_plot_pos = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_C_pos)) + 
  geom_tile() +
scale_fill_distiller(palette = "RdYlBu" , limits = c(0,3), name= "Cases\ndetected (%)") +
  labs(x = "", y = y_label[1]) +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.5, size=14),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
#  scale_x_continuous(limits = c(0.00, 0.05), breaks = seq(0.00, 0.05, 0.005))  +
  geom_contour(aes(z = Prob_genome_detection_C_pos), 
               breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
  geom_contour(aes(z = Prob_genome_detection_C_pos), 
               breaks = 1, col = 'black', size = .7, linetype = 1)+
  geom_contour(aes(z = Prob_genome_detection_C_pos), 
               breaks = 2, col = 'black', size = .7, linetype = 3)



seq_dis_plot_pos = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_D_pos)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", limits = c(0,3), name= "Cases\ndetected (%)") +
  labs(x = "", y = y_label[2]) +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.5, size=14),
    legend.position = leg
  ) +
  ggtitle(title[1]) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  # scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
  geom_contour(aes(z = Prob_genome_detection_D_pos), 
               breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
  geom_contour(aes(z = Prob_genome_detection_D_pos), 
               breaks = 1, col = 'black', size = .7, linetype = 1)+
  geom_contour(aes(z = Prob_genome_detection_D_pos), 
               breaks =2, col = 'black', size = .7, linetype = 3)




seq_con_plot_neg = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_C_neg)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu" , limits = c(0,3), name= "Cases\ndetected (%)") +
  labs(x = "", y = "") +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.5, size=14),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  # scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
  geom_contour(aes(z = Prob_genome_detection_C_neg), 
               breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
  geom_contour(aes(z = Prob_genome_detection_C_neg), 
               breaks = 1, col = 'black', size = .7, linetype = 1)+
  geom_contour(aes(z = Prob_genome_detection_C_neg), 
               breaks = 2, col = 'black', size = .7, linetype = 3)



seq_dis_plot_neg = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_D_neg)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", limits = c(0,3), name= "Cases\ndetected (%)") +
  labs(x = "", y = "") +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.5, size=14),
    legend.position =  "none"
  ) +
  ggtitle(title[2])+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  #  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
  geom_contour(aes(z = Prob_genome_detection_D_neg), 
               breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
  geom_contour(aes(z = Prob_genome_detection_D_neg), 
               breaks = 1, col = 'black', size = .7, linetype = 1)+
  geom_contour(aes(z = Prob_genome_detection_D_neg), 
               breaks = 2, col = 'black', size = .7, linetype = 3)



### sequence both pos and neg 

seq_con_plot_both = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_C_both)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu" , limits = c(0,3), name= "Cases\ndetected (%)") +
  labs(x = "", y = "") +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.5, size=14),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  # scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
  geom_contour(aes(z = Prob_genome_detection_C_both), 
               breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
  geom_contour(aes(z = Prob_genome_detection_C_both), 
               breaks = 1, col = 'black', size = .7, linetype = 1)+
  geom_contour(aes(z = Prob_genome_detection_C_both), 
               breaks = 2, col = 'black', size = .7, linetype = 3)



seq_dis_plot_both = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_D_both)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", limits = c(0,3), name= "Cases\ndetected (%)n") +
  labs(x = "", y = "") +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.5, size=14),
    # axis.text.x = element_blank(), 
    # axis.text.y = axis.text.y, 
    legend.position =  "none"
  ) +
  ggtitle(title[3]) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  #  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
  geom_contour(aes(z = Prob_genome_detection_D_both), 
               breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
  geom_contour(aes(z = Prob_genome_detection_D_both), 
               breaks = 1, col = 'black', size = .7, linetype = 1)+
  geom_contour(aes(z = Prob_genome_detection_D_both), 
               breaks = 2, col = 'black', size = .7, linetype = 3)




return(list(out, seq_con_plot_pos,seq_con_plot_neg,seq_con_plot_both,
            seq_dis_plot_pos,seq_dis_plot_neg, seq_dis_plot_both))
}


#################################################################################
# Function to calculate prob variant detection if molecular follows ANCOV + / - 
#################################################################################

calc_variant_detection_prob_AN_and_DN = function( 
  PCR_sens = 0.92 ,
  Ag_sens_C = 0.643,
  Ag_sens_D = 0.0,
  Ag_spec = 0.99,
  Legend = F,
  y_axis_ticks = F,
  x_axis_ticks = F, 
  y_axis_label = F,
  percent_follow_up = 0.5,
  title) {
  
  
  
  prevelance = seq(0,1, by = 0.005)
  
  
  seq_perc = seq(0,0.03, by = 0.0001)
  
  prev_seq_perc = expand.grid(prevelance, seq_perc)
  

  
  # True positive Ag positive samples 
  
  TP_C_ag = (Ag_sens_C * prev_seq_perc$Var1) 
  TP_D_ag = (Ag_sens_D * prev_seq_perc$Var1)
  
  FN_C_ag = ((1- Ag_sens_C) * prev_seq_perc$Var1) 
  FN_D_ag = ((1- Ag_sens_D) * prev_seq_perc$Var1)

  
  
  Prob_genome_detection_C_pos = TP_C_ag * percent_follow_up * PCR_sens * prev_seq_perc$Var2
  Prob_genome_detection_D_pos = TP_D_ag * percent_follow_up * PCR_sens * prev_seq_perc$Var2
  

  Prob_genome_detection_C_neg = FN_C_ag * percent_follow_up * PCR_sens * prev_seq_perc$Var2 
  Prob_genome_detection_D_neg = FN_D_ag * percent_follow_up * PCR_sens * prev_seq_perc$Var2
  
  

  out = data.frame(prev_seq_perc, Prob_genome_detection_C_pos, Prob_genome_detection_D_pos 
                    ,  Prob_genome_detection_C_neg, Prob_genome_detection_D_neg) * 100
  
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
  
  if (y_axis_label == F){
    y_label = c(" ", "")
  } else (
    y_label = c("Cases\ndetected (%)")
  )
  
  seq_con_plot_pos = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_C_pos)) + 
    geom_tile() +
    scale_fill_distiller(palette = "RdYlBu" , limits = c(0,3), name= "Cases detected (%)") +
    labs(x = "", y = y_label[1]) +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position = leg
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    #  scale_x_continuous(limits = c(0.00, 0.05), breaks = seq(0.00, 0.05, 0.005))  +
    geom_contour(aes(z = Prob_genome_detection_C_pos), 
                 breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = Prob_genome_detection_C_pos), 
                 breaks = 1, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = Prob_genome_detection_C_pos), 
                 breaks = 2, col = 'black', size = .7, linetype = 3)
  
  
  
  seq_dis_plot_pos = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_D_pos)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdYlBu", limits = c(0,3), name= "Cases\ndetected (%)n") +
    labs(x = "", y = y_label[2]) +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position = "none"
    ) +
    ggtitle(title[1]) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    # scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
    geom_contour(aes(z = Prob_genome_detection_D_pos), 
                 breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = Prob_genome_detection_D_pos), 
                 breaks = 1, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = Prob_genome_detection_D_pos), 
                 breaks = 2, col = 'black', size = .7, linetype = 3)
  
  
  
  
  seq_con_plot_neg = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_C_neg)) + 
    geom_tile() +
    scale_fill_distiller(palette = "RdYlBu" , limits = c(0,3), name= "Cases\ndetected (%)") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    # scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
    geom_contour(aes(z = Prob_genome_detection_C_neg), 
                 breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = Prob_genome_detection_C_neg), 
                 breaks = 1, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = Prob_genome_detection_C_neg), 
                 breaks = 2, col = 'black', size = .7, linetype = 3)
  
  
  
  seq_dis_plot_neg = ggplot(out, aes(Var2 ,Var1 , fill= Prob_genome_detection_D_neg)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdYlBu", limits = c(0,3), name= "Cases\ndetected (%)") +
    labs(x = "", y = "") +
    theme_bw() + theme(
      text = element_text(size = 14),
      axis.title.y = element_text(angle = 90, vjust = 0.7),
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position =  "none"
    ) +
    ggtitle(title[2])+
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    #  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.005))  +
    geom_contour(aes(z = Prob_genome_detection_D_neg), 
                 breaks = 0.1, col = 'black', size = .7, linetype = 2)+ 
    geom_contour(aes(z = Prob_genome_detection_D_neg), 
                 breaks = 1, col = 'black', size = .7, linetype = 1)+
    geom_contour(aes(z = Prob_genome_detection_D_neg), 
                 breaks = 2, col = 'black', size = .7, linetype = 3)
  
  

  

  return(list(out, seq_con_plot_pos,seq_con_plot_neg,seq_dis_plot_pos,seq_dis_plot_neg))
}


