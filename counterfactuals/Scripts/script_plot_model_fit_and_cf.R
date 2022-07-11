
# script to plot counterfactual results (Figure 5) -----------------------------  


# Set up -----------------------------------------------------------------------


# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

library("wesanderson")
library(cowplot) 
library(tidyverse)
library(scales)
source("R/plot_model_fit.R")


#  Import  data  -----------------------------------------------
file_names_main = dir("MF_results/baseline") 
list_main = lapply(paste0("MF_results/baseline/",file_names_main),read.csv, )
names(list_main) = file_names_main

file_names_cf = dir("CF_Results") 
list_cf = lapply(paste0("CF_Results/",file_names_cf),read.csv, )
names(list_cf) = file_names_cf


# Plot main models -------------------------------------------------------------

 Veneto_plot = plot_model_fit(
  list_main$Veneto_post_df.csv,
  start_date = "01-05-2020",
  end_date = "31-05-2021",
  location = "Veneto"
)

Italy_plot = plot_model_fit(
  list_main$Italy_post_df.csv,
  start_date = "01-05-2020",
  end_date = "31-05-2021",
  location = "Italy"
)


#  Plot proportion of incidence reported in Veneto and Italy ---------------------------------------------------------

hex_codes1 = hue_pal()(4)   
hex_codes2 = c("#F8766D", "#7CAE00" , "#00BFC4" ,"#C77CFF", "#000000")


Prob_det_plot = list_main$Veneto_ratio.csv %>%
  bind_rows(list_main$Italy_ratio.csv) %>%
  filter(grepl('ratio', variant)) %>%
  pivot_longer(cols = -c(variant, Location)) %>%
  separate(variant, into = c("ratio", "variant")) %>%
  mutate(variant = ifelse(
    variant == "A",
    "A220V model",
    ifelse(variant == "M",
      "M234I-A376T model",
      ifelse(variant == "O",
        "Other model",
        ifelse(variant == "Al", "Alpha model", "NA"))) )) %>%
  mutate(variant = factor(
    variant, levels = c(
      "M234I-A376T model",
      "A220V model",
      "Other model",
      "Alpha model"
    ) ),  value = value * 100) %>%
  pivot_wider(
    id_cols = c(Location, variant),
    names_from = name,
    values_from = value  )  %>%
  ggplot(aes(x = Location, y = mean)) +
  geom_point(
    aes(color = variant),  shape = 19, size = 4,
    position = position_dodge(width = 0.5)  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, color = variant),
    position = position_dodge(width = 0.5),
    width =  0.4,  size = 1) +
  labs(x = " ",
       y = paste0("Mean probability of case detection \n during the modelling study period")) +
  theme_bw() + theme(
    text = element_text(size = 18),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c(0.22, 0.88), legend.title = element_blank(),
    axis.text.x = element_text( size =20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(
    name = '',values =
      c("A220V model" = hex_codes2[1],
        "Alpha model" = hex_codes2[2],
        "M234I-A376T model" = hex_codes2[3],
        "Other model" = hex_codes2[4]) )

list_main$Italy_ratio.csv
list_main$Veneto_ratio.csv
################  Plot cumulative incidence by test strategy ################  

Testing_cf_plot =  list_main$Veneto_ratio.csv %>%  
  bind_rows(list_cf$Veneto_ratio_cf4.csv,
            list_cf$Veneto_ratio_cf2.csv,
            list_cf$Veneto_ratio_cf3.csv,
            list_cf$Veneto_ratio_cf5.csv) %>%  
  filter(grepl('true', variant)) %>%
  mutate(model = factor(rep(c("molecular follows \nantigen+ \n(baseline)",
                              "molecular follows \nantigen-", 
                              "only antigen \nmass testing \n 68.9% test sensitivity",
                              "only antigen \nmass testing \n 87.5% test sensitivity", 
                              "Italy testing"), each = 5),
                       levels = c("molecular follows \nantigen+ \n(baseline)",
                                  "molecular follows \nantigen-",
                                  "Italy testing",
                                  "only antigen \nmass testing \n 68.9% test sensitivity",
                                  "only antigen \nmass testing \n 87.5% test sensitivity") )) %>% 
  separate(variant, into = c("ratio", "variant")) %>%
  mutate(variant = ifelse(variant == "A","A220V model",
    ifelse(variant == "M","M234I-A376T model",
      ifelse(variant == "O","Other model",
        ifelse(variant == "Al", "Alpha model", "Total model"))) )) %>%
  mutate(variant = factor(
    variant, levels = c(
      "M234I-A376T model",
      "A220V model",
      "Other model",
      "Alpha model",
      "Total model"
    ))) %>%  
  ggplot(aes(x = model, y = mean)) +
  geom_point(
    aes(color = variant),
    shape = 19, size = 4,
    position = position_dodge(width = 0.5) ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, color = variant),
    position = position_dodge(width = 0.5),
    width =  0.4, size = 1 ) + 
  labs(x = " ", y = paste0("Cumulative incidence over \n the modelling period (%)")) +
  theme_bw() + theme(
    text = element_text(size = 18),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c(0.07, 0.77),
    legend.title=element_blank(), 
    legend.text=element_text(size=14) ,
        axis.text.x = element_text(size = 20))+
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  scale_color_manual(
    name = '',values =
      c("A220V model" = hex_codes2[1],
        "Alpha model" = hex_codes2[2],
        "M234I-A376T model" = hex_codes2[3],
        "Other model" = hex_codes2[4],
        "Total model"  = hex_codes2[5]) )


################  Plot cumulative incidence by R0m ################  

Transmission_cf_plot = list_cf$Veneto_ratio_cf6.csv %>%  
  bind_rows(list_cf$Veneto_ratio_cf7.csv,list_cf$Veneto_ratio_cf8.csv ) %>%
  filter(grepl('true', variant)) %>%
  mutate(model = factor(rep(c("molecular follows \nantigen+ \n(baseline)",
                              "only antigen \nmass testing", "Italy testing"), each = 25), 
                        levels = c("molecular follows \nantigen+ \n(baseline)","Italy testing",
                                   "only antigen \nmass testing") )) %>% 
  separate(variant, into = c("metric", "variant")) %>%
  filter(variant  == "M" | variant == "tot") %>%  
  mutate(Variant =  ifelse(variant == "M","M234I-A376T \nmodel", "total \nmodel")) %>%
  mutate(Variant = factor(
    Variant, levels = c(
      "M234I-A376T \nmodel",
      "total \nmodel"
    ))) %>% 
  mutate("R0m scaling factor" = factor(R0_scale)) %>%  
  ggplot(aes(x = model, y = mean)) +
  geom_point(shape = 19, size = 4,aes(group = interaction(Variant,`R0m scaling factor`),
                                       color = `R0m scaling factor`, shape = Variant), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper,group = interaction(Variant,`R0m scaling factor`),
                    color = `R0m scaling factor`, linetype = Variant),
                position = position_dodge(width = 0.5),  width =  0.4, size = 1) +
  labs(x = " ", y = paste0("Cumulative incidence over \n the modelling period (%)")) +
  theme_bw() + theme(
    text = element_text(size = 18),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    legend.position = c(0.06,0.54),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14),
    axis.text.x = element_text(size=20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(values = wes_palette(name = "Darjeeling1")) 




################# Save ################# 

Baseline_fig = plot_grid(
  Veneto_plot,
  Italy_plot,
  Prob_det_plot,
  ncol = 3,
  labels = c("a", "b", "c"),
  align = c("h")
)

Fig5 = plot_grid(
  Baseline_fig,
  Testing_cf_plot,
  Transmission_cf_plot,
  ncol = 1,
  labels = c("", "d", "e"),
  rel_heights = c(1 /3, 1 / 4, 1 / 4),
  align = "v,h"
)


ggsave(
  plot = Fig5,
  filename = "Fig5.jpg",
  height = 40,
  width = 50,
  units = "cm",
  dpi = 800
)


