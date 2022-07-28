
# Script to plot R0 as a function of proportion of molecular tests (Figure S7) -----


# Set up -----------------------------------------------------------------------


# rm(list = ls())

setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")
library(tidyverse)
library(ggpubr)
library(cowplot)


# Import data ------------------------------------------------------------------



R0_Veneto = read.csv("MF_results/baseline/Veneto_R0.csv")[, -1]
R0_Veneto$Date = as.Date.character(R0_Veneto$Date, format = "%Y-%m-%d")

R0_Italy = read.csv("MF_results/baseline/Italy_R0.csv")[, -1]
R0_Italy$Date = as.Date.character(R0_Italy$Date, format = "%Y-%m-%d")

Veneto_pMO = read.csv("MF_results/baseline/Veneto_pMO.csv")
Veneto_pMO$Date = as.Date.character(Veneto_pMO$Date, format = "%Y-%m-%d")
Veneto_pMO$variable = "pMO"
Italy_pMO = read.csv("MF_results/baseline/Italy_pMO.csv")
Italy_pMO$Date = as.Date.character(Italy_pMO$Date, format = "%Y-%m-%d")


# Extract R0 at t0 for each variant --------------------------------------------

R0_Veneto %>%  
  filter(time == 1) %>%  
  mutate(lower = round(lower, 2), 
         mean = round(mean,2),
         upper = round(upper,2))


R0_Italy %>%  
  filter(time == 1) %>%  
  mutate(lower = round(lower, 2), 
         mean = round(mean,2),
         upper = round(upper,2))

# Plot R0 for Veneto -----------------------------------------------------------


R0_Veneto_plot =R0_Veneto %>% 
  mutate(variant = ifelse(
    variant == "R0_A",
    "A220V R0",
    ifelse(variant == "R0_M",
           "M234I-A376T R0",
           ifelse(variant == "R0_O",
                  "Other R0",
                  ifelse(variant == "R0_Al", "Alpha R0", "NA"))) )) %>%  
ggplot(aes(x = Date, y = mean    )) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variant),
              alpha = 0.5) +
  geom_line(aes(color = variant), size = 1) +
  geom_line(data = Veneto_pMO, aes(y = mean, linetype = pMO), color = "black") +
  labs(x = " ",
       y = paste0(" ")) +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    legend.position = c(0.16, 0.88),
    legend.key.height = unit(.5, 'cm'),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm")
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Veneto") + 
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_linetype_manual(values = c("pMO" = "solid")) 


# Plot R0 for Italy  -----------------------------------------------------------



R0_Italy_plot = ggplot(R0_Italy, aes(x = Date, y = mean    ))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variant),
              alpha = 0.5) +
  geom_line(aes(color = variant), size = 1) +
  geom_line(data = Italy_pMO, aes(y = mean), color = "black") +
  labs(x = " ",
       y = paste0(" ")) +
  theme_bw() + theme(
    text = element_text(size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.7),
    plot.title = element_text(hjust = 0.4),
    legend.position = c("none")
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months") +
  ggtitle("Rest of Italy") + 
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) 




################# save ################# 



R0_fig = plot_grid(R0_Veneto_plot, R0_Italy_plot, ncol = 2, labels = c("a", "b"))
R0_fig_anot = annotate_figure(R0_fig, bottom = "Date", element_text(size = 30))

ggsave(
  plot = R0_fig_anot,
  filename = "figS7.jpg",
  height = 15,
  width = 30,
  units = "cm",
  dpi = 300
)



