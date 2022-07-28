################################################################################
# Script to plot all data used for model fitting (i.e., data on number of      #
# antigen and molecular tests, the reported incidence, the reported GISAID     #
# variant prevalence and the reconstructed variant reported incidence in       #
# Veneto and the rest of Italy (Figure 3))                                     #                                                        #
################################################################################


# Set up -----------------------------------------------------------------------

# setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")
# rm(list = ls())

# packages 
library(cowplot)
library(tidyverse)
dir.create("Figures")

# define variables 
n_pop_it = 59257566 - 4847026
n_recov_it = 1482377 - 93401
n_pop_veneto = 4847026
n_recov_veneto = 93401
S0_it = n_pop_it - n_recov_it
S0_ven = n_pop_veneto - n_recov_veneto

# Import Italy data ------------------------------------------------------------
A_data_it = read.csv("data/Dataset_italy_A_v5.csv")
M_data_it = read.csv("data/Dataset_italy_M_v5.csv")
O_data_it = read.csv("data/Dataset_italy_O_v1.csv")
Al_data_it = read.csv("data/Dataset_italy_Alpha_v1.csv")

Italy_incidence = read.csv("data/dailyReportedIncidence_italy.csv")

Italy_test_data = read.csv("data/Italy_daily_test_data.csv")


average_daily_vaccination_it  =  read.csv("data/data_vac_italy_day.csv")$prop_vac

# Import Veneto data -----------------------------------------------------------
A_data_veneto = read.csv("data/Dataset_Veneto_A_v5.csv")
M_data_veneto  = read.csv("data/Dataset_Veneto_M_v5.csv")
O_data_veneto  = read.csv("data/Dataset_Veneto_O_v1.csv")
Al_data_veneto  = read.csv("data/Dataset_Veneto_Alpha_v1.csv")
Al_data_veneto$Perc = as.numeric(Al_data_veneto$Perc)
Veneto_incidence = read.csv("data/dailyReportedIncidence_veneto.csv")

Veneto_test_data = read.csv("data/Veneto_daily_test_data.csv")


# combine all variant data ----------------------------------------------------
Variant_data = bind_rows(
  M_data_it,
  A_data_it,
  O_data_it,
  Al_data_it,
  M_data_veneto,
  A_data_veneto,
  O_data_veneto,
  Al_data_veneto
)

Variant_data = Variant_data %>%
  mutate(
    Reported_incidence = ifelse(
      Country == "Italy",
      new_reported_cases_daily_new / S0_it  * 100000,
      new_reported_cases_daily     / S0_ven * 100000
    ),
    percentage = ifelse(Country == "Italy",
                        Freq_new  / TotSeq_new,
                        Freq / TotSeq),
    ,
    Var_rep_inc = percentage * Reported_incidence
  ) %>%
  mutate(Country = ifelse(Country == "Italy", "Rest of Italy", "Veneto")) %>%
  mutate(Date = as.Date.character(Date, format = "%d/%m/%Y")) %>%
  filter(Date > "2020-04-01")



# Plot prevalence in GISAID ----------------------------------------------------
M_prev = Variant_data %>%
  select(Date, Mut, Country, percentage) %>%
  filter(Mut == "M234I-A376T") %>%
  ggplot(aes(x = Date, y = percentage * 100)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "Prevalence in GISAID (%)") +
  theme_bw() + theme(
    text = element_text(size = 16),
    legend.position = c("none"),
    legend.title = element_blank(),
    legend.key.height = unit(.5, 'cm')
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("")

A_prev = Variant_data %>%
  select(Date, Mut, Country, percentage) %>%
  filter(Mut == "A220V") %>%
  ggplot(aes(x = Date, y = percentage * 100)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "  ") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("")

O_prev = Variant_data %>%
  select(Date, Mut, Country, percentage) %>%
  filter(Mut == "Other") %>%
  ggplot(aes(x = Date, y = percentage * 100)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("")


Al_prev = Variant_data %>%
  select(Date, Mut, Country, percentage) %>%
  filter(Mut == "Alpha") %>%
  ggplot(aes(x = Date, y = percentage * 100)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("")


# Plot variant specific reported incidence -------------------------------------


M_inc = Variant_data %>%
  select(Date, Mut, Country, Var_rep_inc) %>%
  filter(Mut == "M234I-A376T") %>%
  ggplot(aes(x = Date, y = Var_rep_inc)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "Reported incidence \n per 100,000 population", x = "") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("M234I-A376T")

A_inc = Variant_data %>%
  select(Date, Mut, Country, Var_rep_inc) %>%
  filter(Mut == "A220V") %>%
  ggplot(aes(x = Date, y = Var_rep_inc)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = " ", x = "") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("A220V")

O_inc = Variant_data %>%
  select(Date, Mut, Country, Var_rep_inc) %>%
  filter(Mut == "Other") %>%
  ggplot(aes(x = Date, y = Var_rep_inc)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "", x = "") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("Other")


Al_inc = Variant_data %>%
  select(Date, Mut, Country, Var_rep_inc) %>%
  filter(Mut == "Alpha") %>%
  ggplot(aes(x = Date, y = Var_rep_inc)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(y = "" , x = "") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("Alpha")

# plot grid of variant specific data -------------------------------------------

variant_plot = plot_grid(
  M_inc,
  A_inc,
  O_inc,
  Al_inc,
  M_prev,
  A_prev,
  O_prev,
  Al_prev,
  labels = c("c", "d", "e", "f", "g", "h", "i", "j"),
  ncol = 4,
  align = "h"
)

# plot total incidence --------------------------------------------------------------

total_inc = Variant_data %>%
  select(Date, Mut, Country, Reported_incidence) %>%
  filter(Mut == "Alpha") %>%
  ggplot(aes(x = Date, y = Reported_incidence)) +
  geom_line(aes(colour = Country), size = 1.5) +
  labs(x = "", y = "Reported incidence \n per 100,000 population") +
  theme_bw() + theme(text = element_text(size = 16),
                     legend.position = c("none")) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  ggtitle("")

# Plot proportion diagnostic tests that are antigen ----------------------------
Italy_test_data$total = Italy_test_data$pcr_daily_average + Italy_test_data$antigen_daily_average
Italy_test_data$Country = "Rest of Italy"

Veneto_test_data$total = Veneto_test_data$pcr_daily_average + Veneto_test_data$antigen_daily_average
Veneto_test_data$Country = "Veneto"


Test_conf = Italy_test_data %>%
  bind_rows(Veneto_test_data) %>%
  mutate(Date = as.Date.character(Date),
         perc_ant = antigen_daily_average / total * 100) %>%
  filter(Date >= "2020-05-01" & Date <  "2021-05-02") %>%
  ggplot(aes(x = Date, y = perc_ant)) +
  geom_line(aes(color = Country), size = 1.5) +
  labs(x = "", y = "Proportion of COVID-19 patients \n receiving antigen diagnostic tests") +
  theme_bw() + theme(
    text = element_text(size = 16),
    legend.position = c(0.08, 0.85),
    legend.title = element_blank()
  ) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))


# Save -------------------------------------------------------------------------

total_plots = plot_grid(Test_conf,
                        total_inc,
                        labels = c("a", "b"),
                        align = "h")

grid_out = plot_grid(total_plots,
                     variant_plot,
                     ncol = 1,
                     rel_heights = c(1 / 4, 1 / 2))


ggsave(
  plot = grid_out,
  filename = "Figures/Figure3.png",
  height = 35,
  width = 55,
  units = "cm",
  dpi = 300
)
