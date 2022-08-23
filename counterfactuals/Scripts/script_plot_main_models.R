################################################################################
# Script to run and plot variants of model testing (sup figure 4)              #
################################################################################


# Set up -----------------------------------------------------------------------


# rm(list = ls())
# setwd("C:/Users/bnc19/Desktop/COV_Italy_multistrain/counterfactuals")

# packages 
library("wesanderson")
library(cowplot) 
library(scales)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# functions 
source("R/sample_posterior_chains.R")
source("R/run_cf.R")
source("R/plot_model_fit.R")

# data 
posterior_chains_symp = read.csv("MF_results/symp_test/1/posterior_chains.csv")
posterior_chains_asymp = read.csv("MF_results/asymp_test/1/posterior_chains.csv")
posterior_chains_asymp2 = read.csv("MF_results/asymp_test_2/1/posterior_chains.csv")
posterior_chains_asymp_symp = read.csv("MF_results/asymp_and_symp_test/1/posterior_chains.csv")

# output 
file_path = "Figures"
dir.create(paste0(file_path))

# models 
symp_model = stan_model("models/est_test_symp_cf.stan")
asymp_model = stan_model("models/est_test_asymp_cf.stan")
asymp2_model = stan_model("models/est_test_asymp2_cf.stan")
asymp_symp_model = stan_model("models/est_test_asymp_and_symp_cf.stan")

list_models = list(symp_model,
                   asymp_model,
                   asymp2_model,
                   asymp_symp_model)
# Sample from posterior distribution -------------------------------------------

posterior_samples_symp = sample_posterior_chains(posterior_chains = posterior_chains_symp,
                                            number_of_samples = 100)
posterior_samples_asymp = sample_posterior_chains(posterior_chains = posterior_chains_asymp,
                                                 number_of_samples = 100)
posterior_samples_asymp2 = sample_posterior_chains(posterior_chains = posterior_chains_asymp2,
                                                 number_of_samples = 100)
posterior_samples_asymp_symp = sample_posterior_chains(posterior_chains = posterior_chains_asymp_symp,
                                                 number_of_samples = 100)

list_samples = list(posterior_samples_symp,
                    posterior_samples_asymp, 
                    posterior_samples_asymp2,
                    posterior_samples_asymp_symp)
# Run model across posterior samples -------------------------------------------


model_posts = lapply(1:length(list_samples), function(j){
  sapply(1:nrow(posterior_samples_symp), function(i) {
  replicate_rstan_fixed(
    model = list_models[[j]],
    posterior_sample_row = list_samples[[j]][i,],
    scale_time_step = 2,
    start_date_veneto =  "01-05-2020",
    time_seed_M_veneto = "01-08-2020"
  )
})})

Italy_posts = lapply(1:length(model_posts), function(j){
  sapply(1:nrow(posterior_samples_symp), function(i) {
  extract_fit_results(posts = model_posts[[j]][, i],
                      location = "Italy",
                      model = "SA")
}, simplify = "array")})


Veneto_posts =  lapply(1:length(model_posts), function(j){
  sapply(1:nrow(posterior_samples_symp), function(i) {
  extract_fit_results(posts = model_posts[[j]][, i],
                      location = "Veneto",
                      model = "SA")
}, simplify = "array")})


# summarise results and plot ---------------------------------------------------

Veneto_post_df = lapply(1:length(Veneto_posts), function(j){
      summarise_results(Veneto_posts[[j]],
                                   start_date = "01-05-2020",
                                   end_date = "31-05-2021",
                                   S0 = 4753625,
                                   no.col = 4)})

Italy_post_df = lapply(1:length(Italy_posts), function(j){
          summarise_results(Italy_posts[[j]],
                                  start_date = "01-05-2020",
                                  end_date = "31-05-2021",
                                  S0 = 53021564,
                                  no.col = 4)})



# plot model fits --------------------------------------------------------------


Veneto_data_plot = lapply(1:length(Veneto_post_df), function(j){
  plot_model_fit(
  Veneto_post_df[[j]],
  start_date = "01-05-2020",
  end_date = "31-05-2021",
  location = "Veneto"
)[[2]]})

Italy_data_plot = lapply(1:length(Italy_post_df), function(j){
  plot_model_fit(
  Italy_post_df[[j]],
  start_date = "01-05-2020",
  end_date = "31-05-2021",
  location = "Rest of Italy"
)[[2]]})


y = "Reported incidence per 100,000 population"

# top row 
vp1 = ggplot(Veneto_data_plot[[1]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
    ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = paste0(y), x = " ") +  theme_bw() +
  ggtitle("Veneto") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position =c(0.11,0.82), 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")

ip1 = ggplot(Italy_data_plot[[1]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = "", x = " ") +  theme_bw() +
  ggtitle("Rest of Italy") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")


# second row
vp2 = ggplot(Veneto_data_plot[[2]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = paste0(y), x = " ") +  theme_bw() +
  ggtitle(" ") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")

ip2 = ggplot(Italy_data_plot[[2]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = "", x = " ") +  theme_bw() +
  ggtitle(" ") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")

# third row
vp3 = ggplot(Veneto_data_plot[[3]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = paste0(y), x = " ") +  theme_bw() +
  ggtitle(" ") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")

ip3 = ggplot(Italy_data_plot[[3]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = "", x = " ") +  theme_bw() +
  ggtitle(" ") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")


# fourth row
vp4 = ggplot(Veneto_data_plot[[4]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = paste0(y)) +  theme_bw() +
  ggtitle(" ") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")

ip4 = ggplot(Italy_data_plot[[4]], aes(x = Date , y = mean_rep)) +
  geom_line(aes(color = Variant, linetype = Model)) +
  geom_ribbon(aes(ymin = lower_rep,
                  ymax = upper_rep, fill = Variant), alpha = 0.4) +
  geom_point(aes(y = PointEst , color = variant.x, shape = Data)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = variant.x )) +
  labs(y = "") +  theme_bw() +
  ggtitle(" ") +  ylim(c(0, 150)) + 
  scale_shape_manual(values = c('Data' = 16)) +
  scale_linetype_manual(values = c("Model" = "solid")) +
  theme(text = element_text(size = 18), legend.position = "none", 
        legend.margin = margin(0, 0, 0, 0),  legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), legend.title = element_blank()) +
  scale_x_date(date_labels = "%b/%Y", breaks = "3 months")



FigS6 = plot_grid(vp1,ip1,vp2,ip2,vp3,ip3,vp4,ip4, 
                ncol = 2, 
                labels = c("a","b","c","d","e","f","g","h"),
                align = "v,h")


ggsave(
  plot = FigS6,
  filename = "Figures/FigureS6.pdf",
  height = 65,
  width = 50,
  units = "cm",
  dpi = 800
)



