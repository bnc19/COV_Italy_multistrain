
################  Plot data on number of antigen and molecular tests in Veneto and Italy (Figure 1) ################   

################    Set up ################   


# rm(list = ls())

################  Import data ###############

library(Hmisc)
library(tidyverse)

  
Italy_test_data = read.csv("model_fitting/data/Italy_daily_test_data.csv") 
Italy_test_data$total = Italy_test_data$pcr_daily_average + Italy_test_data$antigen_daily_average

Veneto_test_data = read.csv("model_fitting/data/Veneto_daily_test_data.csv") 
Veneto_test_data$total = Veneto_test_data$pcr_daily_average + Veneto_test_data$antigen_daily_average


################  Format and plot ###############


Veneto_test_conf = data.frame(Date = Veneto_test_data$Date, 
                              Hmisc::binconf(Veneto_test_data$antigen_daily_average, Veneto_test_data$total) * 100 ) 

Italy_test_conf =  data.frame(Date = Italy_test_data$Date, 
                              Hmisc::binconf(Italy_test_data$antigen_daily_average, Italy_test_data$total)* 100)


Veneto_test_conf %>% 
  filter(Date >= "2020-09-01") %>% 
  summarise(mean = mean (PointEst))

Italy_test_conf %>% 
  filter(Date >= "2020-09-01") %>% 
  summarise(mean = mean (PointEst))

Test_conf = Veneto_test_conf %>% 
  bind_rows(Italy_test_conf) %>% 
  mutate(Location = factor(rep(c("Veneto", "Rest of Italy"), each = 426), levels = c("Veneto", "Rest of Italy")),
         Date = as.Date.character(Date)) %>%  
  ggplot(aes(x = Date, y = PointEst )) +
  geom_line(aes(color = Location)) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Location), alpha =0.2) + 
  labs(y = "Antigenic diagnostic tests conducted (%)") +
  theme_bw() + theme(
    text = element_text(size = 16),
    legend.position = c(0.15, 0.9),
    legend.title = element_blank(),
    legend.key.height = unit(.5, 'cm')) +
  scale_x_date(date_labels = "%b/%Y", breaks = "2 months")+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))


################  Save ###############

ggsave(
  plot = Test_conf,
  filename = "fig1.jpg",
  height = 15,
  width = 20,
  units = "cm",
  dpi = 300
)
