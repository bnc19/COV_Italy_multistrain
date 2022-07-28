
library(usethis)
use_description(check_name = F,
                roxygen = F,
                fields = list(
                  Title = " Impact of antigen test target failure and testing strategies on the transmission of SARS-CoV-2 variants",
                  Description = "Code to reproduce all results from the manuscript Impact of antigen test target failure and testing strategies on the transmission of SARS-CoV-2 variants. 
                  To learn more about this project, start with README.md",
                  `Authors@R` =  person("Bethan", "Cracknell Daniels", email = "bnc19@ic.ac.uk", 
                                        role = c("aut", "cre")),
                  Version = "1.0"
                )
)
usethis::use_package("bayesplot", min_version ="1.8.1")
usethis::use_package("cowplot" , min_version = "1.1.1")
usethis::use_package("deSolve" , min_version = "1.2.8")
usethis::use_package("dplyr", min_version = "1.0.6")
usethis::use_package("ggplot2" , min_version = "3.3.6")
usethis::use_package("ggpubr" , min_version = "0.4.0")
usethis::use_package("Hmisc" , min_version = "4.5-0")
usethis::use_package("rstan" , min_version = "2.21.5")
usethis::use_package("scales" , min_version = "1.1.1")
usethis::use_package("tidyr" , min_version = "1.1.3")
usethis::use_package("wesanderson" , min_version = "0.3.6")
usethis::use_package("R" , min_version = "4.1.1", "Depends")