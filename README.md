# COV_Italy_multistrain

## Repository structure

This repository is divided into 2 sections: 

* *model_fitting* - code to fit a multivariant model to reconstruct the transmission dynamics of the dominant variants carrying the A220V and M234I-A376T mutations, all other co-circulating variants, and the alpha variant (B.1.1.7) for Veneto (between July 2020-May 2021) and the rest of Italy (between May 2020-2021), separately. 
* *counterfactuals*- code to run models of the baseline and counterfactual testing scenarios using posterior samples obtained from *model_fitting* as input. Also includes code to run the sensitivitiy analyses, code to calculate test performance metrics and estimate varaint detection by genomic surveillance, and code to produce figures from the manuscript. 

### model_fitting

This folder contains

Data

* *Dataset_Italy_A_v5.csv*, *Dataset_Italy_Alpha_v1.csv*, *Dataset_Italy_M_v5.csv*, *Dataset_Italy_O_v1.csv* : 
For Italy, data on the number of sequences tested, the number of sequences testing positive by variant, obtained from the GISAID databank (1), data only the total monthly reported incidence and mean daily reported incidence in month m, obtained obtained from the Civil Protection (2).  

* *Dataset_Veneto_A_v5.csv*, *Dataset_Veneto_Alpha_v1.csv*, *Dataset_Veneto_M_v5.csv*, *Dataset_Veneto_O_v1.csv*:
For Veneto, data on the number of sequences tested, the number of sequences testing positive by variant, obtained from the GISAID databank (1), data only the total monthly reported incidence and mean daily reported incidence in month m, obtained obtained from the Civil Protection (2).  

* *Italy_daily_test_data.csv*, *Italy_monthly_test_data.csv* : 
For Italy,  data on the daily reported number of antigen and molecular tests and mean daily reported number of antigen and molecular tests in month m, obtained from the Civil Protection (2). 

* *Veneto_daily_test_data.csv*, *Veneto_monthly_test_data.csv* : 
For Veneto,  data on the daily reported number of antigen and molecular tests and mean daily reported number of antigen and molecular tests in month m, obtained from the Civil Protection (2). 

* *dailyReportedIncidence_italy.csv* :
For Italy, data on the daily reported incidence, obtained obtained from the Civil Protection (2).  

* *dailyReportedIncidence_veneto.csv* :
For Veneto, data on the daily reported incidence, obtained obtained from the Civil Protection (2).  

* *Vac_Italy_For_Month.csv*
For Italy, data on the mean daily per-capita number of second doses administered obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).

* *Vac_Veneto_For_Month.csv*
For Veneto, data on the mean daily per-capita number of second doses administered obtained from the Extraordinary Commissioner for the Covid-19 emergency (3).



### References 
(1) Elbe, S. & Buckland‐Merrett, G. Data, disease and diplomacy: GISAID’s innovative contribution to global health. Global challenges 1, 33-46 (2017).
(2) Civile, P.d.C.d.M.-D.d.P. Dati COVID-19 Italia. 2021.
(3) Commissario straordinario per l'emergenza Covid-19 - Presidenza del Consiglio dei Ministri. Open Data su consegna e somministrazione dei vaccini anti COVID-19 in Italia.  2021 16th November 2021 [cited  3 March 2022]Available from: https://github.com/italia/covid19-opendata-vaccini
