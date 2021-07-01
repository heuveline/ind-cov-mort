#################################
setwd("C:/Users/Heuveline/Box Sync/Documents/Covid/R files")
getwd()
library(dplyr)
library(foreign)
library(reshape2)
library(vcd)
######
###This routine uses extant life table estimates for 2020
###and calculate average differences between actual and potential age at death
##for cohorts male, females and both-sex deaths in the last 3 quarters of 2020
##and the first quarter of 2021
####
###naming conventions
###
###death refers to all deaths included covid deaths
###covid refers to covid-19 deaths
###prev refers to deaths previously projected for 2020 (without covid deaths)
###
##Q2, Q3 and Q4 refer to the last 3 quarters of 2020; Q1 refers to first quarter of 2021
####
#####Step 1: get data
###step 1.1:get estimates of covid deaths for Q2, Q3, Q4 & Q1
##cumulative numbers at quarter end were scraped on 7/3, 10/2, 1/1 & 3/31 from
##https://github.com/CSSEGISandData/COVID-19
##cumulative numbers at the beginning of Q2 were enterred from WHO Situation Report 72 for 4/1/2020
##https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200401-sitrep-72-covid-19.pdf
##Note: sub-national estimates were only available for China
##those for other countries were prorated from national and subnational estimates at a later date 
##
all_Q2_covid <- read.csv("ccdr_200703.csv") [,c(3:4,11)]
all_Q3_covid <- read.csv("ccdr_201002.csv") [,c(3:4,11)]
all_Q4_covid <- read.csv("ccdr_210101.csv") [,c(3:4,11)]
all_Q1_covid <- read.csv("ccdr_210331.csv") [,c(3:4,11,18)]
##
all_2Q_covid <- na.omit(merge(all_Q2_covid,all_Q3_covid[,-c(2)],"LocID"))
all_3Q_covid <- na.omit(merge(all_2Q_covid,all_Q4_covid[,-c(2)],"LocID"))
all_4Q_covid <- na.omit(merge(all_3Q_covid,all_Q1_covid[,-c(2)],"LocID"))
colnames(all_4Q_covid)[1:7] <- c("LocID", "Location","Covid_JHU_Q2","Covid_JHU_Q3",
                                 "Covid_JHU_Q4", "Covid_JHU_Q1", "Covid_WHO_Q1")
##
###
##step 1.2: get covid mortality by age and sex from CDC
##retrieved on 1/2/21 (corresponding to CDC 12/30/20 update) from
##https://data.cdc.gov/NCHS/Provisional-COVID-19-Death-Counts-by-Sex-Age-and-S/9bhg-hcku
##
usa_covid_11 <- read.csv("ccmr_4Q.csv")
##
##step 1.3: get population data by age and sex
##get population data by age and sex
##
##data for country population derived from UN Population Division data at
##https://population.un.org/wpp/Download/Standard/CSV/
##data for U.S. state population derived from Census Bureau data at
##https://data.census.gov
##data for China provinces dervied from China Statistical Office data at
##http://data.stats.gov.cn
##data for Brazil state population derived from Brazil Statistical Office
##https://www.ibge.gov.br/en/cities-and-states.html?view=municipio
## data for Italian regions and provinces from Istat 
##http://demo.istat.it/tvm2016/index.php
###data for Spanish autonomous communities and provinces
##from Instituto Nacional Estadistica (INE)
##https://www.ine.es/dyngs/INEbase/en/categoria.htm?c=Estadistica_P&cid=1254734710984
##data for Mexican States from Instituto Nacional de Estadistica y Geografía (INEGI)
##https://en.www.inegi.org.mx/temas/estructura/default.html#Tabulados
###data for Peruvian Departmientos from El Instituto Nacional de Estadísticañ 
###e Informática (INEI)
##https://www.inei.gob.pe
###
all_20_pop_11 <- read.csv("ccmr_fixed.csv")
#
##
###step 1.4:get estimated life table values for 2020 in all countries
##these life tables were derived from UN Population Division life tables
##for 2015-20 & 2020-25 retrieved from
##https://population.un.org/wpp/Download/Standard/CSV/
##
all_20_lt_11 <- read.csv("diff_e0_fixed.csv")
##                 
##step 2: calculate US age-specific covid death rates by age and sex
##
#get US population
#
usa_pop_11 <- all_20_pop_11[ which(all_20_pop_11$LocID=='USA'), ]
#
#get US covid data
#
usa_4Q_covid <- all_4Q_covid[ which(all_4Q_covid$LocID=='USA'), ]
#
#calculate mid-point of Covid infection to adjust exposure
#tbar is fraction of 2020 missed by Covid death
#
usa_4Q_covid$tbar = (usa_4Q_covid$Covid_WHO_Q1+usa_4Q_covid$Covid_JHU_Q2+
  usa_4Q_covid$Covid_JHU_Q3+(usa_4Q_covid$Covid_JHU_Q4/2))/(4*usa_4Q_covid$Covid_JHU_Q4)
#
#reshape CDC data by sex
#
usa_covid_11_m <- reshape(usa_covid_11[,c(1,3,7)], timevar = "age_group",
                          idvar = c("state"),direction = "wide")
usa_covid_11_f <- reshape(usa_covid_11[,c(1,4,7)], timevar = "age_group",
                          idvar = c("state"),direction = "wide")
usa_covid_11 <- cbind(usa_covid_11_m,usa_covid_11_f[,-c(1)])
#
#calculate adjustment to CDC data to get annual number of Covid deaths in USA
#
usa_4Q_covid$mult <- usa_4Q_covid$Covid_JHU_Q4/rowSums(usa_covid_11[,2:23])
#
#adjust population size to get better measure of exposure
#
usa_exp_11 <- as.vector(usa_pop_11[5:26])-(as.vector(usa_covid_11[,2:23])*
                    usa_4Q_covid$mult*usa_4Q_covid$tbar/1000)
#
##divide covid deaths by exposure for each age group to get deaths per person-year
##
usa_asdr_11 <- as.vector(usa_covid_11[,2:23])*usa_4Q_covid$mult/
  as.vector(usa_exp_11[,1:22])
##
##
##step 3:calculate Chiang ratios of death estimate with to without covid death
##for Q2, Q3, Q4 & Q1
##
##calculate estimated number of covid deaths in Q2, Q3, Q4 & Q1
#
all_4Q_covid$Covid_Q2 <- all_4Q_covid$Covid_JHU_Q2 - all_4Q_covid$Covid_WHO_Q1
all_4Q_covid$Covid_Q3 <- all_4Q_covid$Covid_JHU_Q3 - all_4Q_covid$Covid_JHU_Q2
all_4Q_covid$Covid_Q4 <- all_4Q_covid$Covid_JHU_Q4 - all_4Q_covid$Covid_JHU_Q3
all_4Q_covid$Covid_Q1 <- all_4Q_covid$Covid_JHU_Q1 - all_4Q_covid$Covid_JHU_Q4
#
##calculate number of deaths in absence of Covid in first quarter of 2020 (in thousands)
##for each sex and age group from nmx and population size
##quarter referred to as Q0 below
#
#in the first quarter of 2020, exposure without covid is 1/4 of population size
#
all_20_nmx_11 <- reshape(all_20_lt_11[,c(2:3,5:16)], 
                             timevar = "Sex",
                             idvar = c("LocID", "Location"),
                             direction = "wide")
#
#need to drop populations file <90,000 from population file
#
#make sure both files have same location names
#
all_Q0_dNoC_11 <- merge(all_20_nmx_11,all_20_pop_11[,-c(1,3:4)],"LocID")
all_Q0_dNoC_11_M <- as.matrix(all_Q0_dNoC_11[,c(3:24)])*
  (as.matrix(all_Q0_dNoC_11[,25:46])/4)
all_Q0_dNoC_11 <- cbind(all_Q0_dNoC_11[,1:2],all_Q0_dNoC_11_M)
colnames(all_Q0_dNoC_11)[3:24] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                        "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                        "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                        "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
##distribute covid deaths by sex and age group
##first using the asdr for the US and the population by age and sex
##to get the expected distribution in all areas
##
usa_asdr_11_M <- diag(usa_asdr_11,22,22)
#
#getting exposure by age in each location
#
all_Q0_covid_11 <- merge(all_4Q_covid[,c(1:2,7)],all_20_pop_11[,-c(1,3:4)],"LocID")
#
####################
#check the number of locations above for input in diagonal matrix below
####################
#multiplying by covid asdr
#
all_Q0_covid_11_M <- as.matrix(all_Q0_covid_11[,4:25]) %*% usa_asdr_11_M 
#
#prorate to get the total number projected for each quarter in each country
#
all_Q0_covid_11 <- cbind(all_Q0_covid_11[,c(1:3)],all_Q0_covid_11_M,na.rm = TRUE)
all_Q0_ratio <-  as.vector(all_Q0_covid_11$Covid_WHO_Q1) / 
  as.vector(rowSums(all_Q0_covid_11_M))
#
####################################################
##may need update if number of units changes
##current number of countries/states this quarter is 397
#includes Mexico, Peru & US States, plus Italian and Spanish Provinces/Areas
##
all_Q0_ratio_M <- diag(all_Q0_ratio,397,397)
##
all_Q0_covid_11_M <- as.matrix(all_Q0_ratio_M) %*% 
  as.matrix(all_Q0_covid_11[4:25])
all_Q0_covid_11 <- cbind(all_Q0_covid_11[,c(1:3)],all_Q0_covid_11_M)
all_Q0_covid_11$Time <- "2020_Q1"
colnames(all_Q0_covid_11)[4:25] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                      "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                      "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                      "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
##to calculate number of non-covid deaths, need to adjust exposure
##assuming covid deaths are mid-way through each quarter
#reduces exposure in Q1 by 1/2 of Q1 covid deaths/1,000
##
all_Q0_exp_11 <- merge(all_Q0_covid_11,all_20_pop_11[,-c(1,3:4)],"LocID")
all_Q0_exp_11_M <- (as.matrix(all_Q0_exp_11[,27:48]) -
  (as.matrix(all_Q0_exp_11[4:25])/2000) )/4
all_Q0_exp_11 <- cbind(all_Q0_covid_11[,c(1:3)],all_Q0_exp_11_M)
#
#multiply by mortality rates to get new number of deaths from non-covid causes
#
all_Q0_other_11 <- merge(all_20_nmx_11, all_Q0_exp_11[,-c(2:3)],"LocID")
all_Q0_other_11_M <- as.matrix(all_Q0_other_11[,c(3:24)]) *
  as.matrix(all_Q0_other_11[,c(25:46)])
all_Q0_other_11 <- cbind(all_Q0_other_11[,c(1:2)],all_Q0_other_11_M)
#
#calculate the new number of total deaths in first quarter of 2020
#
all_Q0_death_11 <- merge(all_Q0_other_11,all_Q0_covid_11[,-c(2:3)],"LocID")
all_Q0_death_11_M <- (as.matrix(all_Q0_death_11[,25:46])/1000) +
  as.matrix(all_Q0_death_11[,3:24])
all_Q0_death_11 <- cbind(all_Q0_death_11[,1:2],all_Q0_death_11_M)
#
#calculate the number of excess deaths in first quarter
#
all_Q0_dXS_11 <- merge(all_Q0_death_11,all_Q0_dNoC_11[,-c(2)],"LocID")
all_Q0_dXS_11_M <- as.matrix(all_Q0_dXS_11[,3:24]) - 
  as.matrix(all_Q0_dXS_11[,25:46])
all_Q0_dXS_11 <- cbind(all_Q0_dXS_11[,1:2],all_Q0_dXS_11_M)
#
#calculate Chiang's age-specific ratios in first quarter
#by dividing new total by previously projected
#
all_Q0_chratio_11 <- merge(all_Q0_death_11,all_Q0_dNoC_11[,-c(2)],"LocID")
all_Q0_chratio_11_M <- as.matrix(all_Q0_chratio_11[,3:24]) / 
  as.matrix(all_Q0_chratio_11[,25:46])
all_Q0_chratio_11 <- cbind(all_Q0_chratio_11[,1:2],all_Q0_chratio_11_M)
colnames(all_Q0_chratio_11)[3:24] <- c("0R1.M","4R1.M","10R5.M","10R15.M","10R25.M",
                    "10R35.M","10R45.M","10R55.M","10R65.M","10R75.M","R85+.M",
                   "0R1.F","4R1.F","10R5.F","10R15.F","10R25.F","10R35.F",
                   "10R45.F","10R55.F","10R65.F","1R75.F","R85+.F")
#
##second quarter
#the exposure in the absence of Covid has changed due to excess mortality
#in first quarter
#
all_Q2_expNoC_11 <- merge(all_Q0_dXS_11,all_20_pop_11[,-c(1,3:4)],"LocID")
all_Q2_expNoC_11_M <- (as.matrix(all_Q2_expNoC_11[,25:46]) -
                      as.matrix(all_Q2_expNoC_11[3:24]) )/4
all_Q2_expNoC_11 <- cbind(all_Q2_expNoC_11[,c(1:2)],all_Q2_expNoC_11_M)
#
all_Q2_dNoC_11 <- merge(all_20_nmx_11,all_Q2_expNoC_11[,-c(2)],"LocID")
all_Q2_dNoC_11_M <- as.matrix(all_Q2_dNoC_11[,c(3:24)])*
  as.matrix(all_Q2_dNoC_11[,25:46])
all_Q2_dNoC_11 <- cbind(all_Q2_dNoC_11[,1:2],all_Q2_dNoC_11_M)
colnames(all_Q2_dNoC_11)[3:24] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                    "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                    "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                    "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
##distribute covid deaths by sex and age group
#
#getting exposure by age and locations 
#
all_Q2_covid_11 <- merge(all_4Q_covid[,c(1,2,8)],all_Q2_expNoC_11[,-c(2)],"LocID")
#
#multiplying by covid asdr
#
all_Q2_covid_11_M <- as.matrix(all_Q2_covid_11[,4:25]) %*% usa_asdr_11_M 
#
#prorate to get the total number projected for each quarter in each location
#
all_Q2_covid_11 <- cbind(all_Q2_covid_11[,c(1:3)],all_Q2_covid_11_M,na.rm = TRUE)
all_Q2_ratio <-  as.vector(all_Q2_covid_11$Covid_Q2) / 
  as.vector(rowSums(all_Q2_covid_11_M))
##
####################################################
##number of locations has changed because now use expNoC from dXS in previous
##quarter instead of population size which had more locations
##current number of locations is 339
##
all_Q2_ratio_M <- diag(all_Q2_ratio,339,339)
##
all_Q2_covid_11_M <- as.matrix(all_Q2_ratio_M) %*% 
  as.matrix(all_Q2_covid_11[4:25])
all_Q2_covid_11 <- cbind(all_Q2_covid_11[,c(1:3)],all_Q2_covid_11_M)
all_Q2_covid_11$Time <- "2020_Q2"
colnames(all_Q2_covid_11)[4:25] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                     "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                     "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                     "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
##to calculate number of non-covid deaths, need to adjust exposure
##assuming covid deaths are mid-way through each quarter
#reduces exposure in Q2 by 1/2 of Q2 covid deaths/1,000
##
all_Q2_exp_11 <- merge(all_Q2_covid_11,all_Q2_expNoC_11[,-c(2)],"LocID")
all_Q2_exp_11_M <- as.matrix(all_Q2_exp_11[,27:48]) -
                      (as.matrix(all_Q2_exp_11[4:25])/8000)
all_Q2_exp_11 <- cbind(all_Q2_covid_11[,c(1:3)],all_Q2_exp_11_M)
#
all_Q2_other_11 <- merge(all_20_nmx_11, all_Q2_exp_11[,-c(2:3)],"LocID")
all_Q2_other_11_M <- as.matrix(all_Q2_other_11[,c(3:24)]) *
  as.matrix(all_Q2_other_11[,c(25:46)])
all_Q2_other_11 <- cbind(all_Q2_other_11[,c(1:2)],all_Q2_other_11_M)
#
all_Q2_death_11 <- merge(all_Q2_other_11,all_Q2_covid_11[,-c(2:3)],"LocID")
all_Q2_death_11_M <- (as.matrix(all_Q2_death_11[,25:46])/1000) +
  as.matrix(all_Q2_death_11[,3:24])
all_Q2_death_11 <- cbind(all_Q2_death_11[,1:2],all_Q2_death_11_M)
#
all_Q2_dXS_11 <- merge(all_Q2_death_11,all_Q2_dNoC_11[,-c(2)],"LocID")
all_Q2_dXS_11_M <- as.matrix(all_Q2_dXS_11[,3:24]) - 
  as.matrix(all_Q2_dXS_11[,25:46])
all_Q2_dXS_11 <- cbind(all_Q2_dXS_11[,1:2],all_Q2_dXS_11_M)
#
all_Q2_chratio_11 <- merge(all_Q2_death_11,all_Q2_dNoC_11[,-c(2)],"LocID")
all_Q2_chratio_11_M <- as.matrix(all_Q2_chratio_11[,3:24]) / 
  as.matrix(all_Q2_chratio_11[,25:46])
all_Q2_chratio_11 <- cbind(all_Q2_chratio_11[,1:2],all_Q2_chratio_11_M)
colnames(all_Q2_chratio_11)[3:24] <- c("0R1.M","4R1.M","10R5.M","10R15.M","10R25.M",
                                       "10R35.M","10R45.M","10R55.M","10R65.M","10R75.M","R85+.M",
                                       "0R1.F","4R1.F","10R5.F","10R15.F","10R25.F","10R35.F",
                                       "10R45.F","10R55.F","10R65.F","1R75.F","R85+.F")
#
##repeat for third and fourth quarter of 2020 and first quarter of 2021
#
all_Q3_expNoC_11 <- merge(all_Q2_dXS_11,all_Q2_expNoC_11[,-c(2)],"LocID")
all_Q3_expNoC_11_M <- as.matrix(all_Q3_expNoC_11[,25:46]) -
                         (as.matrix(all_Q3_expNoC_11[3:24]) /4)
all_Q3_expNoC_11 <- cbind(all_Q3_expNoC_11[,c(1:2)],all_Q3_expNoC_11_M)
#
all_Q3_dNoC_11 <- merge(all_20_nmx_11,all_Q3_expNoC_11[,-c(2)],"LocID")
all_Q3_dNoC_11_M <- as.matrix(all_Q3_dNoC_11[,c(3:24)])*
  as.matrix(all_Q3_dNoC_11[,25:46])
all_Q3_dNoC_11 <- cbind(all_Q3_dNoC_11[,1:2],all_Q3_dNoC_11_M)
colnames(all_Q3_dNoC_11)[3:24] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                    "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                    "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                    "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
#
all_Q3_covid_11 <- merge(all_4Q_covid[,c(1,2,9)],all_Q3_expNoC_11[,-c(2)],"LocID")
#
#multiplying by covid asdr
#
all_Q3_covid_11_M <- as.matrix(all_Q3_covid_11[,4:25]) %*% usa_asdr_11_M 
#
all_Q3_covid_11 <- cbind(all_Q3_covid_11[,c(1:3)],all_Q3_covid_11_M,na.rm = TRUE)
all_Q3_ratio <-  as.vector(all_Q3_covid_11$Covid_Q3) / 
  as.vector(rowSums(all_Q3_covid_11_M))
##
all_Q3_ratio_M <- diag(all_Q3_ratio,339,339)
##
all_Q3_covid_11_M <- as.matrix(all_Q3_ratio_M) %*% 
  as.matrix(all_Q3_covid_11[4:25])
all_Q3_covid_11 <- cbind(all_Q3_covid_11[,c(1:3)],all_Q3_covid_11_M)
all_Q3_covid_11$Time <- "2020_Q3"
colnames(all_Q3_covid_11)[4:25] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                     "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                     "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                     "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
all_Q3_exp_11 <- merge(all_Q3_covid_11,all_Q3_expNoC_11[,-c(2)],"LocID")
all_Q3_exp_11_M <- as.matrix(all_Q3_exp_11[,27:48]) -
  (as.matrix(all_Q3_exp_11[4:25])/8000)
all_Q3_exp_11 <- cbind(all_Q3_covid_11[,c(1:3)],all_Q3_exp_11_M)
#
all_Q3_other_11 <- merge(all_20_nmx_11, all_Q3_exp_11[,-c(2:3)],"LocID")
all_Q3_other_11_M <- as.matrix(all_Q3_other_11[,c(3:24)]) *
  as.matrix(all_Q3_other_11[,c(25:46)])
all_Q3_other_11 <- cbind(all_Q3_other_11[,c(1:2)],all_Q3_other_11_M)
#
all_Q3_death_11 <- merge(all_Q3_other_11,all_Q3_covid_11[,-c(2:3)],"LocID")
all_Q3_death_11_M <- (as.matrix(all_Q3_death_11[,25:46])/1000) +
  as.matrix(all_Q3_death_11[,3:24])
all_Q3_death_11 <- cbind(all_Q3_death_11[,1:2],all_Q3_death_11_M)
#
all_Q3_dXS_11 <- merge(all_Q3_death_11,all_Q3_dNoC_11[,-c(2)],"LocID")
all_Q3_dXS_11_M <- as.matrix(all_Q3_dXS_11[,3:24]) - 
  as.matrix(all_Q3_dXS_11[,25:46])
all_Q3_dXS_11 <- cbind(all_Q3_dXS_11[,1:2],all_Q3_dXS_11_M)
#
all_Q3_chratio_11 <- merge(all_Q3_death_11,all_Q3_dNoC_11[,-c(2)],"LocID")
all_Q3_chratio_11_M <- as.matrix(all_Q3_chratio_11[,3:24]) / 
  as.matrix(all_Q3_chratio_11[,25:46])
all_Q3_chratio_11 <- cbind(all_Q3_chratio_11[,1:2],all_Q3_chratio_11_M)
colnames(all_Q3_chratio_11)[3:24] <- c("0R1.M","4R1.M","10R5.M","10R15.M","10R25.M",
                                       "10R35.M","10R45.M","10R55.M","10R65.M","10R75.M","R85+.M",
                                       "0R1.F","4R1.F","10R5.F","10R15.F","10R25.F","10R35.F",
                                       "10R45.F","10R55.F","10R65.F","1R75.F","R85+.F")
#
###
#
all_Q4_expNoC_11 <- merge(all_Q3_dXS_11,all_Q3_expNoC_11[,-c(2)],"LocID")
all_Q4_expNoC_11_M <- as.matrix(all_Q4_expNoC_11[,25:46]) -
  (as.matrix(all_Q4_expNoC_11[3:24]) /4)
all_Q4_expNoC_11 <- cbind(all_Q4_expNoC_11[,c(1:2)],all_Q4_expNoC_11_M)
#
all_Q4_dNoC_11 <- merge(all_20_nmx_11,all_Q4_expNoC_11[,-c(2)],"LocID")
all_Q4_dNoC_11_M <- as.matrix(all_Q4_dNoC_11[,c(3:24)])*
  as.matrix(all_Q4_dNoC_11[,25:46])
all_Q4_dNoC_11 <- cbind(all_Q4_dNoC_11[,1:2],all_Q4_dNoC_11_M)
colnames(all_Q4_dNoC_11)[3:24] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                    "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                    "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                    "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
#
all_Q4_covid_11 <- merge(all_4Q_covid[,c(1,2,10)],all_Q4_expNoC_11[,-c(2)],"LocID")
#
#multiplying by covid asdr
#
all_Q4_covid_11_M <- as.matrix(all_Q4_covid_11[,4:25]) %*% usa_asdr_11_M 
#
all_Q4_covid_11 <- cbind(all_Q4_covid_11[,c(1:3)],all_Q4_covid_11_M,na.rm = TRUE)
all_Q4_ratio <-  as.vector(all_Q4_covid_11$Covid_Q4) / 
  as.vector(rowSums(all_Q4_covid_11_M))
##
all_Q4_ratio_M <- diag(all_Q4_ratio,339,339)
##
all_Q4_covid_11_M <- as.matrix(all_Q4_ratio_M) %*% 
  as.matrix(all_Q4_covid_11[4:25])
all_Q4_covid_11 <- cbind(all_Q4_covid_11[,c(1:3)],all_Q4_covid_11_M)
all_Q4_covid_11$Time <- "2020_Q4"
colnames(all_Q4_covid_11)[4:25] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                     "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                     "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                     "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
all_Q4_exp_11 <- merge(all_Q4_covid_11,all_Q4_expNoC_11[,-c(2)],"LocID")
all_Q4_exp_11_M <- as.matrix(all_Q4_exp_11[,27:48]) -
  (as.matrix(all_Q4_exp_11[4:25])/8000)
all_Q4_exp_11 <- cbind(all_Q4_covid_11[,c(1:3)],all_Q4_exp_11_M)
#
all_Q4_other_11 <- merge(all_20_nmx_11, all_Q4_exp_11[,-c(2:3)],"LocID")
all_Q4_other_11_M <- as.matrix(all_Q4_other_11[,c(3:24)]) *
  as.matrix(all_Q4_other_11[,c(25:46)])
all_Q4_other_11 <- cbind(all_Q4_other_11[,c(1:2)],all_Q4_other_11_M)
#
all_Q4_death_11 <- merge(all_Q4_other_11,all_Q4_covid_11[,-c(2:3)],"LocID")
all_Q4_death_11_M <- (as.matrix(all_Q4_death_11[,25:46])/1000) +
  as.matrix(all_Q4_death_11[,3:24])
all_Q4_death_11 <- cbind(all_Q4_death_11[,1:2],all_Q4_death_11_M)
##
all_Q4_dXS_11 <- merge(all_Q4_death_11,all_Q4_dNoC_11[,-c(2)],"LocID")
all_Q4_dXS_11_M <- as.matrix(all_Q4_dXS_11[,3:24]) - 
  as.matrix(all_Q4_dXS_11[,25:46])
all_Q4_dXS_11 <- cbind(all_Q4_dXS_11[,1:2],all_Q4_dXS_11_M)
#
all_Q4_chratio_11 <- merge(all_Q4_death_11,all_Q4_dNoC_11[,-c(2)],"LocID")
all_Q4_chratio_11_M <- as.matrix(all_Q4_chratio_11[,3:24]) / 
  as.matrix(all_Q4_chratio_11[,25:46])
all_Q4_chratio_11 <- cbind(all_Q4_chratio_11[,1:2],all_Q4_chratio_11_M)
colnames(all_Q4_chratio_11)[3:24] <- c("0R1.M","4R1.M","10R5.M","10R15.M","10R25.M",
                                       "10R35.M","10R45.M","10R55.M","10R65.M","10R75.M","R85+.M",
                                       "0R1.F","4R1.F","10R5.F","10R15.F","10R25.F","10R35.F",
                                       "10R45.F","10R55.F","10R65.F","1R75.F","R85+.F")
#
###
#
all_Q1_expNoC_11 <- merge(all_Q4_dXS_11,all_Q4_expNoC_11[,-c(2)],"LocID")
all_Q1_expNoC_11_M <- as.matrix(all_Q1_expNoC_11[,25:46]) -
                         (as.matrix(all_Q1_expNoC_11[3:24]) /4)
all_Q1_expNoC_11 <- cbind(all_Q1_expNoC_11[,c(1:2)],all_Q1_expNoC_11_M)
#
all_Q1_dNoC_11 <- merge(all_20_nmx_11,all_Q1_expNoC_11[,-c(2)],"LocID")
all_Q1_dNoC_11_M <- as.matrix(all_Q1_dNoC_11[,c(3:24)])*
  as.matrix(all_Q1_dNoC_11[,25:46])
all_Q1_dNoC_11 <- cbind(all_Q1_dNoC_11[,1:2],all_Q1_dNoC_11_M)
colnames(all_Q1_dNoC_11)[3:24] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                    "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                    "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                    "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
#
all_Q1_covid_11 <- merge(all_4Q_covid[,c(1,2,11)],all_Q1_expNoC_11[,-c(2)],"LocID")
#
#multiplying by covid asdr
#
all_Q1_covid_11_M <- as.matrix(all_Q1_covid_11[,4:25]) %*% usa_asdr_11_M 
#
all_Q1_covid_11 <- cbind(all_Q1_covid_11[,c(1:3)],all_Q1_covid_11_M,na.rm = TRUE)
all_Q1_ratio <-  as.vector(all_Q1_covid_11$Covid_Q1) / 
  as.vector(rowSums(all_Q1_covid_11_M))
##
all_Q1_ratio_M <- diag(all_Q1_ratio,339,339)
##
all_Q1_covid_11_M <- as.matrix(all_Q1_ratio_M) %*% 
  as.matrix(all_Q1_covid_11[4:25])
all_Q1_covid_11 <- cbind(all_Q1_covid_11[,c(1:3)],all_Q1_covid_11_M)
all_Q1_covid_11$Time <- "2021_Q1"
colnames(all_Q1_covid_11)[4:25] <- c("0d1.M","4d1.M","10d5.M","10d15.M","10d25.M",
                                     "10d35.M","10d45.M","10d55.M","10d65.M","10d75.M","d85+.M",
                                     "0d1.F","4d1.F","10d5.F","10d15.F","10d25.F","10d35.F",
                                     "10d45.F","10d55.F","10d65.F","10d75.F","d85+.F")
##
all_Q1_exp_11 <- merge(all_Q1_covid_11,all_Q1_expNoC_11[,-c(2)],"LocID")
all_Q1_exp_11_M <- as.matrix(all_Q1_exp_11[,27:48]) -
  (as.matrix(all_Q1_exp_11[4:25])/8000)
all_Q1_exp_11 <- cbind(all_Q1_covid_11[,c(1:3)],all_Q1_exp_11_M)
#
all_Q1_other_11 <- merge(all_20_nmx_11, all_Q1_exp_11[,-c(2:3)],"LocID")
all_Q1_other_11_M <- as.matrix(all_Q1_other_11[,c(3:24)]) *
  as.matrix(all_Q1_other_11[,c(25:46)])
all_Q1_other_11 <- cbind(all_Q1_other_11[,c(1:2)],all_Q1_other_11_M)
#
all_Q1_death_11 <- merge(all_Q1_other_11,all_Q1_covid_11[,-c(2:3)],"LocID")
all_Q1_death_11_M <- (as.matrix(all_Q1_death_11[,25:46])/1000) +
  as.matrix(all_Q1_death_11[,3:24])
all_Q1_death_11 <- cbind(all_Q1_death_11[,1:2],all_Q1_death_11_M)
##
all_Q1_dXS_11 <- merge(all_Q1_death_11,all_Q1_dNoC_11[,-c(2)],"LocID")
all_Q1_dXS_11_M <- as.matrix(all_Q1_dXS_11[,3:24]) - 
  as.matrix(all_Q1_dXS_11[,25:46])
all_Q1_dXS_11 <- cbind(all_Q1_dXS_11[,1:2],all_Q1_dXS_11_M)
#
all_Q1_chratio_11 <- merge(all_Q1_death_11,all_Q4_dNoC_11[,-c(2)],"LocID")
all_Q1_chratio_11_M <- as.matrix(all_Q1_chratio_11[,3:24]) / 
  as.matrix(all_Q1_chratio_11[,25:46])
all_Q1_chratio_11 <- cbind(all_Q1_chratio_11[,1:2],all_Q1_chratio_11_M)
colnames(all_Q1_chratio_11)[3:24] <- c("0R1.M","4R1.M","10R5.M","10R15.M","10R25.M",
                                       "10R35.M","10R45.M","10R55.M","10R65.M","10R75.M","R85+.M",
                                       "0R1.F","4R1.F","10R5.F","10R15.F","10R25.F","10R35.F",
                                       "10R45.F","10R55.F","10R65.F","1R75.F","R85+.F")
#
##
##step 4: apply Chiang's formula to calculate new values of the npx 
##
all_20_noCnpx_10 <- reshape(all_20_lt_11[,c(2:3,5,17:26)], 
        timevar = "Sex",
        idvar = c("LocID", "Location"),
        direction = "wide")
#
all_Q2_noCnpx_10 <- merge(all_20_noCnpx_10,all_Q2_chratio_11[,-c(2,13,24)],
                          "LocID")
all_Q2_npx_10_M <- as.matrix(all_Q2_noCnpx_10[,3:22]) ^ 
  as.matrix(all_Q2_noCnpx_10[,23:42])
all_Q2_npx_10 <- cbind(all_Q2_noCnpx_10[,1:2],all_Q2_npx_10_M)
##
all_Q3_noCnpx_10 <- merge(all_20_noCnpx_10,all_Q3_chratio_11[,-c(2,13,24)],
                          "LocID")
all_Q3_npx_10_M <- as.matrix(all_Q3_noCnpx_10[,3:22]) ^ 
  as.matrix(all_Q3_noCnpx_10[,23:42])
all_Q3_npx_10 <- cbind(all_Q3_noCnpx_10[,1:2],all_Q3_npx_10_M)
#
all_Q4_noCnpx_10 <- merge(all_20_noCnpx_10,all_Q4_chratio_11[,-c(2,13,24)],
                            "LocID")
all_Q4_npx_10_M <- as.matrix(all_Q4_noCnpx_10[,3:22]) ^ 
  as.matrix(all_Q4_noCnpx_10[,23:42])
all_Q4_npx_10 <- cbind(all_Q4_noCnpx_10[,1:2],all_Q4_npx_10_M)
#
all_Q1_noCnpx_10 <- merge(all_20_noCnpx_10,all_Q1_chratio_11[,-c(2,13,24)],
                          "LocID")
all_Q1_npx_10_M <- as.matrix(all_Q1_noCnpx_10[,3:22]) ^ 
  as.matrix(all_Q1_noCnpx_10[,23:42])
all_Q1_npx_10 <- cbind(all_Q1_noCnpx_10[,1:2],all_Q1_npx_10_M)
##
##step 3.4: calculate new values of nax
##
#first derive previous values of nax from npx and nmx in 2020 life table
#
#second part of formula first
#
all_20_noCnax_10_M <- (10 * as.matrix(all_20_noCnpx_10[,3:22])) / 
  (1 - as.matrix(all_20_noCnpx_10[,3:22]))
all_20_noCnax_10 <- cbind(all_20_noCnpx_10[,1:2],all_20_noCnax_10_M)
#
#need to modify value of n in first 2 age groups
#
all_20_noCnax_10[,3] <- all_20_noCnax_10[,3]/10
all_20_noCnax_10[,4] <- all_20_noCnax_10[,4]*4/10
all_20_noCnax_10[,13] <- all_20_noCnax_10[,13]/10
all_20_noCnax_10[,14] <- all_20_noCnax_10[,14]*4/10
#
#subtract from first part of formula
#
all_20_noCnax_11 <- merge(all_20_noCnax_10,all_20_nmx_11[,-c(2)],"LocID")
all_20_noCnax_11_M <- (1/as.matrix(all_20_noCnax_11[,c(23:32,34:43)])) - 
  as.matrix(all_20_noCnax_11[,3:22])
#
#different treatment for open age-interval
#
all_20_noCnax_op_M <- (1/as.matrix(all_20_noCnax_11[,c(33,44)]))
all_20_noCnax_11 <- cbind(all_20_noCnax_11[,1:2],all_20_noCnax_11_M)
all_20_noCnax_11 <- cbind(all_20_noCnax_11,all_20_noCnax_op_M)
#
#note the open age groups male and female are at the end
#from old values of nax to new ones
#add Chiang ratios
#
all_Q2_nax_11 <- merge(all_20_noCnax_11,all_Q2_chratio_11[,-c(2)],"LocID")
all_Q3_nax_11 <- merge(all_20_noCnax_11,all_Q3_chratio_11[,-c(2)],"LocID")
all_Q4_nax_11 <- merge(all_20_noCnax_11,all_Q4_chratio_11[,-c(2)],"LocID")
all_Q1_nax_11 <- merge(all_20_noCnax_11,all_Q1_chratio_11[,-c(2)],"LocID")
#
#subtract n from nax
#
all_Q2_nax_11[,3:22] <- all_Q2_nax_11[,3:22]-10
all_Q3_nax_11[,3:22] <- all_Q3_nax_11[,3:22]-10
all_Q4_nax_11[,3:22] <- all_Q4_nax_11[,3:22]-10
all_Q1_nax_11[,3:22] <- all_Q1_nax_11[,3:22]-10
#
#different values of n for the first two age groups
#
all_Q2_nax_11[,3] <- all_Q2_nax_11[,3]+9
all_Q2_nax_11[,4] <- all_Q2_nax_11[,4]+6
all_Q2_nax_11[,13] <- all_Q2_nax_11[,13]+9
all_Q2_nax_11[,14] <- all_Q2_nax_11[,14]+6
#
all_Q3_nax_11[,3] <- all_Q3_nax_11[,3]+9
all_Q3_nax_11[,4] <- all_Q3_nax_11[,4]+6
all_Q3_nax_11[,13] <- all_Q3_nax_11[,13]+9
all_Q3_nax_11[,14] <- all_Q3_nax_11[,14]+6
#
all_Q4_nax_11[,3] <- all_Q4_nax_11[,3]+9
all_Q4_nax_11[,4] <- all_Q4_nax_11[,4]+6
all_Q4_nax_11[,13] <- all_Q4_nax_11[,13]+9
all_Q4_nax_11[,14] <- all_Q4_nax_11[,14]+6
#
all_Q1_nax_11[,3] <- all_Q1_nax_11[,3]+9
all_Q1_nax_11[,4] <- all_Q1_nax_11[,4]+6
all_Q1_nax_11[,13] <- all_Q1_nax_11[,13]+9
all_Q1_nax_11[,14] <- all_Q1_nax_11[,14]+6
#
#multiply the difference by Chiang ratios
#
all_Q2_nax_10_M <- as.matrix(all_Q2_nax_11[,3:22]) * 
  as.matrix(all_Q2_nax_11[,c(25:34,36:45)])
#
all_Q3_nax_10_M <- as.matrix(all_Q3_nax_11[,3:22]) * 
  as.matrix(all_Q3_nax_11[,c(25:34,36:45)])
#
all_Q4_nax_10_M <- as.matrix(all_Q4_nax_11[,3:22]) * 
  as.matrix(all_Q4_nax_11[,c(25:34,36:45)])
#
all_Q1_nax_10_M <- as.matrix(all_Q1_nax_11[,3:22]) * 
  as.matrix(all_Q1_nax_11[,c(25:34,36:45)])
#
#get the ratio of the old to new nqx
#
all_Q2_qratio_10 <-merge(all_20_noCnpx_10[,1:22],all_Q2_npx_10[,-c(2)],"LocID")
all_Q2_qratio_10_M <- (1-as.matrix(all_Q2_qratio_10[,3:22]))/
  (1-as.matrix(all_Q2_qratio_10[,23:42]))
#
all_Q3_qratio_10 <-merge(all_20_noCnpx_10[,1:22],all_Q3_npx_10[,-c(2)],"LocID")
all_Q3_qratio_10_M <- (1-as.matrix(all_Q3_qratio_10[,3:22]))/
  (1-as.matrix(all_Q3_qratio_10[,23:42]))
#
all_Q4_qratio_10 <-merge(all_20_noCnpx_10[,1:22],all_Q4_npx_10[,-c(2)],"LocID")
all_Q4_qratio_10_M <- (1-as.matrix(all_Q4_qratio_10[,3:22]))/
                           (1-as.matrix(all_Q4_qratio_10[,23:42]))
#
all_Q1_qratio_10 <-merge(all_20_noCnpx_10[,1:22],all_Q1_npx_10[,-c(2)],"LocID")
all_Q1_qratio_10_M <- (1-as.matrix(all_Q1_qratio_10[,3:22]))/
  (1-as.matrix(all_Q1_qratio_10[,23:42]))
#
#multiply the two products
#
all_Q2_nax_10_M <- as.matrix(all_Q2_nax_10_M) * as.matrix(all_Q2_qratio_10_M)
all_Q3_nax_10_M <- as.matrix(all_Q3_nax_10_M) * as.matrix(all_Q3_qratio_10_M)
all_Q4_nax_10_M <- as.matrix(all_Q4_nax_10_M) * as.matrix(all_Q4_qratio_10_M)
all_Q1_nax_10_M <- as.matrix(all_Q1_nax_10_M) * as.matrix(all_Q1_qratio_10_M)
#
#combine with values needed for open age groups
#
all_Q2_nax_11 <- cbind(all_Q2_nax_11[,c(1:2,23:24,35,46)],all_Q2_nax_10_M)
all_Q3_nax_11 <- cbind(all_Q3_nax_11[,c(1:2,23:24,35,46)],all_Q3_nax_10_M)
all_Q4_nax_11 <- cbind(all_Q4_nax_11[,c(1:2,23:24,35,46)],all_Q4_nax_10_M)
all_Q1_nax_11 <- cbind(all_Q1_nax_11[,c(1:2,23:24,35,46)],all_Q1_nax_10_M)
#
#add n to resulting values for closed age-groups
#
all_Q2_nax_11[,7:26] <- all_Q2_nax_11[,7:26]+10
all_Q3_nax_11[,7:26] <- all_Q3_nax_11[,7:26]+10
all_Q4_nax_11[,7:26] <- all_Q4_nax_11[,7:26]+10
all_Q1_nax_11[,7:26] <- all_Q1_nax_11[,7:26]+10
#
#different values of n for the first two age groups
#
all_Q2_nax_11[,7] <- all_Q2_nax_11[,7]-9
all_Q2_nax_11[,8] <- all_Q2_nax_11[,8]-6
all_Q2_nax_11[,17] <- all_Q2_nax_11[,17]-9
all_Q2_nax_11[,18] <- all_Q2_nax_11[,18]-6
#
all_Q3_nax_11[,7] <- all_Q3_nax_11[,7]-9
all_Q3_nax_11[,8] <- all_Q3_nax_11[,8]-6
all_Q3_nax_11[,17] <- all_Q3_nax_11[,17]-9
all_Q3_nax_11[,18] <- all_Q3_nax_11[,18]-6
#
all_Q4_nax_11[,7] <- all_Q4_nax_11[,7]-9
all_Q4_nax_11[,8] <- all_Q4_nax_11[,8]-6
all_Q4_nax_11[,17] <- all_Q4_nax_11[,17]-9
all_Q4_nax_11[,18] <- all_Q4_nax_11[,18]-6
#
all_Q1_nax_11[,7] <- all_Q1_nax_11[,7]-9
all_Q1_nax_11[,8] <- all_Q1_nax_11[,8]-6
all_Q1_nax_11[,17] <- all_Q1_nax_11[,17]-9
all_Q1_nax_11[,18] <- all_Q1_nax_11[,18]-6
#
#for open-ended age groups divide previous value by Chiang ratio
#
all_Q2_nax_11[,27] <- all_Q2_nax_11[,3]/all_Q2_nax_11[,5]
all_Q2_nax_11[,28] <- all_Q2_nax_11[,4]/all_Q2_nax_11[,6]
#
all_Q3_nax_11[,27] <- all_Q3_nax_11[,3]/all_Q3_nax_11[,5]
all_Q3_nax_11[,28] <- all_Q3_nax_11[,4]/all_Q3_nax_11[,6]
#
all_Q4_nax_11[,27] <- all_Q4_nax_11[,3]/all_Q4_nax_11[,5]
all_Q4_nax_11[,28] <- all_Q4_nax_11[,4]/all_Q4_nax_11[,6]
#
all_Q1_nax_11[,27] <- all_Q1_nax_11[,3]/all_Q1_nax_11[,5]
all_Q1_nax_11[,28] <- all_Q1_nax_11[,4]/all_Q1_nax_11[,6]
#
#only keep necessary values and reorder
#
all_Q2_nax_11 <- all_Q2_nax_11[,c(1:2,7:16,27,17:26,28)]
colnames(all_Q2_nax_11)[3:24] <- c("0a1.M","4a1.M","10a5.M","10a15.M","10a25.M",
                                   "10a35.M","10a45.M","10a55.M","10a65.M","10a75.M","a85+.M",
                                   "0a1.F","4a1.F","10a5.F","10a15.F","10a25.F","10a35.F",
                                   "10a45.F","10a55.F","10a65.F","1a75.F","a85+.F")
all_Q3_nax_11 <- all_Q3_nax_11[,c(1:2,7:16,27,17:26,28)]
colnames(all_Q3_nax_11)[3:24] <- colnames(all_Q2_nax_11)[3:24]
all_Q4_nax_11 <- all_Q4_nax_11[,c(1:2,7:16,27,17:26,28)]
colnames(all_Q4_nax_11)[3:24] <- colnames(all_Q2_nax_11)[3:24]
all_Q1_nax_11 <- all_Q1_nax_11[,c(1:2,7:16,27,17:26,28)]
colnames(all_Q2_nax_11)[3:24] <- colnames(all_Q2_nax_11)[3:24]
#
###
##step 3.5: calculate original values of life expectancies from npx and nmx
#
#chain calculation from oldest age groups but format of life table file
#allows to do for both sexes in single loop
#
all_20_lt_11$e85 <- 1/all_20_lt_11[,16]
for (i in 1:10) {all_20_lt_11[,27+i] <- (all_20_lt_11[,27-i]*all_20_lt_11[,26+i])+
  ((1-all_20_lt_11[,27-i])/all_20_lt_11[,16-i])}
colnames(all_20_lt_11)[28:37] <- c("e75","e65","e55","e45","e35","e25","e15","e5",
                                     "e1","e0")
#
all_20_noCex_11 <- reshape(all_20_lt_11[,c(2:3,5,27:37)], 
                           timevar = "Sex",
                           idvar = c("LocID", "Location"),
                           direction = "wide")
#
##
###step 4: calculate Mean Unfulfilled length of Life/Lifespan/Lifetime/Longevity 
##in each quarter and semester
##
##for deaths in age group x to x+n, it is:
##the difference between the original and new average ages at death nax 
##plis, for excess deaths, addo the difference between the original life expectantcy  
##at age x and the original average age at death nax
#
##for excess deaths in the open age group x+ it is:
##the original life expectancy at age x
#
##first prepare the files with these differences
##
all_20_mul_21 <- merge(all_20_noCex_11[,c(1:3,14,13,12,11,10,9,8,7,6,5,4,
            24,23,22,21,20,19,18,17,16,15)],all_20_noCnax_11[,-c(2,23:24)],"LocID")
all_20_mul_10_M <- as.matrix(all_20_mul_21[,c(5:24)])-
  as.matrix(all_20_mul_21[,c(25:44)])
all_20_mul_21 <- cbind(all_20_mul_21[,c(1:4,25:44)],all_20_mul_10_M)
#
all_Q2_mul_21 <- merge(all_20_mul_21,all_Q2_nax_11[,-c(2,13,24)],"LocID")
all_Q2_mul_10_M <- as.matrix(all_Q2_mul_21[,c(5:24)])-
  as.matrix(all_Q2_mul_21[,c(45:64)])
all_Q2_mul_21 <- cbind(all_Q2_mul_21[,c(1:4,25:44)],all_Q2_mul_10_M)
##
#arrange the number of deaths in each corresponding groups
##
all_Q2_death_21 <- merge(all_Q2_dXS_11[,c(1:2,13,24,3:12,14:23)],
                         all_Q2_death_11[,-c(2,13,24)],"LocID")
#
#multiply the difference in each sex-and-age group 
#by number of deaths in that group
#
all_Q2_tul_21 <- merge(all_Q2_mul_21,all_Q2_death_21[,-c(2)],"LocID")
all_Q2_tul_21_M <- as.matrix(all_Q2_tul_21[,c(3:44)]) *
  as.matrix(all_Q2_tul_21[,c(45:86)])
all_Q2_tul_21 <- cbind(all_Q2_tul_21[,c(1:2)],all_Q2_tul_21_M)
#
#calculate the total and mean unfulfilled lifespan in Q1 death cohort
#
all_Q2_mul_21$Q2_TUL_M <- rowSums(all_Q2_tul_21[,c(3,5:14,25:34)],na.rm=TRUE)
all_Q2_mul_21$Q2_TUL_F <- rowSums(all_Q2_tul_21[,c(4,15:24,35:44)],na.rm=TRUE)
all_Q2_mul_21$Q2_death_M <- rowSums(all_Q2_death_11[,c(3:13)],na.rm=TRUE)
all_Q2_mul_21$Q2_death_F <- rowSums(all_Q2_death_11[,c(14:24)],na.rm=TRUE)
all_Q2_mul_21$Q2_MUL_M <- all_Q2_mul_21$Q2_TUL_M/all_Q2_mul_21$Q2_death_M
all_Q2_mul_21$Q2_MUL_F <- all_Q2_mul_21$Q2_TUL_F/all_Q2_mul_21$Q2_death_F
all_Q2_mul_21$Q2_MUL_B <- (all_Q2_mul_21$Q2_TUL_M+all_Q2_mul_21$Q2_TUL_F)/
  (all_Q2_mul_21$Q2_death_M+all_Q2_mul_21$Q2_death_F)
##
##also provide an average per covid death
#
all_4Q_mul_1 <- merge(all_4Q_covid[,-c(3:6)],all_Q2_mul_21[,-c(2:44)],
                      "LocID")
all_4Q_mul_1$Q2_PAYLL <- (all_4Q_mul_1$Q2_TUL_M+all_4Q_mul_1$Q2_TUL_F)*1000/
  all_4Q_mul_1$Covid_Q2
###
#
##
#repeat for Q3
#
all_Q3_mul_21 <- merge(all_20_mul_21,all_Q3_nax_11[,-c(2,13,24)],"LocID")
all_Q3_mul_10_M <- as.matrix(all_Q3_mul_21[,c(5:24)])-
  as.matrix(all_Q3_mul_21[,c(45:64)])
all_Q3_mul_21 <- cbind(all_Q3_mul_21[,c(1:4,25:44)],all_Q3_mul_10_M)
##
all_Q3_death_21 <- merge(all_Q3_dXS_11[,c(1:2,13,24,3:12,14:23)],
                         all_Q3_death_11[,-c(2,13,24)],"LocID")
#
all_Q3_tul_21 <- merge(all_Q3_mul_21,all_Q3_death_21[,-c(2)],"LocID")
all_Q3_tul_21_M <- as.matrix(all_Q3_tul_21[,c(3:44)]) *
  as.matrix(all_Q3_tul_21[,c(45:86)])
all_Q3_tul_21 <- cbind(all_Q3_tul_21[,c(1:2)],all_Q3_tul_21_M)
#
all_Q3_mul_21$Q3_TUL_M <- rowSums(all_Q3_tul_21[,c(3,5:14,25:34)],na.rm=TRUE)
all_Q3_mul_21$Q3_TUL_F <- rowSums(all_Q3_tul_21[,c(4,15:24,35:44)],na.rm=TRUE)
all_Q3_mul_21$Q3_death_M <- rowSums(all_Q3_death_11[,c(3:13)],na.rm=TRUE)
all_Q3_mul_21$Q3_death_F <- rowSums(all_Q3_death_11[,c(14:24)],na.rm=TRUE)
all_Q3_mul_21$Q3_MUL_M <- all_Q3_mul_21$Q3_TUL_M/all_Q3_mul_21$Q3_death_M
all_Q3_mul_21$Q3_MUL_F <- all_Q3_mul_21$Q3_TUL_F/all_Q3_mul_21$Q3_death_F
all_Q3_mul_21$Q3_MUL_B <- (all_Q3_mul_21$Q3_TUL_M+all_Q3_mul_21$Q3_TUL_F)/
  (all_Q3_mul_21$Q3_death_M+all_Q3_mul_21$Q3_death_F)
##
all_4Q_mul_1 <- merge(all_4Q_mul_1,all_Q3_mul_21[,-c(2:44)],
                      "LocID")
all_4Q_mul_1$Q3_PAYLL <- (all_4Q_mul_1$Q3_TUL_M+all_4Q_mul_1$Q3_TUL_F)*1000/
  all_4Q_mul_1$Covid_Q3
###
##
#repeat for Q4
#
all_Q4_mul_21 <- merge(all_20_mul_21,all_Q4_nax_11[,-c(2,13,24)],"LocID")
all_Q4_mul_10_M <- as.matrix(all_Q4_mul_21[,c(5:24)])-
  as.matrix(all_Q4_mul_21[,c(45:64)])
all_Q4_mul_21 <- cbind(all_Q4_mul_21[,c(1:4,25:44)],all_Q4_mul_10_M)
##
all_Q4_death_21 <- merge(all_Q4_dXS_11[,c(1:2,13,24,3:12,14:23)],
                         all_Q4_death_11[,-c(2,13,24)],"LocID")
#
all_Q4_tul_21 <- merge(all_Q4_mul_21,all_Q4_death_21[,-c(2)],"LocID")
all_Q4_tul_21_M <- as.matrix(all_Q4_tul_21[,c(3:44)]) *
  as.matrix(all_Q4_tul_21[,c(45:86)])
all_Q4_tul_21 <- cbind(all_Q4_tul_21[,c(1:2)],all_Q4_tul_21_M)
#
all_Q4_mul_21$Q4_TUL_M <- rowSums(all_Q4_tul_21[,c(3,5:14,25:34)],na.rm=TRUE)
all_Q4_mul_21$Q4_TUL_F <- rowSums(all_Q4_tul_21[,c(4,15:24,35:44)],na.rm=TRUE)
all_Q4_mul_21$Q4_death_M <- rowSums(all_Q4_death_11[,c(3:13)],na.rm=TRUE)
all_Q4_mul_21$Q4_death_F <- rowSums(all_Q4_death_11[,c(14:24)],na.rm=TRUE)
all_Q4_mul_21$Q4_MUL_M <- all_Q4_mul_21$Q4_TUL_M/all_Q4_mul_21$Q4_death_M
all_Q4_mul_21$Q4_MUL_F <- all_Q4_mul_21$Q4_TUL_F/all_Q4_mul_21$Q4_death_F
all_Q4_mul_21$Q4_MUL_B <- (all_Q4_mul_21$Q4_TUL_M+all_Q4_mul_21$Q4_TUL_F)/
  (all_Q4_mul_21$Q4_death_M+all_Q4_mul_21$Q4_death_F)
##
all_4Q_mul_1 <- merge(all_4Q_mul_1,all_Q4_mul_21[,-c(2:44)],
                      "LocID")
all_4Q_mul_1$Q4_PAYLL <- (all_4Q_mul_1$Q4_TUL_M+all_4Q_mul_1$Q4_TUL_F)*1000/
  all_4Q_mul_1$Covid_Q4
###
#repeat for Q1
#
all_Q1_mul_21 <- merge(all_20_mul_21,all_Q1_nax_11[,-c(2,13,24)],"LocID")
all_Q1_mul_10_M <- as.matrix(all_Q1_mul_21[,c(5:24)])-
  as.matrix(all_Q1_mul_21[,c(45:64)])
all_Q1_mul_21 <- cbind(all_Q1_mul_21[,c(1:4,25:44)],all_Q1_mul_10_M)
##
all_Q1_death_21 <- merge(all_Q1_dXS_11[,c(1:2,13,24,3:12,14:23)],
                         all_Q1_death_11[,-c(2,13,24)],"LocID")
#
all_Q1_tul_21 <- merge(all_Q1_mul_21,all_Q1_death_21[,-c(2)],"LocID")
all_Q1_tul_21_M <- as.matrix(all_Q1_tul_21[,c(3:44)]) *
  as.matrix(all_Q1_tul_21[,c(45:86)])
all_Q1_tul_21 <- cbind(all_Q1_tul_21[,c(1:2)],all_Q1_tul_21_M)
#
all_Q1_mul_21$Q1_TUL_M <- rowSums(all_Q1_tul_21[,c(3,5:14,25:34)],na.rm=TRUE)
all_Q1_mul_21$Q1_TUL_F <- rowSums(all_Q1_tul_21[,c(4,15:24,35:44)],na.rm=TRUE)
all_Q1_mul_21$Q1_death_M <- rowSums(all_Q1_death_11[,c(3:13)],na.rm=TRUE)
all_Q1_mul_21$Q1_death_F <- rowSums(all_Q1_death_11[,c(14:24)],na.rm=TRUE)
all_Q1_mul_21$Q1_MUL_M <- all_Q1_mul_21$Q1_TUL_M/all_Q1_mul_21$Q1_death_M
all_Q1_mul_21$Q1_MUL_F <- all_Q1_mul_21$Q1_TUL_F/all_Q1_mul_21$Q1_death_F
all_Q1_mul_21$Q1_MUL_B <- (all_Q1_mul_21$Q1_TUL_M+all_Q1_mul_21$Q1_TUL_F)/
  (all_Q1_mul_21$Q1_death_M+all_Q1_mul_21$Q1_death_F)
##
all_4Q_mul_1 <- merge(all_4Q_mul_1,all_Q1_mul_21[,-c(2:44)],
                      "LocID")
all_4Q_mul_1$Q1_PAYLL <- (all_4Q_mul_1$Q1_TUL_M+all_4Q_mul_1$Q1_TUL_F)*1000/
  all_4Q_mul_1$Covid_Q1
###
#also calculate 12-month values
#
all_4Q_mul_1$L12m_TUL_M <- all_4Q_mul_1$Q2_TUL_M+
  all_4Q_mul_1$Q3_TUL_M+all_4Q_mul_1$Q4_TUL_M+all_4Q_mul_1$Q1_TUL_M
all_4Q_mul_1$L12m_TUL_F <- all_4Q_mul_1$Q2_TUL_F+
  all_4Q_mul_1$Q3_TUL_F+all_4Q_mul_1$Q4_TUL_F+all_4Q_mul_1$Q1_TUL_F
all_4Q_mul_1$L12m_TUL_B <- all_4Q_mul_1$L12m_TUL_F+all_4Q_mul_1$L12m_TUL_M
all_4Q_mul_1$L12m_MUL_M <- all_4Q_mul_1$L12m_TUL_M/(all_4Q_mul_1$Q2_death_M+
  all_4Q_mul_1$Q3_death_M+all_4Q_mul_1$Q4_death_M+all_4Q_mul_1$Q1_death_M)
all_4Q_mul_1$L12m_MUL_F <- all_4Q_mul_1$L12m_TUL_F/(all_4Q_mul_1$Q2_death_F+
  all_4Q_mul_1$Q3_death_F+all_4Q_mul_1$Q4_death_F+all_4Q_mul_1$Q1_death_F)
all_4Q_mul_1$L12m_MUL_B <- all_4Q_mul_1$L12m_TUL_B/(all_4Q_mul_1$Q2_death_F+
  all_4Q_mul_1$Q3_death_F+all_4Q_mul_1$Q4_death_F+all_4Q_mul_1$Q1_death_F+
    all_4Q_mul_1$Q2_death_M+all_4Q_mul_1$Q3_death_M+all_4Q_mul_1$Q4_death_M+
    all_4Q_mul_1$Q1_death_M)
#
#add differences in life expectancies for the year
#
all_l12m_diff <- read.csv('output_210331.csv')
colnames(all_l12m_diff)[1] <- c("LocID")
all_4Q_mort <- merge(all_4Q_mul_1[,c(1:2,14,22,30,38,45,15,23,31,39,12:13,20:21,
   28:29,36:37,43:44,8:9,16:17,24:25,32:33,40:41,4:7)],all_l12m_diff[,c(1,19:21)],
   "LocID")
#
colnames(all_4Q_mort)[c(32:38)] <- c("Q2_death_CoV","Q3_death_Cov",
  "Q4_death_Cov","Q1_death_CoV","L12m_Diff_M","L12m_Diff_F","L12m_Diff_B")
#
#sort on average loss for both sexes over the entire 12-month period
#
all_4Q_mort <- all_4Q_mort[order(all_4Q_mort$"L12m_MUL_B", na.last=TRUE,
                                 decreasing=TRUE), ]
#
#to match quarterly value, make 12-month differences positive
#
all_4Q_mort$`L12m_Diff_M` <- 0-all_4Q_mort$`L12m_Diff_M`
all_4Q_mort$`L12m_Diff_F` <- 0-all_4Q_mort$`L12m_Diff_F`
all_4Q_mort$`L12m_Diff_B` <- 0-all_4Q_mort$`L12m_Diff_B`
#
write.csv(all_4Q_mort,file="all_4Q_mort.csv")
#
#add a summary file with fewer results
all_4Q_mort_summ <- all_4Q_mort[,c(1:7,38,8:11,32:35,12:21,36:37)]
all_4Q_mort_summ <- all_4Q_mort_summ[order(all_4Q_mort_summ$"LocID", na.last=TRUE,
                                 decreasing=FALSE), ]
#
write.csv(all_4Q_mort_summ,file="all_4Q_mort_summ.csv")
#
####
#####