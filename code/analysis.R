# final analysis for the paper

set.seed(1234) # for reproducible results with the bayesian prior
library(stringr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(DescTools)
library(epiR)
library(dampack)

# Data preparation
source("~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/code/data_prep.R")

# 1. Population characteristics 
# Table 1:
referral.data.raw[! is.na(referral.data.raw$baseline.group),] %>%
  group_by(baseline.group, age.group ) %>% summarise(n=n())

referral.data.raw[! is.na(referral.data.raw$baseline.group),] %>%
  group_by(baseline.group, genotype.data ) %>% summarise(n=n())

referral.data.raw[! is.na(referral.data.raw$baseline.group),] %>%
group_by(hist.groups, first.round.hist, baseline.group ) %>% summarise(n=n()) %>% print(n=51)

round(c(66, 1, 109, 22, 1, 21, 157)/sum(c(66, 1, 109, 22, 1, 21, 157))*100, 1)
round(c(30, 1, 84, 40, 1, 60, 510)/sum(c(30, 1, 84, 40, 1, 60, 510))*100, 1)
round(c(2, 0, 3, 7, 0, 12, 189)/sum(c(2, 0, 3, 7, 0, 12, 189))*100, 1)

# 2. Colposcopy referrals in the first and second screening round
#   - Figure 2 and Table S1
source("~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/code/referral_rates.R") 

# 3. PPV analysis
#   -  Table 2: PPV for CIN3+ of 14 different strategies in hrHPV-positive women stratified by baseline cytology.
source("~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/code/PPV_analysis.R") 

# 4. NPV analysis
#   - Table 3: NPV for CIN3+ of 14 different strategies in hrHPV-positive women.
source("~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/code/NPV_analysis.R") 

# 5. Incremental analysis 
#   - Figure 3: Efficient frontier of the fourteen referral strategies. 
source("~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/code/INNR_analysis.R") 
