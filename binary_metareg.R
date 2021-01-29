###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for IPD meta-regression under one-step and two-step approaches
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. stratified intercepts and random effects
#    for between study heterogeneity, approx ML fit
# OR. One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
#    (uncorrelated intercept and slope). approx ML fit
# OR. Two-step beta-binomial GLMM, dispersion param for study labels. Laplace approximate ML estimation

###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(mltools)
library(data.table)
library(metafor)
library(dmetar)
library(aod)
library(ggplot2)
library(influence.ME)
library(kableExtra)
source('generalpurpose_funcs.R')

# One-step GLMM accounting for clustering of studies using a statified intercept (
# random slope, correlated intercept)
CalcOnestepBiStrat <- function(data){
  model <- glmer(multiple.founders ~  factor(publication) + (1| publication),
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  return(model)
}


###################################################################################################
###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()

set.seed(4472)

###################################################################################################
