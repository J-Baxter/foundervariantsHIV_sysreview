###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for univariate IPD meta-regression under a one-step approach
# One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies

###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(mltools)
library(ggplot2)
library(influence.ME)
library(kableExtra)
library(parallel)
library(performance)
library(reshape2)
library(cowplot)
library(stringr)
library(data.table)
library(meta)
source('~/foundervariantsHIV_sysreview/generalpurpose_funcs.R')


# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula, opt = NULL){
  
  if(is.character(opt)){
    cntrl <- glmerControl(optCtrl = list(maxfun = 5000000),
                          check.nobs.vs.nlev = 'ignore',
                          check.nobs.vs.nRE = 'ignore',
                          optimizer = opt)
  }else{
    cntrl <- glmerControl(optCtrl = list(maxfun = 5000000),
                          check.nobs.vs.nlev = 'ignore',
                          check.nobs.vs.nRE = 'ignore')
  }
  
  options(warn = 1)
  f <- as.formula(formula)
  environment(f) <- environment()
  model <- lme4::glmer(f,
                       data = data,
                       family = binomial(link = "logit"),
                       nAGQ = 7,
                       control = cntrl)
  return(model)
}


CalcUniVar <- function(data,...){
  props <- CalcProps(data,...)
  var <- colnames(props)[2]
  # metaprop wraps rma.glmm which in turn wraps lme4 glmm
  # advantage of using specific meta analysis packages is ease of calculating
  # heterogeneity statistics
  require(meta)
  model <- meta::metaprop(
    event = props[,4], 
    n = props[,3], 
    data = props, 
    studlab = props[,1],
    byvar = props[,var], 
    comb.fixed = F, 
    method.tau = 'ML', 
    method = 'GLMM' , 
    addincr = 0.0005, 
    allincr = F,
    hakn = T, 
    backtransf = F)
  return(model)
}


GetEstimates <- function(model){
  df = cbind.data.frame(
    covar = colnames(model$data)[2] %>% gsub('[_]', '', .),
    level = model$bylevs,
    k = model$k.w,
    n = model$n.w,
    beta = model$TE.random.w,
    beta.lb = model$lower.random.w,
    beta.ub = model$upper.random.w,
    tau2 = model$tau2.w,
    i2 = model$I2.w,
    Q = model$Q.w,
    pval.Q = model$pval.Q.w
  )
  qtest <- cbind.data.frame(
    covar = colnames(model$data)[2] %>% gsub('[_]', '', .),
    Q = model$Q.b.random, 
    df = model$df.Q.b, 
    pval.Q = model$pval.Q.b.random)
  
  out <- list(qtest,df)
  return(out)
}
###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% 
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  droplevels()


##################################################################################################
###################################################################################################
# STAGE 2: Univariate meta-regression of individual covariates against founder variant multiplicity
# Initial regression models with one fixed effect covariate with random effects for publication and cohort
# Equivalent to a subgroup analysis with random effects for subgroup and cohort

fixeff_uni.forms <- c(f0 = "multiple.founders_ ~  1 + (1 | publication_)",
                      f1 = "multiple.founders_ ~  riskgroup_ - 1 + (1 | publication_)",
                      f2 = "multiple.founders_ ~ reported.exposure_ - 1 +  ( 1 | publication_)",
                      f3 = "multiple.founders_ ~ grouped.method_ - 1 +  ( 1 | publication_)",
                      f5 = "multiple.founders_ ~ sequencing.gene_ - 1 + (1 | publication_)",
                      f6 = "multiple.founders_ ~ sampling.delay_ - 1 +  ( 1 | publication_)")

fixeff_uni.effectstruct <- GetName(fixeff_uni.forms, effects = 'fixed')

# Run models
fixeff_uni.models <- RunParallel(CalcRandMetaReg, fixeff_uni.forms, df)

# Check model convergence and singularity
lapply(fixeff_uni.models, check_singularity)
lapply(fixeff_uni.models, check_convergence)

###################################################################################################
# Test for differences between levels within subgroups


###################################################################################################
# Tabulate outputs


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################