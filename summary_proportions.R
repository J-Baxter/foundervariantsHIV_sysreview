# IPD meta analysis of HIV founder variant multiplicity
# calculates summary proportion of infections initiated by multiple founder variants
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. stratified intercepts and random effects
#    for between study heterogeneity, approx ML fit
# 3. One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
#    (uncorrelated intercept and slope). approx ML fit
# 4. Two-step beta-binomial GLMM, dispersion param for study labels. Laplace approximate ML estimation

# estimations of mean effect size, confidence intervals and heterogeneity are presented 


#dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(mltools)
library(data.table)
library(metafor)
library(dmetar)
library(aod)

# Formats data spreadsheet for analysis. Removes duplicates and NAs.
formatDF <-  function(df){
  #create dummy variables in founder multiplicity col
  if (class(df$multiple.founders)=="factor"){
    df$multiple.founders = 1 - (as.numeric(df$multiple.founders)-1)
  }else{
    stop('founder multiplicity is not a factor')
  }
  df_nona <- df[!is.na(df$multiple.founders),]
  df_nodups <- df_nona[(df_nona$include.main == '') & (df_nona$exclude.repeatstudy == ''),]
  df_labelled <- unite(df_nodups, "publication", c(author ,year), sep = '_')
  return(df_labelled)
}


# Sums number of patients within each study and number of infections initiated by multiple founders. 
# Can be stratified with additional covariates
CalcProps <- function(.data, ...){
    .data %>% 
      group_by(publication, ...) %>%
      summarise(subjects = n(), multiplefounders = sum(multiple.founders))
  }


# Creates dummy variables for clustering random effects of studies
onehotEncode <- function(data, covar, names){
  onehot <- cbind.data.frame(data[, covar]) %>%
    data.table::as.data.table() %>%
    mltools::one_hot()
  
  colnames(onehot) <- names[order(names)]
  
  testset_onehot <- cbind.data.frame(data,onehot)
  
  return(testset_onehot)
}


# First step glm of two-step binomial/normal model
BNstepOne <- function(data, study_id){
  df_subset <- subset(data, publication==study_id)
  events <- nrow(df_subset)
  success <-  sum(df_subset$multiple.founders)
  prop.multifounders <- success/events
  incr = 0.0005 
  incr.event = incr/events
  
  if (prop.multifounders == 0){
    lr <- glm(multiple.founders+incr.event ~ 1, df_subset , family = binomial(link = "logit"))
    log_or <- summary(lr)$coefficients[,1]
    se <- summary(lr)$coefficients[, 2]
    
  }else{
    lr <- glm(multiple.founders ~ 1, df_subset , family = binomial(link = "logit"))
    log_or <- summary(lr)$coefficients[,1]
    se <- summary(lr)$coefficients[, 2]
  }
  
  agg.results <- cbind.data.frame(study_id, log_or, se)
  
  return(agg.results)
}


###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()
  
# Set test data
publist <- unique(df$publication)
testlist <- c("Keele_2008", "Abrahams_2009", "Haaland_2009","Li_2010", "Janes_2015")#
testset_df <- lapply(testlist, function(x,y) subset(x, publication == y), x = df) %>% do.call(rbind.data.frame,.)

#set seed
set.seed(4472)


###################################################################################################

# 1. Two-step binomial-normal model
# step 1: pooling within studies using binomial model/logit 
# step 2: pooling across studies using random effects (normal model). Inverse variance method used to pool. 
# REML estimator of Tau. 
step1_bn <- lapply(publist, BNstepOne , data = df) %>% do.call(rbind.data.frame,.)

step2_bn <- rma.uni(log_or, se, 
                    data = step1_withinstudy,
                    method = "REML",
                    knha = TRUE, 
                    measure = "PLO")

twostep.sum <- summary(step2_bn)
twostep.sum


###################################################################################################

# 2. One-step binomial GLMM allowing for clustering by study. 
# random slope of x within group with correlated intercept
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution

df_onestage <- onehotEncode(df, covar = "publication", names = publist)

onestep_bi_strat <- glmer(multiple.founders ~  1 +  (1 |publication),
                           data = df,
                           family = binomial(link = "logit"),
                           control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

onestep_bi_strat.sum <- summary(onestep_bi_strat)
onestep_bi_strat.sum

onestep_bi_strat.tau2 <- VarCorr(onestep_bi_strat)[[1]][1] 
onestep_bi_strat.tau2_se <- 

###################################################################################################

# 3. One-step binomial GLMM allowing for clustering by study. 
# uncorrelated random intercept and random slope within group 
# approx ML fit
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution
# (1 | random.factor) + (0 + fixed.factor | random.factor) = fixed.factor + (fixed.factor || random.factor)

onestep_bi_rand <- glmer(multiple.founders ~  1 + (1|publication) + (0+1|publication)  ,
                        data = df,
                        family = binomial(link = "logit"),
                        control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

onestep_bi_rand.sum <- summary(onestep_bi_rand)
onestep_bi_rand.sum

onestep_bi_rand.tau2 <- VarCorr(onestep_bi_rand)[[1]][1] 

###################################################################################################

# 4. Beta-binomial model
# the event rate pi for the ith study is drawn from a beta distirbution. conditional on pi, the number of 
# individuals xi in the ith study of size ni who experience the events of interest follows a binomial 
# distribution B(ni;pi) (From Chuang-Stein 1993).
# Laplace approximate ML estimation

df_props <- CalcProps(df, reported.exposure)
betabinom <- betabin(cbind(multiplefounders, subjects - multiplefounders) ~ 1, ~ 1, data = df_props)

betabinom.sum <- summary(betabinom)
betabinom.sum

binom.ci <- varbin(subjects,multiplefounders, data = df_props)@tab[5,c(3,4)] #Bootstrapped Binomial CIs-check Chuang-Stein 1993
###################################################################################################

# Model comparison: Estimated sumary effects (prop, CI), between study variance (tau, I^)
summary_props <- c(step2_bn$beta,
                   onestep_bi_strat.sum$coefficients[1,1],
                   onestep_bi_ind.sum$coefficients[1,1], 
                   betabinom@param[1]) %>% 
  as.numeric() %>% 
  transf.ilogit()
  
summary_props.ci95 <- list(c(step2_bn$ci.lb , step2_bn$ci.ub),
                           c(confint(onestep_bi_strat)[2,c(1,2)]),
                           c(confint(onestep_bi_rand)[2,c(1,2)]),
                           binom.ci) %>% lapply(.,transf.ilogit)

tau2 <- c(step2_bn$tau2, onestep_bi_strat.tau2, onestep_bi_rand.tau2 )

tau2.ci <- 
  
I2 <- 
  
I2.ci <- 
  

modelcomp_df <- cbind.data.frame(c('meta.ran.reml.hk', 'onestep_bn', 'onestep_bb'),
                                 summary_props, sapply(summary_props.ci95, rbind.data.frame) %>% t())

colnames(modelcomp_df) <- c('model', 'proportions', 'props.ci95_lower', 'props.ci95_upper')

comp.table <- as_tibble(modelcomp_df)


# Model Comparison: Does beta-binomial better represent overdispersion than binomial mode
AIC(onestep_bi,onestep_bb)
lrtest(onestep_bi,onestep_bb)

# Influence analyses

 

