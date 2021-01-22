# IPD meta analysis of HIV founder variant multiplicity
# calculates summary proportion of infections initiated by multiple founder variants
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. stratified intercepts and random effects
#    for between study heterogeneity, approx ML fit
# 3. One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
#    (uncorrelated intercept and slope). approx ML fit
# 4. One-step beta-binomial GLMM, dispersion param for study labels. Laplace approximate ML estimation

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

# 2. One-step binomial GLMM allowing for clustering by study. stratified intercepts and random effects for between 
#    study heterogeneity, approx ML fit 
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution
df_onestage <- onehotEncode(df, covar = "publication", names = publist)

onestep_bi_strat <- glmer(multiple.founders ~ 1 + Brooks_2020 + Leda_2020 + Liu_2020 + Macharia_2020 + 
                            Martinez_2020 + Rolland_2020 + VillabonaArenas_2020 + Sivay_2019 + 
                            Todesco_2019 + Tovanabutra_2019 + Ashokkumar_2018 + Dukhovlinova_2018 + 
                            deCamp_2017 + Kijak_2017 + Lyer_2017  + AbigailSmith_2016 + Chaillon_2016 + 
                            Novitsky_2016 +Oberle_2016 + SalazarGonzalez_2016 + Tully_2016 + 
                            Chen_2015 + Danaviah_2015 + Deymier_2015 + Gounder_2015 + Janes_2015 + 
                            Le_2015 + Zanini_2015 + Chaillon_2014 + Sterrett_2014  + Wagner_2014 +
                            Baalwa_2013 + Frange_2013 + Henn_2012 + Kiwelu_2012 + Rossenkhan_2012 + 
                            Sturdevant_2012 + CollinsFairclough_2011 + Cornelissen_2011 + Herbeck_2011 +
                            Kishko_2011 + Nofemela_2011 + Novitsky_2011 + Rachinger_2011 + Rieder_2011 + 
                            Rolland_2011 + Bar_2010 + Fischer_2010 + Li_2010 + Masharsky_2010 + 
                            Zhang_2010 + Abrahams_2009 + Haaland_2009 + Kearney_2009 + Novitsky_2009 + 
                            SalazarGonzalez_2009 + Gottlieb_2008 + Keele_2008 + Kwiek_2008 + 
                            SalazarGonzalez_2008 + Sagar_2006 + Derdeyn_2004 + Ritola_2004 + Sagar_2004 + 
                            Renjifo_2003 + Sagar_2003 + Verhofstede_2003 +  Delwart_2002 + Learn_2002 + 
                            Long_2002 + Nowak_2002 + Dickover_2001 + Long_2000 +Wade_1998 + Briant_1995 + 
                            Poss_1995 + Wolinsky_1992 + (1 | publication),
                          data = df_onestage,
                          family = binomial,
                          control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))

# known warning 1: fixed-effect model matrix is rank deficient so dropping 1 column / coefficient boundary 
# (singular) fit: see ?isSingular
# known warning 2: maxfun < 10 * length(par)^2 is not recommended.

onestep_bi_strat.sum <- summary(onestep_bi_strat)
# known warning: In vcov.merMod(object, correlation = correlation, sigm = sig) :
# variance-covariance matrix computed from finite-difference Hessian is
# not positive definite or contains NA values: falling back to var-cov estimated from RX
# onestep_bi_strat.sum


###################################################################################################

# 3. One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies 
#    (uncorrelated intercept and slope). approx ML fit
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution
# (1 | random.factor) + (0 + fixed.factor | random.factor) = fixed.factor + (fixed.factor || random.factor)

onestep_bi_ind <- glmer(multiple.founders ~ (1|publication) + ((0+1|publication)) ,
                        data = df,
                        family = binomial,
                        control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))

onestep_bi_ind.sum <- summary(onestep_bi_ind)
onestep_bi_ind.sum


###################################################################################################

# 4. Beta-binomial model
# the event rate pi for the ith study is drawn from a beta distirbution. conditional on pi, the number of 
# individuals xi in the ith study of size ni who experience the events of interest follows a binomial 
# distribution B(ni;pi) (From Chuang-Stein 1993).
# Laplace approximate ML estimation

#current issue NAs for all values other than estimate.
df_props <- CalcProps(df)
onestep_bb <- betabin(cbind(multiplefounders, subjects - multiplefounders) ~ 1, ~ 1, data = df_props)

onestep_bb.sum <- summary(onestep_bb)
onestep_bb.sum


###################################################################################################

# Model comparison: Estimated sumary effects (prop, CI), between study variance (tau, I^)
summary_props <- c(step2_bn$beta,
                   onestep_bi_strat.sum$coefficients[1,1],
                   onestep_bi_ind.sum$coefficients[1,1], 
                   onestep_bb.sum$coefficients$cond[1,1]) %>% transf.ilogit()
  
summary_props.ci95 <- list(c(meta.ran.reml.hk$lower.random, meta.ran.reml.hk$upper.random),
                             c(confint(onestep_bi)[2,c(1,2)]),
                             c(confint(onestep_bb, component = "cond")[1,c(1,2)])) %>% lapply(.,transf.ilogit)


modelcomp_df <- cbind.data.frame(c('meta.ran.reml.hk', 'onestep_bn', 'onestep_bb'),
                                 summary_props, sapply(summary_props.ci95, rbind.data.frame) %>% t())

colnames(modelcomp_df) <- c('model', 'proportions', 'props.ci95_lower', 'props.ci95_upper')

comp.table <- as_tibble(modelcomp_df)


# Model Comparison: Does beta-binomial better represent overdispersion than binomial mode
AIC(onestep_bi,onestep_bb)
lrtest(onestep_bi,onestep_bb)

# Influence analyses

 

