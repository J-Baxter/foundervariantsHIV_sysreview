###################################################################################################
###################################################################################################
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
# Sensitivity analyses are conducted as follows:
# SA1: leave-one-out cross validation of all models
# SA2: exclude studies with fewer than 10 participants
# SA3: excluding studies that only report single founder infectioins (ie prop multifounder =0)
# SA4: alternative inclusion criteria for repeated participants
# SA5: bootstrap/resampling of participants for whom multiple measurements are available

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
library(reshape2)
source('generalpurpose_funcs.R')


# Add increment to studies with only single founder infections (for one-step models)
AddIncr <- function(df, incr){
  split_df <- split.data.frame(df , df$publication)
  
  list_incr = list()
  for (i in 1:length(split_df)){
    data = split_df[[i]]
    if(!(1 %in% data$multiple.founders)){
      data$multiple.founders = incr
    }else{
      data = data}
    list_incr[[i]] = data
  }
  df_incr = do.call(rbind.data.frame,list_incr)
  stopifnot(nrow(df) == nrow(df_incr))
  
  return(df_incr)
}


# Two-step binomial/normal model, pooling studies using Inverse Variance method,
# random effects, REML estimator of tau2.
CalcTwostepBiNorm <- function(data){
  step1 <- escalc(xi = multiplefounders , ni = subjects , data= data , add = 0.0005, measure = "PLO")
  step2 <- rma.uni(yi, vi, 
                   data = step1,
                   method = "REML",
                   knha = TRUE, 
                   measure = "PLO")
  
  list(step1, step2) %>%
    return()
}


# One-step GLMM accounting for clustering of studies using a statified intercept (
# random slope, correlated intercept)
CalcOnestepBiStrat <- function(data){
  model <- glmer(multiple.founders ~  1 + factor(publication) + (1| publication),
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optCtrl = list(maxfun = 1000000)))
  return(model)
}


# One-step GLMM accounting for clustering of studies using a random intercept (
# random slope, random & uncorrelated intercept)
CalcOnestepBiRand <- function(data){
  model <- lme4::glmer(multiple.founders ~  1 + (1|publication),
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  return(model)
}


# Two-step (Compound dist) GLMM accounting for clustering of studies using dispersion param Phi
# of beta distriubtion to capture between study heterogeneity of p
CalcTwostepBetaBi <- function(proportions){
  model <- betabin(cbind(multiplefounders, subjects - multiplefounders) ~ 1, ~ 1, data = proportions)
  return(model)
}


CalcCI <- function(u,se,threshold){
  value <- 1-(threshold/2)
  upper <- u + se*qnorm(value)
  lower <- u - se*qnorm(value)
  ci <- c(lower,upper)
  return(ci)
}

# Extracts estimates of summary effect from models
# CURRENTLY FULLY FUNCTIONAL ONLY FOR METAFOR MODELS  
CalcEstimates <- function(model , analysis = "original"){
  if (class(model) == "rma" || class(model) == "rma.uni"){
    beta <- model$beta
    ci.lb <- model$ci.lb
    ci.ub <- model$ci.ub
  }
  else if (class(model) =="glmerMod"){
    beta <- summary(model)$coefficients[1,1]
    ci <- confint(model)
    ci.lb <- ci[nrow(ci),1]
    ci.ub <- ci[nrow(ci),2]
  }
  else if (class(model) =="glimML"){
    beta <- model@param[1]
    ci <- varbin(subjects,multiplefounders, data = model@data)@tab[5,c(3,4)]
    ci.lb <- as.numeric(ci[1]) %>% transf.logit()
    ci.ub <- as.numeric(ci[2]) %>% transf.logit()
  }

  results <- cbind.data.frame("estimate" = beta,
                              "estimate.lb" = ci.lb,
                              "estimate.ub" = ci.ub) %>%
    mapply(transf.ilogit, ., SIMPLIFY = FALSE)
  
  results_lab <- cbind.data.frame("model" = substitute(model) %>% deparse(),
                                  "analysis" = analysis,
                                  results)
  
  return(results_lab)
}


# Extracts estimates of heterogeneity (tau2) from models
CalcTau2 <- function(model , analysis = "original"){
  if (class(model) == "rma" || class(model) == "rma.uni"){
    tau2 <- model$tau2
    tau2.lb <- model$tau2.lb
    tau2.ub <- model$tau2.ub
  }
  else if (class(model) =="glmerMod"){
    tau2 <- VarCorr(model)[[1]][1] 
    tau2.lb <- 
    tau2.ub <- c
  }
  else if (class(model) =="glimML"){
    tau2 <- 
    tau2.lb <- 
    tau2.ub <- 
  }
  
  results <- cbind.data.frame("model" = substitute(model) %>% deparse(),
                              "tau2" = tau2,
                              "tau2.lb" = tau2.lb,
                              "tau2.ub" = tau2.ub,
                              "analysis" = analysis)
  
  return(results)
}


# Create list of dataframes for leave-one-out cross validation
LOOCV.dat <- function(data){
  pubs <- unique(data$publication)
  loo <- list()
  loo.pubs <- list()
  
  for (i in pubs){
    loo[[i]] <- data[data$publication != i, ]
    loo.pubs[[i]] <- pubs[pubs != i]
  }
  out <- list(loo, loo.pubs)
  stopifnot(length(loo) == length(loo.pubs))
  return(out)
}


# Extract estimates from LOOCV to create dataframe (input for influence plot)
DFInfluence <- function(model,labs){
  
  names <- paste("Omitting" , labs %>% names(), sep = " ") %>% as.factor()
  beta = numeric()
  ci.lb = numeric()
  ci.ub = numeric() 
  
  if (class(model[[1]]) == "rma" || class(model[[1]]) == "rma.uni"){
    for (i in 1:length(model)){
      beta[i] <- model[[i]]$beta
      ci.lb[i] <-model[[i]]$ci.lb
      ci.ub[i] <- model[[i]]$ci.ub
    }
  }else if (class(model[[1]]) =="glmerMod"){
    for (i in 1:length(model)){
      ci <- confint(model[[i]])
      beta[i] <- summary(model[[i]])$coefficients[1,1]
      ci.lb[i] <-ci[nrow(ci),1]
      ci.ub[i] <- ci[nrow(ci),2]
    }
  }else if (class(model[[1]]) =="glimML"){
    #Bootstrapped Binomial CIs-check Chuang-Stein 1993
    for (i in 1:length(model)){
      binom.ci <- varbin(subjects,multiplefounders, data = model[[i]]@data)@tab[5,c(3,4)] %>%
        as.numeric() %>%
        transf.logit()
      beta[i] <- model[[i]]@param[1]
      ci.lb[i] <- binom.ci[1]
      ci.ub[i] <- binom.ci[2]
    }
  }else{
    stop('no valid model detected.')
  }
  
  influence_out <- cbind.data.frame('trial'= names,
                                    "estimate" = transf.ilogit(beta),
                                    "ci.lb" = transf.ilogit(ci.lb),
                                    "ci.ub"= transf.ilogit(ci.ub)) 
  

  return(influence_out)
}

# Generate resampled datasets and calculate model estimates for psuedo-bootstrap 
# sensitivity analysis of inclusion/exclusion criteria
BootParticipant <- function(data, replicates){
  require(parallel)
  require(lme4)
  require(metafor)
  require(dplyr)
  
  resampled <- lapply(1:replicates, function(x,y) {y %>% group_by(participant.id) %>% slice_sample(n=1)},
                      y = data)
  
  resampled_props <- lapply(resampled , CalcProps)
  
  cl <- detectCores() %>%
    `-` (2) %>%
    makeCluster()
  
  clusterEvalQ(cl, c(library(lme4), library(metafor), library(aod), library(dplyr),
                     set.seed(4472),'CalcOnestepBiRand', 'CalcTwostepBiNorm',
                     'CalcOnestepBiStrat','CalcTwostepBetaBi'))
  
  start <- Sys.time()
  print(start)
  
  twostep_boot <- parLapply(cl = cl, resampled_props, CalcTwostepBiNorm)
  rand_boot <- parLapply(cl = cl, resampled, CalcOnestepBiRand)
  strat_boot <- parLapply(cl = cl, resampled, CalcOnestepBiStrat)
  beta_boot <- parLapply(cl = cl, resampled_props, CalcTwostepBetaBi)
  
  end <- Sys.time()
  elapsed <- end-start
  print(elapsed)
  
  stopCluster(cl)
  remove(cl)
  
  twostep_boot.est <- lapply(twostep_boot, function(model) model[[2]]$beta) %>% 
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  strat_boot.est <- lapply(strat_boot, function(model) summary(model)$coefficients[1,1]) %>% 
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  rand_boot.est <- lapply(rand_boot, function(model) summary(model)$coefficients[1,1]) %>% 
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  beta_boot.est <- lapply(beta_boot, function(model) model@param[1]) %>% 
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  boot_estimates <- list('twostep' = twostep_boot.est,
                         'rand' = rand_boot.est,
                         'strat' = strat_boot.est,
                         'beta' = beta_boot.est)
  
  return(boot_estimates)
}


###################################################################################################
###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF(., noreps = TRUE)
df_props <- CalcProps(df)  

# Set seed
set.seed(4472)


###################################################################################################

# 1. Two-step binomial-normal model
# step 1: pooling within studies using binomial model/logit 
# step 2: pooling across studies using random effects (normal model). Inverse variance method used to pool. 
# REML estimator of Tau. 

twostep_binorm.all <- CalcTwostepBiNorm(df_props)
twostep_binorm.step1 <- twostep_binorm[[1]]
twostep_binorm <- twostep_binorm[[2]]

twostep_binorm.sum <- summary(twostep_binorm.step2)
twostep_binorm.sum

twostep_binorm.est <- CalcEstimates(twostep_binorm)
###################################################################################################

# 2. One-step binomial GLMM allowing for clustering by study. 
# random slope of x within group with correlated intercept
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution

df_onestage <- AddIncr(df, incr = 0.0005)

onestep_bi_strat <- CalcOnestepBiStrat(df_onestage )
onestep_bi_strat.sum <- summary(onestep_bi_strat)
onestep_bi_strat.sum

onestep_bi_rand.est <- CalcEstimates(onestep_bi_strat)
onestep_bi_strat.tau2 <- 



###################################################################################################

# 3. One-step binomial GLMM allowing for clustering by study. 
# uncorrelated random intercept and random slope within group 
# approx ML fit
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution
# (1 | random.factor) + (0 + fixed.factor | random.factor) = fixed.factor + (fixed.factor || random.factor)

onestep_bi_rand <- CalcOnestepBiRand(df_onestage) 
onestep_bi_rand.sum <- summary(onestep_bi_rand)
onestep_bi_rand.sum

onestep_bi_rand.est <- CalcEstimates(onestep_bi_rand)
onestep_bi_rand.tau2 <- 

###################################################################################################

# 4. Beta-binomial model
# the event rate pi for the ith study is drawn from a beta distirbution. conditional on pi, the number of 
# individuals xi in the ith study of size ni who experience the events of interest follows a binomial 
# distribution B(ni;pi) (From Chuang-Stein 1993).
# Laplace approximate ML estimation


df_props$multiplefounders[df_props$multiplefounders == 0 ] <- 0.0005 # refactor to CalcProps function

twostep_betabi <- CalcTwostepBetaBi(df_props)

twostep_betabi.sum <- summary(twostep_betabi)
twostep_betabi.sum
twostep_betabi.est <- CalcEstimates(twostep_betabi)
 
###################################################################################################
###################################################################################################
# Model comparison: Estimated sumary effects (prop, CI), between study variance (tau, I^)
# Tau2 = Var(theta_i), theta_i = E[theta_i]
estimates <- rbind.data.frame(twostep_binorm.est, onestep_bi_rand.est,
                                    onestep_bi_rand.est,twostep_betabi.est)


tau2 <- rbind.data.frame(twostep_binorm.tau2, onestep_bi_rand.tau2,
                         onestep_bi_rand.tau2,twostep_betabi.tau2 )

  
###################################################################################################
###################################################################################################
# Sensitivity Analyses
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)


# SA1. Influence of Individual Studies (LOOCV)
df_loocv <- LOOCV.dat(df)[[1]]
publist_loocv <- LOOCV.dat(df)[[2]]
dfp_loocv <- LOOCV.dat(df_props)[[1]]

twostep_binorm.influence <- lapply(dfp_loocv, function(x) CalcTwostepBiNorm(data = x)[[2]]) %>%
  DFInfluence(., labs = publist_loocv) %>% {cbind.data.frame(.,'model' = 'twostep_binorm')}

onestep_bi_strat.influence <- lapply(df_loo ,CalcOnestepBiStrat) %>%
  DFInfluence()

onestep_bi_rand.influence <-  lapply(df_loocv ,CalcOnestepBiRand) %>%
  DFInfluence(., labs =publist_loocv) %>% {cbind.data.frame(.,'model' = 'onestep_rand')}

twostep_betabi.influence <-  lapply(dfp_loocv ,CalcTwostepBetaBi) %>%
  DFInfluence(., labs = publist_loocv) %>% {cbind.data.frame(.,'model' = 'twostep_betabi')}

influence_df <- rbind.data.frame(twostep_binorm.influence,onestep_bi_rand.influence,twostep_betabi.influence)


# SA2. Exclusion of small sample sizes (less than n = 10)
publist.nosmallsample <- subset(df_props , subjects > 9 , select = publication) %>%
  pull(.,var=publication) %>%
  unique()

df.nosmallsample <- lapply(publist.nosmallsample, function(x,y) subset(x, publication == y), x = df) %>%
  do.call(rbind.data.frame,.)

df_props.nosmallsample <- subset(df_props , subjects > 9)

twostep_binorm.nosamllsample <- CalcTwostepBiNorm(df_props.nosmallsample)[[2]] %>% 
  CalcEstimates(., analysis = "no_small")

onestep_bi_strat.nosamllsample <- CalcOnestepBiStrat(df.nosmallsample) %>%
  CalcEstimates(., analysis = "no_small")

onestep_bi_rand.nosamllsample <- CalcOnestepBiRand(df.nosmallsample) %>% 
  CalcEstimates(., analysis = "no_small")

twostep_betabi.nosamllsample <- CalcTwostepBetaBi(df_props.nosmallsample) %>%
  CalcEstimates(., analysis = "no_small")

SA2_results <-  rbind.data.frame(twostep_binorm.nosamllsample, onestep_bi_rand.nosamllsample,
                                 onestep_bi_rand.nosamllsample,twostep_betabi.nosamllsample)


# SA3. Exclusion of studies with 0 multiple founder variants
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication) %>%
  pull(.,var=publication) %>%
  unique()

df.nozeros <- lapply(publist.nozeros, function(x,y) subset(x, publication == y), x = df) %>%
  do.call(rbind.data.frame,.)

df_props.nozeros <- subset(df_props , multiplefounders != 0)

twostep_binorm.nozeros <- CalcTwostepBiNorm(df_props.nozeros)[[2]] %>%
  CalcEstimates(., analysis = "no_zeros")

onestep_bi_strat.nozeros <- CalcOnestepBiStrat(df.nozeros) %>%
  CalcEstimates(., analysis = "no_zeros")

onestep_bi_rand.nozeros <- CalcOnestepBiRand(df.nozeros) %>%
  CalcEstimates(., analysis = "no_zeros")

twostep_betabi.nozeros <- CalcTwostepBetaBi(df_props.nozeros) %>%
  CalcEstimates(., analysis = "no_zeros") #expectation is this is no better than binomial models

SA3_results <- rbind.data.frame(twostep_binorm.nozeros, onestep_bi_rand.nozeros,
                                onestep_bi_rand.nozeros,twostep_betabi.nozeros)


# SA4. Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)
resampling_df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF(., noreps = FALSE)

boot_participant <- BootParticipant(resampling_df , 500) #Out to file: save as df -> csv with indexes denoting list structure

###################################################################################################
###################################################################################################
# Outputs
# CSV of pooling model results (estimates only) with SAs 2 + 3
pooled_mods <- rbind.data.frame(estimates,
                                SA2_results,
                                SA3_results)

write.csv(pooled_mods, file = 'bp_estsa2sa3.csv')

# CSV of pooling model results (estimates and tau2)
model_comp_df <- cbind.data.frame(estimates, tau2)

write.csv(model_comp_df , file = 'bp_esttau.csv')

# CSV study influence
write.csv(influence_df, file = 'bp_sa1.csv')


# Forest Plot 2-step BN
pdf("testplot.pdf" , width = 14 , height = 20)
forest(twostep_binorm.step2 , showweights = TRUE, slab = sort(publist) , transf = transf.ilogit, header = TRUE, digits = 3,
       refline = 0.283, ilab = cbind(df_props$multiplefounders, df_props$subjects), ilab.xpos = c(-0.9, -0.4), psize = 1, cex = 0.75,
       )
text(c(-0.9,-0.4), 79, c('Multiple Founders' , 'Subjects') , font = 2)

dev.off()

###################################################################################################
###################################################################################################
# END #
###################################################################################################
###################################################################################################