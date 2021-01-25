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
  df_labelled <- unite(df_nodups, "publication", c(author ,year), sep = ' ')
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


# First step glm of two-step binomial/normal model. Zero value cells + 0.0005
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


# Second step of two-step binomial/normal model, pooling studies using Inverse Variance method,
# random effects, REML estimator of tau2.
CalcTwostepBiNorm <- function(data, study_list){
  step1 <- lapply(study_list, BNstepOne , data = data) %>% do.call(rbind.data.frame,.)
  step2 <- rma.uni(log_or, se, 
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
  model <- glmer(multiple.founders ~  1 +  ( 1 |publication),
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  return(model)
}


# One-step GLMM accounting for clustering of studies using a random intercept (
# random slope, random & uncorrelated intercept)
CalcOnestepBiRand <- function(data){
  model <- glmer(multiple.founders ~  1 + (1|publication) + (0 + 1|publication),
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


#extracts estimates of summary effect from models
#can potentially refactor around DFInfluence
CalcEstimates <- function(model1, model2, model3, model4){
  summary_estimates <- c(model1$beta,
                         summary(model2)$coefficients[1,1],
                         summary(model3)$coefficients[1,1], 
                         model4@param[1]) %>% 
    as.numeric() %>% 
    transf.ilogit() %>% cbind.data.frame()
  
  binom.ci <- varbin(subjects,multiplefounders, data = df_props)@tab[5,c(3,4)] 
  #Bootstrapped Binomial CIs-check Chuang-Stein 1993
  
  summary_props.ci95 <- list(c(model1$ci.lb , model1$ci.ub),
                             c(confint(model2)[2,c(1,2)]),
                             c(confint(model3)[2,c(1,2)]),
                             binom.ci) %>% 
    lapply(.,as.numeric) %>% 
    lapply(.,transf.ilogit) %>% do.call(rbind.data.frame,.)
  
  results <- cbind.data.frame(summary_props,
                              summary_props.ci95[,1],
                              summary_props.ci95[,2])
  
  colnames(results) <- c("Estimate", "95% CI Lower", "95% CI Upper" )
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


#extract estimates from LOOCV to create dataframe (input for influence plot)
DFInfluence <- function(model){
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
      beta[i] <- summary(model[[i]])$coefficients[1,1]
      ci.lb[i] <-confint(model[[i]])[2,1]
      ci.ub[i] <- confint(model[[i]])[2,2]
    }
  }else if (class(model[[1]]) =="glmerMod"){
    for (i in 1:length(model)){
      beta[i] <- model[[i]]@param[1]
      #ci.lb[i] <-
      #ci.ub[i] <- 
    }
  }else{
    stop('no valid model detected.')
  }
  
  influence_out <- cbind.data.frame('trial'= paste("Omitting" ,unique(df$publication), sep = " ") %>% as.factor(),
                                    "estimate" = transf.ilogit(beta),
                                    "ci.lb" = transf.ilogit(ci.lb),
                                    "ci.ub"= transf.ilogit(ci.ub)) 
  return(influence_out)
}


#print influence plot of leave one out cross validation
PlotInfluence <- function(df, original){
  influence.labs <- paste0(round(df$estimate, digits = 3), " " ,"[",round(df$ci.lb, digits= 3), "-" ,round(df$ci.lb,digits= 3), "]")
  
  orig.estimate <- original[1]
  orig.ci.lb <- original[2]
  orig.ci.ub <- original[3]
  
  plt <- ggplot(df,aes(x = trial , y = estimate) ) +
    geom_point() + 
    scale_y_continuous(limits=c(0,0.75),expand = c(0, 0)) +
    scale_x_discrete()+
    theme_classic() + 
    coord_flip()+
    geom_errorbar(aes(x = trial ,ymin=ci.lb, ymax=ci.ub))+
    annotate( "rect", xmin=0, xmax=Inf, ymin=orig.ci.lb, ymax=orig.ci.ub ,alpha = .2, fill = 'blue') +
    geom_hline(yintercept = orig.estimate,linetype="dashed", 
               color = "blue", size=0.5)+
    annotate("text", label = influence.labs, x = influence_out$trial , y = 0.7, size = 2.5, colour = "black", hjust = 1)+
    theme(
      axis.line.y =element_blank(),
      axis.title.y =element_blank(),
      axis.ticks.y=element_blank()
    )
  
  print(plt)
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

twostep_binorm <- CalcTwostepBiNorm(df, publist)
twostep_binorm.step1 <- twostep_binorm[[1]]
twostep_binorm.step2 <- twostep_binorm[[2]]

twostep_binorm.sum <- summary(twostep_binorm.step2)
twostep_binorm.sum


###################################################################################################

# 2. One-step binomial GLMM allowing for clustering by study. 
# random slope of x within group with correlated intercept
# Laplace approximate ML estimation
# assumes conditional independence and follow binomial distribution

df_onestage <- onehotEncode(df, covar = "publication", names = publist)

onestep_bi_strat <- CalcOnestepBiStrat(df)
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

onestep_bi_rand <- CalcOnestepBiRand(df) 
onestep_bi_rand.sum <- summary(onestep_bi_rand)
onestep_bi_rand.sum

onestep_bi_rand.tau2 <- VarCorr(onestep_bi_rand)[[1]][1] 

###################################################################################################

# 4. Beta-binomial model
# the event rate pi for the ith study is drawn from a beta distirbution. conditional on pi, the number of 
# individuals xi in the ith study of size ni who experience the events of interest follows a binomial 
# distribution B(ni;pi) (From Chuang-Stein 1993).
# Laplace approximate ML estimation

df_props <- CalcProps(df)
twostep_betabi <- CalcTwostepBetaBi(df_props)

twostep_betabi.sum <- summary(twostep_betabi)
twostep_betabi.sum

 
###################################################################################################

# Model comparison: Estimated sumary effects (prop, CI), between study variance (tau, I^)
summary_results <- CalcEstimates(twostep_binorm.step2,
                                 onestep_bi_strat,
                                 onestep_bi_ind,
                                 twostep_betabi)


##Errors with bi_rand and BB confidence intervals (values are completely wrong!)

tau2 <- c(twostep_binorm$tau2, onestep_bi_strat.tau2, onestep_bi_rand.tau2 )

tau2.ci <- 
  
I2 <- 
  
I2.ci <- 
  

modelcomp_df <- cbind.data.frame(c('2step.binorm', '1step_bn_strat', '1step_bn_rand', 'betabinom'),
                                 summary_props, sapply(summary_props.ci95, rbind.data.frame) %>% t())

colnames(modelcomp_df) <- c('model', 'proportions', 'props.ci95_lb', 'props.ci95_ub', 'tau2' , 'tau2.ci95_lb', 
                            'tau2.ci95_ub', I2, I2.ci95_lb, I2.ci95_ub, Q)


###################################################################################################
# Sensitivity Analyses
# Influence of Individual Studies
# 1. Exclusion of small sample sizes (less than n = 10)
# 2. Exclusion of studies with 0 multiple founder variants
# 3. Flipped exclusion criteria for repeat measurements
#  Resampling of participants for which we have multiple measurments 


# Influence of Individual Studies (LOOCV)
df_loocv <- LOOCV.dat(df)[[1]]
publist_loocv <- LOOCV.dat(df)[[2]]

twostep_binorm.influence <- mapply(function(x,y) CalcTwostepBiNorm(data = x, study_list = y)[[2]] , x = df_loocv , y = publist_loocv, SIMPLIFY = FALSE) %>%
  DFInfluence()

pdf("testinfluence.pdf" , width = 14 , height = 20) 
PlotInfluence(twostep_binorm.influence)
dev.off()

onestep_bi_strat.influence <- lapply(df_loo ,CalcOnestepBiStrat) %>%
  DFInfluence()

pdf("testinfluence.pdf" , width = 14 , height = 20) 
PlotInfluence(onestep_bi_strat.influence)
dev.off()

onestep_bi_rand.influence <-  lapply(df_loo ,CalcOnestepBiRand) %>%
  DFInfluence()##TBC due to CIs

pdf("testinfluence.pdf" , width = 14 , height = 20) 
PlotInfluence(onestep_bi_rand.influence)
dev.off()

twostep_betabi.influence <-  lapply(df_loo ,CalcOnestepBiRand) %>%
  DFInfluence() ##TBC due to CIs

pdf("testinfluence.pdf" , width = 14 , height = 20) 
PlotInfluence(twostep_betabi.influence)
dev.off()
  

# 1. Exclusion of small sample sizes (less than n = 10)
publist.nosmallsample <- subset(df_props , subjects > 9 , select = publication) %>% pull(.,var=publication) %>% unique()
df.nosmallsample <- lapply(publist.nosmallsample, function(x,y) subset(x, publication == y), x = df) %>% do.call(rbind.data.frame,.)
df_props.nosmallsample <- subset(df_props , subjects > 9)

twostep_binorm.nosamllsample <- CalcTwostepBiNorm(df.nosmallsample, publist.nosmallsample)
onestep_bi_strat.nosamllsample <- CalcOnestepBiStrat(df.nosmallsample)
onestep_bi_rand.nosamllsample <- CalcOnestepBiRand(df.nosmallsample)
twostep_betabi.nosamllsample <- CalcTwostepBetaBi(df_props.nosmallsample)

SA1_results <- CalcEstimates(twostep_binorm.nosamllsample[[2]],
                             onestep_bi_strat.nosamllsample,
                             onestep_bi_rand.nosamllsample,
                             twostep_betabi.nosamllsample)


# 2. Exclusion of studies with 0 multiple founder variants
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication) %>% pull(.,var=publication) %>% unique()
df.nozeros <- lapply(publist.nozeros, function(x,y) subset(x, publication == y), x = df) %>% do.call(rbind.data.frame,.)
df_props.nozeros <- subset(df_props , multiplefounders != 0)

twostep_binorm.nozeros <- CalcTwostepBiNorm(df.nozeros, publist.nozeros)
onestep_bi_strat.nozeros <- CalcOnestepBiStrat(df.nozeros)
onestep_bi_rand.nozeros <- CalcOnestepBiRand(df.nozeros)
twostep_betabi.nozeros <- CalcTwostepBetaBi(df_props.nozeros) #expectation is this is no better than binomial models

SA2_results <- CalcEstimates(twostep_binorm.nozeros[[2]],
                             onestep_bi_strat.nozeros,
                             onestep_bi_rand.nozeros,
                             twostep_betabi.nozeros)


# 3. Flipped exclusion criteria for repeat measurements


#  Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)


###################################################################################################
#Visualisation


# plot log odds of individual studies. should be normally distiributed to satisfy binomial-normal model.
ggplot(twostep_binorm.step1, aes(x=log_or)) + geom_histogram(binwidth = 0.25,color="black", fill="white")+
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

# Forest Plot 2-step BN
pdf("testplot.pdf" , width = 14 , height = 20)
forest(twostep_binorm.step2 , showweights = TRUE, slab = sort(publist) , transf = transf.ilogit, header = TRUE, digits = 3,
       refline = 0.283, ilab = cbind(df_props$multiplefounders, df_props$subjects), ilab.xpos = c(-0.9, -0.4), 
       )
text(c(-0.9,-0.4), 79, c('Multiple Founders' , 'Subjects') , font = 2)
dev.off()


# Influence dot plot against summary proportion of founder variant multiplicity


# Table summarising sensitivity analyses 1 & 2 compared to original estimations of effect size
Models <- c('Two-Step Binomial Normal',
            'One-Step Binomial (random slope) and correlated intercept',
            'One-Step Binomial (uncorrelated random intercept and slope)',
            "Two-Step Beta-Binomial")

sensitivity_df <- cbind.data.frame(summary_results,
                                   SA1_results,
                                   SA2_results, row.names = Models)

tbl <- kbl(sensitivity_df, digits = 3) %>%
  kable_classic(html_font = "Arial") %>%
  add_header_above(c(" " = 1, 
                     'Summary Estimate' = 3,
                     'Exclusion of Small (<10) Studies' = 3,
                     'Exclusion of Studies reporting 0 MF' = 3))
  
  
tbl
#modify to include number of papers in each analyses
#aim is to include bracketed cis in table that include heterogeneity estimates (tau2 and I2) in addition to raw estimate