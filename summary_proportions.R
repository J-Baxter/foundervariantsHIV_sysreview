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


CalcOnestepBiStrat <- function(data){
  model <- glmer(multiple.founders ~  1 +  ( 1 |publication),
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  return(model)
}


CalcOnestepBiRand <- function(data){
  model <- glmer(multiple.founders ~  1 + (1|publication) + (0 + 1|publication),
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  return(model)
}


CalcTwostepBetaBi <- function(proportions){
  model <- betabin(cbind(multiplefounders, subjects - multiplefounders) ~ 1, ~ 1, data = proportions)
  return(model)
}


CalcEstimates <- function(model1, model2, model3, model4, itername){
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
  
  colnames(results) <- c(paste0(itername,".Estimate"), paste0(itername,".95% CI Lower"),paste0(itername,".95% CI Upper") )
  return(results)
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
                                 twostep_betabi,
                                 itername = 'summary')


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


# Influence of Individual Studies
twostep_binorm.influence <- influence.rma.uni(twostep_binorm)
onestep_bi_strat.influence <- influence.ME::influence(onestep_bi_strat , group = "publication")
onestep_bi_rand.influence <- influence.ME::influence(onestep_bi_rand , group = "publication")
twostep_betabi.influence #TBC

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
                             twostep_betabi.nosamllsample,
                             itername = 'SA1')


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
                             twostep_betabi.nozeros,
                             itername = 'SA2')


# Bootstrap replicates of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)


###################################################################################################
#Visualisation


# plot log odds of individual studies. should be normally distiributed to satisfy binomial-normal model.
ggplot(twostep_binorm.step1, aes(x=log_or)) + geom_histogram(binwidth = 0.25,color="black", fill="white")+
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

# Forest Plot 2-step BN
forest()

# Influence dot plot against summary proportion of founder variant multiplicity


# Table summarising sensitivity analyses 1 & 2 compared to original estimations of effect size
models <- c('Two-Step Binomial Normal',
            'One-Step Binomial (random slope) and correlated intercept',
            'One-Step Binomial (uncorrelated random intercept and slope)',
            "Two-Step Beta-Binomial")

sensitivity_df <- cbind.data.frame(models,
                                   summary_results,
                                   SA1_results,
                                   SA2_results)

gt_tbl <- 
  sensitivity_df %>%
  gt(rowname_col = "models") %>%
  tab_stubhead(label = 'Model') %>%
  tab_spanner(label = 'Summary Estimate \n',
              columns = vars(summary.Estimate, `summary.95% CI Lower`, `summary.95% CI Upper`)) %>%
  tab_spanner(label = 'Exclusion of Small (<10) Studies',
              columns = vars(SA1.Estimate, `SA1.95% CI Lower`, `SA1.95% CI Upper`)) %>%
  tab_spanner(label = 'Exclusion of Studies that report 0 MF',
              columns = vars(SA2.Estimate, `SA2.95% CI Lower`, `SA2.95% CI Upper`)) %>%
  cols_label(summary.Estimate = 'Estimate', `summary.95% CI Lower` = '95% CI Lower', `summary.95% CI Upper`= '95% CI Upper',
             SA1.Estimate = 'Estimate', `SA1.95% CI Lower` = '95% CI Lower', `SA1.95% CI Upper` = '95% CI Upper', 
             SA2.Estimate = 'Estimate', `SA2.95% CI Lower` = '95% CI Lower' , `SA2.95% CI Upper` = '95% CI Upper') %>%
  fmt_number(columns = 1:9, decimals = 3) %>%
  tab_options(table_body.hlines.color = 'black',
              table_body.vlines.color = "black",
              column_labels.border.top.color = "black",
              column_labels.border.bottom.color = "black",
              table.border.bottom.color = "black",
              table.border.left.color = 'black',
              )
  
gt_tbl
