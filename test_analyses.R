# IPD meta analysis of HIV founder variant multiplicity
# calculates summary proportion of infections initiated by multiple founder variants
# models implemented:
# 1. One-step glmm (independent intercepts for studies with random effects, ml fit)
# 2. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 3. Two-step beta-binomial model
# estimations of mean effect size, confidence intervals and heterogeneity are presented 


#dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(mltools)
library(data.table)


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
  df_labelled <- unite(df, "publication", c(author ,year), sep = '_')
  return(df_labelled)
}


# Sums number of patients within each study and number of infections initiated by multiple founders. 
# Can be stratified with additional covariates
CalcProps <- function(.data, ...){
    .data %>% 
      group_by(publication, ...) %>%
      summarise(subjects = n(), multiplefounders = sum(multiple.founders))
  }


## Functions for one-step GLMM

# Creates dummy variables for clustering random effects of studies
onehotEncode <- function(data, covar, names){
  onehot <- cbind.data.frame(data[, covar]) %>%
    data.table::as.data.table() %>%
    mltools::one_hot()
  
  colnames(onehot) <- names[order(names)]
  
  testset_onehot <- cbind.data.frame(data,onehot)
  
  return(testset_onehot)
}


# Specifications for fitting of glmm model
glmmREM <- function(data, formula){
  model <- glmer(formula = formula, 
        data = data,
        family = binomial,
        control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))
  return(model)
}


# Calculates confidence intervals at set threshold
CalcCI <- function(u,se,threshold){
  value <- 1-(threshold/2)
  upper <- u + se*qnorm(value)
  lower <- u - se*qnorm(value)
  ci <- c(lower,upper)
  return(ci)
}


## Functions for two-step binomial-normal model

# Logistic regression to pool within study proportions (logit transformed)
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


# Random effects normal model to estimate between study heterogeneity and pool overall proportion (backtransformed)
BNstepTwo <-  function(){
  

}






##########################START##########################

# Import data
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()
  
# Set test data
testlist <- c('Keele_2008' , "Abrahams_2009", "Haaland_2009")
testset_df <- lapply(testlist, function(x,y) subset(x, publication == y), x = df) %>% do.call(rbind.data.frame,.)
 
# 1. One-step  GLMM with clustering (independent intercepts) for studies and random effects for between study heterogeneity
# ML estimation, but could use ASREML (ML estimation may lead to lower E(x) and narrower CIs)
testset_onestage <- onehotEncode(testset_df, covar = "publication", names = testlist)

pool_f <- as.formula('multiple.founders ~ 1 + Keele_2008 + Haaland_2009 + Abrahams_2009 + (1 | publication)')

onestrat.logit <- glmmREM(formula = pool_f,
                           data = testset_onestage) %>% summary()


# 2. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
#step1 pooling within studies
test_summary <- lapply(testlist, BNstepOne, data = testset_df) %>% do.call(rbind.data.frame,.)
  
#step 2 pooling across studies using random effects (normal model). 
meta.ran.reml.hk <- metagen(TE = log_or,
                            seTE= se, 
                            studlab = study_id,
                            data = test_summary,
                            sm = "OR", 
                            comb.random = TRUE,
                            comb.fixed = FALSE,
                            overall = TRUE,
                            method.tau = "REML", 
                            hakn = TRUE, 
                            backtransf = FALSE)

  
#twostep - binomial/normal model coded using metaprop (logit calc + calls metagen internally for random effects normal model)
testset_props <- CalcProps(df)

metaprop.ran <-   meta::metaprop(multiplefounders,
                                 subjects,
                                 studlab = publication,
                                 data = testset_props,
                                 subset = NULL,
                                 exclude = NULL,
                                 method = 'Inverse', #Inverse variance method to pool
                                 sm = "PLOGIT", #log transform proportions
                                 incr = 0.0005,#A numeric which is added to event number and sample size of studies with zero or all events, i.e., studies with an event probability of either 0 or 1
                                 allincr = gs("allincr"),#A logical indicating if incr is considered for all studies if at least one study has either zero or all events. If FALSE (default), incr is considered only in studies with zero or all events
                                 addincr = gs("addincr"),#A logical indicating if incr is used for all studies irrespective of number of events
                                 method.ci = "NAsm",
                                 level = gs("level"),#The level used to calculate confidence intervals for individual studies
                                 level.comb = gs("level.comb"),#he level used to calculate confidence intervals for pooled estimates
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 hakn = TRUE, 
                                 adhoc.hakn = "se",
                                 method.tau = 'REML', #Maximum likelihood estimation of Tau
                                 method.tau.ci = "QP", #Q-Profile method (Viechtbauer, 2007)
                                 prediction = gs("prediction"),#A logical indicating whether a prediction interval should be printed
                                 level.predict = gs("level.predict"),#The level used to calculate prediction interval for a new study.
                                 null.effect = NA,#A numeric value specifying the effect under the null hypothesis.
                                 method.bias = gs("method.bias"), #A character string indicating which test is to be used. Either "rank", "linreg", or "mm", can be abbreviated. 
                                   
                                 #presentation
                                 backtransf = TRUE, pscale = 1, title = gs("title"), complab = gs("complab"), outclab = "",
                                 print.byvar = gs("print.byvar"), byseparator = gs("byseparator"), keepdata = TRUE, warn = TRUE, control = NULL
                                   )
  
  
  forest(metaprop.ran)
  
  
  
  #compare
  model_intercepts <- c( meta.ran.reml.hk$TE.random,
                   metaprop.ran$TE.random,
                   onestrat.logit$coefficients[1,1])
  
  model_se <- c(meta.ran.reml.hk$seTE.random ,
                metaprop.ran$seTE.random,
                onestrat.logit$coefficients[1,2])
  
  precalc_ci <- list(c(meta.ran.reml.hk$lower.random, meta.ran.reml.hk$upper.random),
                  c(metaprop.ran$lower.random, metaprop.ran$upper.random))
  
  onestep_ci <- ci.calc(model_intercepts[3], model_se[3], 0.05) %>% list()
  
  modelcomp_df <- append(precalc_ci , onestep_ci) %>% do.call(rbind.data.frame, . ) %>%
    cbind.data.frame(model_intercepts)
 
  colnames(modelcomp_df) <- c('lower.ci' , 'upper.ci' , 'intercept')
}