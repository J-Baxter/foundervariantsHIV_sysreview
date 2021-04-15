###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for IPD meta-regression under a one-step and two-step approaches
# One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
# Model Building Process:
# STAGE 1: Selecting Random Effects
# STAGE 2: Univariate meta-regression of individual covariates against founder variant multiplicity
# STAGE 3: Selecting Fixed effects to be included in model (bottom up approach)
# STAGE 4: Evaluating the inclusion on interactions

# Sensitivity analyses conducted on final model:
# SA1. 

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
                 nAGQ = 1,
                 control = cntrl)
  return(model)
}


# Plot binned residuals
# Y = average residual, X = Founder Variant Multiplicity, Ribbon = SE
PlotBinned <- function(data){
  
  if (class(data) == "list"){
    plt_list <- list()
    
    for (i in 1:length(data)){
      plt_list[[i]] <- ggplot(data = data[[i]]) + 
        geom_ribbon(aes(x = xbar, ymin = -se, ymax = se), fill = "white", colour = "grey60") + 
        geom_point(aes(x = xbar, y = ybar , colour = group), shape = 4, size = 3)+
        geom_abline(intercept = 0, slope = 0, linetype = "dashed")+
        theme_classic() +
        theme(panel.background = element_rect(fill = 'gray95' )) +
        scale_color_npg() +
        scale_x_continuous(name = element_blank(), 
                           labels = scales::percent,
                           limits = c(0,0.73),
                           expand = c(0, 0.005)) +
        scale_y_continuous(name = element_blank())+
        theme(legend.position = "none")
    }
  }else{
    warning('data supplied is not a list')
  }
  
  plts <- cowplot::plot_grid(plotlist = plt_list , labels = "AUTO")
  print(plts)
  return(plts)
}


# Calculates AIC, BIC for models, and pairwise LRT (interpret only if appropriate)
ModelComp <- function(modellist){
  len <- length(modellist)-1
  lrt <- list()
  n <- names(modellist)
  for (i in 1:len){
    if(is_nested_models(modellist[[i]], modellist[[i+1]])[1]){
      print(paste('1Models' , n[i] ,'and' , n[i+1], 'are nested.'))
      lrt[[i]] <- anova(modellist[[i]], modellist[[i+1]])
      
    }else if(is_nested_models(modellist[[i-1]], modellist[[i+1]])[1]){
      print(paste('3Models' , n[i-1] ,'and' , n[i+1], 'are nested.'))
      lrt[[i]] <- anova(modellist[[i-1]], modellist[[i+1]])
      
    }else if(is_nested_models(modellist[[i-2]], modellist[[i+1]])[1]){
      print(paste('4Models' , n[i-2] ,'and' , n[i+1], 'are nested.'))
      lrt[[i]] <- anova(modellist[[i-2]], modellist[[i+1]])
      
    }else{
      print('model not nested, comparing to base model')
      lrt[[i]] <- anova(modellist[[1]], modellist[[i]])
    }
    
  }
  rt.df <- do.call(rbind.data.frame, lrt) %>% 
    subset(!duplicated(AIC)) #%>%
    #`row.names<-` (names(modellist))
  
  log.loss <- lapply(modellist, performance_logloss) %>% do.call(rbind.data.frame, .)
  r2 <- lapply(modellist, r2_nakagawa) %>% do.call(rbind.data.frame, .)
  ICC <- lapply(modellist, icc) %>% do.call(rbind.data.frame, .)
  temp <- cbind.data.frame(log.loss, r2, ICC) %>% 
    `colnames<-` (c('log_loss', 'R2_conditional', 'R2_marginal' , 'ICC_adjusted', 'ICC_conditional'))
  out <- cbind.data.frame(rt.df , temp)
  
  return(out)
}


# Pipeline for check_singularity, check_convergence and logloss functions
CheckModels <- function(modellist){
  require(performance)
  
  if (class(modellist) == 'list'){
    is.sing <- lapply(modellist, check_singularity) %>% do.call(rbind.data.frame, .)
    is.con <- lapply(modellist, check_convergence) %>% do.call(rbind.data.frame, .)
  }else{
    is.sing <- check_singularity(modellist) 
    is.con <- check_convergence(modellist) 
  }
  
  
  out <- cbind.data.frame(is.sing, is.con) %>% 
  `colnames<-` (c('is.singular', 'converged'))
  
  rownames(out) <- names(modellist)

  return(out)
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
  
  if (class(model[[1]]) == "rma" || class(model[[1]]) == "rma.uni" || class(model[[1]]) == "rma.glmm"){
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
  }else if (class(model) =="metaprop"){
    beta <- model$TE.random
    ci.lb <- model$lower.random
    ci.ub <- model$upper.random
    
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
BootMetaReg <- function(data, replicates){
  require(parallel)
  require(lme4)
  require(dplyr)
  
  resampled <- lapply(1:replicates, function(x,y) {y %>% group_by(participant.id_) %>% slice_sample(n=1)},
                      y = data)
  
  resampled_props <- lapply(resampled , CalcProps)
  
  cl <- detectCores() %>%
    `-` (2)
 
  start <- Sys.time()
  print(start)
  
  boot_reg <- mclapply(CalcRandMetaReg, resampled, 
                       formula = model_selected.form,
                       opt = 'bobyqa',
                       mc.cores = 4,
                       mc.set.seed = FALSE)
  
  end <- Sys.time()
  elapsed <- end-start
  print(elapsed)
  
  remove(cl)
  
  rand_boot.est <- lapply(rand_boot, function(mod) mod$beta) %>%
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  rand_boot.het <- lapply(rand_boot, function(mod) CalcHet(mod, analysis = "metareg")) %>%
    do.call(rbind.data.frame,.)

  
  out <- cbind.data.frame(rand_boot.est, rand_boot.het)
  
  
  return(out)
}


Effects2File <- function(effectslist){
  int <- lapply(effectslist, function(x) x$int) %>% do.call(rbind.data.frame, .)
  fe <- lapply(effectslist, function(x) x$fe) %>% do.call(rbind.data.frame, .)
  re <- lapply(effectslist, function(x) x$re) %>% do.call(rbind.data.frame, .)
  
  out <- list(int, fe, re)
  names(out) <- c('int', 'fe', 're')
  return(out)
}


###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
# Note that this filters the covariates specified and removes levels where n<5
# Also removes unknown exposures

setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()
  

# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_")
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21")

df <- SetBaseline(df, baseline.covar, baseline.level)
df$alignment.length_ <- scale(df$alignment.length_)

###################################################################################################
###################################################################################################
# STAGE 1: Selecting Random Effects

raneff.forms <- c(r0 = "multiple.founders_ ~  1 + (1 | publication_)",
                  r1 = "multiple.founders_ ~  1 + (1 | publication_) + (1 | cohort_)",
                  r2 = "multiple.founders_ ~  1 + (1 | publication_) + (1 | cohort_) + (1 | cohort_:publication_)")

raneff.effectstruct = GetName(raneff.forms, effects = 'random')

# Run models
raneff.models <- RunParallel(CalcRandMetaReg, raneff.forms, df)

# Check model convergence and singularity
raneff.check <- CheckModels(raneff.models) %>% 
  `row.names<-`(raneff.effectstruct)

# Extract random effects
raneff.effects <- RunParallel(GetEffects, raneff.models, raneff.effectstruct)

# Model selection
raneff.selection <- ModelComp(raneff.models) %>% 
  `row.names<-`(raneff.effectstruct)

# RE Selected = "(1 | publication) + (1|cohort)", significantly p(<0.05) better fit than publication only.
# AIC in agreement, BIC between first two models is indistinguishable

###################################################################################################
# STAGE 2: Selecting Fixed effects to be included in model (bottom up approach)
# Random effects as previously specified
# Baseline covariates: HSX:MTF, phylogenetic, unknown seropositivity, B, env

fixeff_modelbuild.forms<- c(f00 = "multiple.founders_ ~  1  + (1 | publication_) + (1 | cohort_)",
                            f01 = "multiple.founders_ ~ reported.exposure_ + (1 | publication_) + (1 | cohort_)",
                            f02 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + (1 | publication_) + (1 | cohort_)",
                            f03 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sampling.delay_ + (1 | publication_) + (1 | cohort_)",
                            f04 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + (1 | publication_) + (1 | cohort_)",
                            f05 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + alignment.length_ + (1 | publication_) + (1 | cohort_)",
                            f06 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + grouped.subtype_ + (1 | publication_)+ (1 | cohort_)",
                            f07 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + alignment.length_ + (1 | publication_) + (1 | cohort_)",
                            f08 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + sampling.delay_ + (1 | publication_) + (1 | cohort_)",
                            f09 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + grouped.subtype_ + (1 | publication_) + (1 | cohort_)",
                            f10 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + sampling.delay_ + alignment.length_ + (1 | publication_) + (1 | cohort_)",
                            f11 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + sampling.delay_ + grouped.subtype_ + (1 | publication_) + (1 | cohort_)",
                            f12 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + sampling.delay_ + alignment.length_ + grouped.subtype_ + (1 | publication_) + (1 | cohort_)")

fixeff_modelbuild.effectstruct <- GetName(fixeff_modelbuild.forms, effects = 'fixed')
fixeff_modelbuild.models <- RunParallel(CalcRandMetaReg, fixeff_modelbuild.forms, df , opt = 'bobyqa') 

# Model diagnostics prior to selection of fixed effects structure
# 1. Identify models that satisfy convergence threshold
# 2. Check Singularity
# 3. Check for multicollinearity between fixed effects
# 4. Binned residuals (ideally >95% within SE, but >90% is satisfactory)

# 1. & 2. Check model convergence and singularity
fixeff_modelbuild.check <- CheckModels(fixeff_modelbuild.models) %>% 
  `row.names<-`(fixeff_modelbuild.effectstruct)

fixeff_modelbuild.models.converged <- fixeff_modelbuild.models[which(fixeff_modelbuild.check$converged)]
fixeff_modelbuild.forms.converged <- fixeff_modelbuild.forms[which(fixeff_modelbuild.check$converged)]
fixeff_modelbuild.effectstruct.converged <- GetName(fixeff_modelbuild.forms.converged, effects = 'fixed')

# 3. Check for multicollinearity between fixed effects
fe_multico <- lapply(fixeff_modelbuild.models.converged, check_collinearity)
fixeff_modelbuild.models.nomultico <- fixeff_modelbuild.models.converged[-c(7,9)]
fixeff_modelbuild.forms.nomultico <- fixeff_modelbuild.forms.converged[-c(7,9)]
fixeff_modelbuild.effectstruct.nomultico <- GetName(fixeff_modelbuild.forms.nomultico, effects = 'fixed')

# 4. Binned residuals (ideally >95% within SE, but >90% is satisfactory)
binned <- lapply(fixeff_modelbuild.models.nomultico, binned_residuals)
binnedplots <- PlotBinned(binned)

# Extract fixed and random effects for models that satisfy model checks and assumptions
fixeff_modelbuild.nomultico.effects <- RunParallel(GetEffects, fixeff_modelbuild.models.nomultico, fixeff_modelbuild.effectstruct.nomultico)

# Models f0-f5 & f8 converge. f7 and f10 is discounted due to a high degree of collinearity between
# gene.segment and alignment.length
# Final selection concluded following evaluation of interaction terms

###################################################################################################
# STAGE 4: Evaluating the inclusion on interactions
interaction_modelbuild.forms <- c(i0 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sampling.delay_ + sequencing.gene_ + (1 | publication_) + (1 | cohort_)",
                                  i1 = "multiple.founders_ ~ reported.exposure_*grouped.subtype_ + grouped.method_ + sampling.delay_ + sequencing.gene_ + (1 | publication_) + (1 | cohort_)",
                                  i2 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sampling.delay_ + sequencing.gene_*alignment.length_ + (1 | publication_) + (1 | cohort_)")
  
interaction_modelbuild.models <- RunParallel(CalcRandMetaReg, interaction_modelbuild.forms, df , opt = 'bobyqa') 
interaction_modelbuild.effectstruct <- GetName(interaction_modelbuild.forms, effects = 'fixed')

interaction_modelbuild.check <- CheckModels(interaction_modelbuild.models)%>% 
  `row.names<-`(interaction_modelbuild.effectstruct)

# No interaction models converge succesfully

###################################################################################################
# Model selection
# function identifies nesting of models to calculate LTR
fixeff_modelbuild.selection <- ModelComp(fixeff_modelbuild.models.nomultico) %>% 
  `row.names<-`(fixeff_modelbuild.effectstruct.nomultico)
fixeff_modelbuild.multico.check <- CheckModels(fixeff_modelbuild.models.nomultico)%>% 
  `row.names<-`(fixeff_modelbuild.effectstruct.nomultico)

# No significant differences between pairwise LTR, negligble chenge in AIC/BIC
# Model selected = Reported Exposure + Grouped Method + Sequencing Gene + Participant Seropositivity
# Model effects sent to file as part of fixeff_modelbuild.nomultico.effects
model_selected <- fixeff_modelbuild.models.nomultico[[7]]
model_selected.form <- fixeff_modelbuild.forms.nomultico[[7]]
model_selected.effectstruct <- GetName(model_selected.form, effects = 'fixed')


###################################################################################################
###################################################################################################
# Sensitivity analyses on selected model
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA6. Inclusion of unknown sampling delay with repeated studies

# SA1. Influence of Individual Studies (LOOCV)
df_loocv <- LOOCV.dat(df)[[1]]
publist_loocv <- LOOCV.dat(df)[[2]]

model_selected.influence <- mclapply(CalcRandMetaReg, df_loocv, 
                                     formula = model_selected.form,
                                     opt = 'bobyqa',
                                     mc.cores = 4,
                                     mc.set.seed = FALSE) %>%
  DFInfluence(., labs = publist_loocv) 


# SA2. Exclusion of small sample sizes (less than n = 10)
publist.nosmallsample <- subset(df_props , subjects > 9 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nosmallsample <- df[df$publication_ %in% publist.nosmallsample,]

model_selected.nosmallsample <- CalcRandMetaReg(df.nosmallsample, model_selected.form, opt = 'bobyqa')
model_selected.nosmallsample.out <- list(CheckModels(model_selected.nosmallsample), 
                                         GetEffects(model_selected.nosmallsample, label = 'no_small')) 


# SA3. Exclusion of studies with 0 multiple founder variants 
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication) %>%
  pull(.,var=publication) %>%
  unique()

df.nozeros <- df[df$publication %in% publist.nozeros,]

model_selected.nozeros <- CalcRandMetaReg(df.nozeros, model_selected.form, opt = 'bobyqa')
model_selected.nozeros.out <- list(CheckModels(model_selected.nozeros), 
                                   GetEffects(model_selected.nozeros, label = 'no_small')) 


# SA4. Exclusion of all studies that do not use SGA
publist.sgaonly <- subset(df , sample.amplification == 'SGA', select = publication) %>%
  pull(.,var=publication) %>%
  unique()

df.sgaonly <- df[df$publication %in% publist.sgaonly,]

model_selected.sgaonly <- CalcRandMetaReg(df.sgaonly, model_selected.form, opt = 'bobyqa')
model_selected.sgaonly.out <- list(CheckModels(model_selected.sgaonly), 
                                   GetEffects(model_selected.sgaonly, label = 'no_small')) 


# SA5. Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)
resampling_df <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

model_selected.boot_participant <- BootMetaReg(resampling_df , 1000)

# SA6. Optimisation Algorithm selected by glmerCrtl
opt.algo <- c('bobyqa', 'Nelder_Mead')
algo <- lapply(opt.algo, function(x) CalcRandMetaReg(df , model_selected.form, opt = x))
lapply(algo, check_convergence)

# SA7. Delay/Repeat permutation tests
sa7_dflist <- list()

sa7_dflist$sa7_unknown.sing <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_dflist$sa7_unknown.rep <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_dflist$sa7_nounknown.sing <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  filter(sampling.delay_ != 'unknown') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_dflist$sa7_nounknown.rep <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  filter(sampling.delay_ != 'unknown') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_mods <- lapply(sa7_dflist, CalcRandMetaReg, formula = model_selected.form)

sa7_effects <- RunParallel(GetEffects, sa7_mods, names(sa7_dflist))
###################################################################################################
###################################################################################################
# Outputs to file
# Effect files contain intercept, fixed and random effects including CI
# Selection file includes AIC, loglikelihood, R2, logloss and LRT (as calculated)
# Sensitivity analyses to follow

t1 <- Effects2File(raneff.effects) # Error 
t1.names <- c('raneff_int.csv', 'raneff_fe.csv', 'raneff_re.csv')
mapply(write.csv, t1, file = t1.names, row.names = T)

t2 <- Effects2File(fixeff_uni.effects) #fixef_univariate
t2.names <- c('fixef_univariate_int.csv', 'fixef_univariate_fe.csv', 'fixef_univariate_re.csv')
mapply(write.csv, t2, file = t2.names, row.names = T)

t3 <- Effects2File(fixeff_modelbuild.nomultico.effects)
t3.names <- c('fixef_modelbuild_int.csv', 'fixef_modelbuild_fe.csv', 'fixef_modelbuild_re.csv')
mapply(write.csv, t3, file = t3.names, row.names = T)

t4 <- rbind.data.frame(raneff.selection, fixeff_modelbuild.selection)
write.csv(t4, 'model_selection.csv')

s7 <- Effects2File(sa7_effects)
s7.names <- c('s7_int.csv', 's7_fe.csv', 's7_re.csv')
mapply(write.csv, s7, file = s7.names, row.names = T)
###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################