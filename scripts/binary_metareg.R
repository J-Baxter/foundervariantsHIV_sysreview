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
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA6. Optimisation Algorithm selected by glmerCrtl
# SA7. Inclusion of unknown sampling delay with repeated studies

###################################################################################################
###################################################################################################
# RUN FROM HERE #
# Dependencies
source('./scripts/load_packages.R')
source('./scripts/generalpurpose_funcs.R')
library(emmeans)

# Wrapper to performance::check_collinearity
# filters according to a tolerance VIF value (default = 5)
CheckCollinearity <- function(model, tol = 5){
  require(performance)
  filtered <- list()
  
  if (class(model) == 'list'){
    check <- lapply(model, performance::check_collinearity %>%
                      as.data.frame())
    
    filtered <- check[which(sapply(check, function(x) max(x[,2])<= tol))]
    
    out <- filtered %>% 
      .[!sapply(.,is.null)] %>%
      do.call(rbind.data.frame, .) %>%
      cbind.data.frame(model = (gsub('\\.[:1,2,3,4:]', '', row.names(.) %>% as.character())))
  }
  
  return(out)
}


# Plot binned residuals
# Y = average residual, X = Founder Variant Multiplicity, Ribbon = SE
PlotBinned <- function(data){
  
  plt<- ggplot(data) + 
    geom_ribbon(aes(x = xbar, ymin = -se, ymax = se), fill = "white", colour = "grey60") + 
    geom_point(aes(x = xbar, y = ybar , colour = group), shape = 4, size = 3)+
    geom_abline(intercept = 0, slope = 0, linetype = "dashed")+
    theme_classic() +
    theme(panel.background = element_rect(fill = 'gray95' )) +
    scale_color_npg() +
    scale_x_continuous(name = 'Probability of Multiple Founders', 
                       labels = scales::percent,
                       limits = c(0,0.73),
                       expand = c(0, 0.005)) +
    scale_y_continuous(name = 'Average Residual', 
                       limits = c(-0.5,0.5),
                       expand = c(0, 0.005))+
    theme(legend.position = "none")
  
  
  return(plt)
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

# Create list of dataframes for leave-one-out cross validation
LOOCV <- function(data, col){
  vars <- unique(data[col]) %>% unlist()
  loo <- list()
  loo.vars <- list()
  
  for (var in vars){
    loo[[var]] <- data[data[col] != var, ]
    loo.vars[[var]] <- vars[vars != var]
  }
  
  out <- list(loo, loo.vars)
  stopifnot(length(loo) == length(loo.vars))
  
  return(out)
}


# Extract estimates from LOOCV to create dataframe (input for influence plot)
DFInfluenceMV <- function(model,labs){
  
  names <- paste("Omitting" , labs %>% names(), sep = " ") %>% 
    as.factor() %>%
    rep(each = 9)
  
  beta = list()
  if (class(model[[1]]) =="glmerMod"){
    for (i in 1:length(model)){
      beta[[i]] <- summary(model[[i]])$coefficients[1:9,] %>% cbind.data.frame()
    }
  }else{
    stop('no valid model detected.')
  }
  
  influence_out <- cbind.data.frame('trial'= names,
                                    "estimate" = do.call(rbind.data.frame,beta))
  
  return(influence_out)
}

# Generate resampled datasets and calculate model estimates for psuedo-bootstrap 
# sensitivity analysis of inclusion/exclusion criteria
BootMetaRegMV <- function(data, replicates){
  require(parallel)
  require(lme4)
  require(dplyr)
  
  resampled <- lapply(1:replicates, function(x,y) {y %>% group_by(participant.id_) %>% sample_n(.,1)},
                      y = data)
  
  resampled_props <- lapply(resampled , CalcProps)
  
  cl <- detectCores() %>%
    `-` (2)
  
  start <- Sys.time()
  print(start)
  
  boot_reg <- mclapply(resampled, CalcRandMetaReg, 
                       formula = model_selected.form,
                       opt = 'bobyqa',
                       mc.cores = cl,
                       mc.set.seed = FALSE)
  
  end <- Sys.time()
  elapsed <- end-start
  print(elapsed)
  
  remove(cl)
  
  boot_reg.coef <- lapply(boot_reg, function(mod) summary(mod)$coefficients[1:9,] %>% cbind.data.frame()) %>%
    do.call(rbind.data.frame,.)
  boot_reg.pred <- lapply(boot_reg, GetPreds, byvar = 'reported.exposure_')%>%
    do.call(rbind.data.frame,.)
  
  
  out <- list(coefs = cbind.data.frame(boot_reg.coef, rep= rep(1:replicates, each=9)),
              preds = boot_reg.pred )
  
  
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


ConcatSA <- function(data){
  int <- list()
  fe <- list()
  re <- list()
  
  for (i in 1:length(data)){
    int[[i]] <- data[[i]]$int
    fe[[i]] <- data[[i]]$fe
    re[[i]] <- data[[i]]$re
  }
  
  int.df <- do.call(rbind.data.frame, int)
  fe.df <- do.call(rbind.data.frame, fe)
  re.df <- do.call(rbind.data.frame, re)
  
  out <- list("int" = int.df,
              "fe" = fe.df,
              "re" = re.df)
  return(out)
}


###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
# Note that this filters the covariates specified and removes levels where n<5
# Also removes unknown exposures
if (!dir.exists('data')){
  Retrieve('data.zip')
}else{
  Sys.sleep(0.2)
}

df <- read.csv("./data/meta_analysis_data.csv",
               na.strings = "NA",
               stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()


# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "env" , "<21", 'NFLG')

df <- SetBaseline(df, baseline.covar, baseline.level)
df$alignment.length_ <- scale(df$alignment.length_)

df_props <- CalcProps(df)


###################################################################################################
###################################################################################################
# STAGE 1: Selecting Random Effects

raneff_forms <- c(r1 = "multiple.founders_ ~  1 + (1 | publication_)",
                  r2 = "multiple.founders_ ~  1 + (1 | publication_) + (1 | cohort_)",
                  r3 = "multiple.founders_ ~  1 + (1 | publication_) + (1 | cohort_) + (1 | cohort_:publication_)")

raneff_effectstruct = GetName(raneff_forms, effects = 'random')

# Run models
raneff_models <- RunParallel(CalcRandMetaReg, raneff_forms, df)

# Check model convergence and singularity
raneff_check <- CheckModels(raneff_models) %>% 
  `row.names<-`(raneff_effectstruct)

# Model selection
raneff_selection <- ModelComp(raneff_models) %>% 
  `row.names<-`(raneff_effectstruct)

# RE Selected = "(1 | publication) + (1|cohort)", significantly p(<0.05) better fit than publication only.
# AIC in agreement, BIC between first two models is indistinguishable


###################################################################################################
# STAGE 2: Selecting Fixed effects to be included in model (bottom up approach)
# Random effects as previously specified
# Baseline covariates: HSX:MTF, phylogenetic, unknown seropositivity, B, env

fixeff_forms<- c(f01 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + (1 | publication_)",
                 f02 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sampling.delay_ + (1 | publication_)",
                 f03 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + (1 | publication_)",
                 f04 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + alignment.bin_ + (1 | publication_)",
                 f05 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + alignment.bin_ + (1 | publication_)",
                 f06 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + sampling.delay_ + (1 | publication_)",
                 f07 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_ + sampling.delay_ + alignment.bin_ + (1 | publication_)")

fixeff_effectstruct <- GetName(fixeff_forms, effects = 'fixed')
fixeff_models <- RunParallel(CalcRandMetaReg, fixeff_forms, df, opt = 'bobyqa') #opt for mac architecture?

# Model diagnostics prior to selection of fixed effects structure
# 1. Identify models that satisfy convergence threshold
# 2. Check Singularity (all values in variance-covariance matrix >0)
# 3. Check multicollinearity between fixed effects (Variance Inflation Factor, VIF >5 removed)

fixeff_check <-  CheckModels(fixeff_models) %>% 
  `row.names<-`(fixeff_effectstruct)

fixeff_multico <- CheckCollinearity(fixeff_models)

fixeff_models.viable <- fixeff_models[which(fixeff_check$is.converged & 
                                              !fixeff_check$is.singular  & 
                                              (names(fixeff_models) %in% unique(fixeff_multico$model)))]

fixeff_forms.viable <- fixeff_forms[which(fixeff_check$is.converged &
                                            !fixeff_check$is.singular & 
                                            (names(fixeff_models) %in% unique(fixeff_multico$model)))]


###################################################################################################
# STAGE 4: Evaluating the inclusion on interactions
# NB A*B = A + B + A:B

interact_forms <- c(i1 = "multiple.founders_ ~ reported.exposure_ + grouped.method_*sequencing.gene_  + sampling.delay_ + (1 | publication_)",
                    i2 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sequencing.gene_*alignment.bin_ + (1 | publication_)",
                    i3 = "multiple.founders_ ~ reported.exposure_ + grouped.method_*sequencing.gene_ + sampling.delay_ + alignment.bin_ + (1 | publication_)",
                    i4 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + sampling.delay_ + sequencing.gene_*alignment.bin_ + (1 | publication_)")

interact_models <- RunParallel(CalcRandMetaReg, interact_forms, df, opt = 'bobyqa') 
interact_effectstruct <- GetName(interact_forms, effects = 'fixed')

interact_check <- CheckModels(interact_models)%>% 
  `row.names<-`(interact_effectstruct)

# Check model convergence and singularity
interact_models.converged <- interact_models[(which(interact_check$is.converged & !interact_check$is.singular))]
interact_forms.converged <- interact_forms[(which(interact_check$is.converged & !interact_check$is.singular))]

# Viable Models
models_viable <- c(fixeff_models.viable , interact_models.converged)
forms_viable <- c(fixeff_forms.viable, interact_forms.converged)
effectstruct_viable <- GetName(forms_viable, effects = 'fixed')

###################################################################################################
# Model selection

# Binned residuals (ideally >95% within SE, but >90% is satisfactory)
binned <- lapply(models_viable, binned_residuals)
binnedplots <- lapply(binned , PlotBinned) %>%
  cowplot::plot_grid(plotlist = ., labels = 'AUTO', nrow = 3)

# function identifies nesting of models to calculate LTR
# Also calculates pseudo R2 and ICC
models_viable.comp <- ModelComp(models_viable) %>% 
  `row.names<-`(effectstruct_viable)

# No significant differences between pairwise LTR, negligble chenge in AIC/BIC
# Model selected = Reported Exposure + Grouped Method + Sequencing Gene + Participant Seropositivity
# Model effects sent to file as part of fixeff_modelbuild.nomultico.effects
model_selected <- models_viable$f06
model_selected.form <- forms_viable[[4]]
model_selected.effectstruct <- GetName(model_selected.form, effects = 'fixed')


###################################################################################################
###################################################################################################
# Extract fixed and random effect coefficients and calculate bootstrapped 95% CIs
start <- Sys.time()
print(start)
model_selected.coef <- GetCoefs(model_selected, model_selected.effectstruct)
end <- Sys.time()
paste('Bootstrapped coefficient estimation was', end-start , 'hours in duration.') %>% print()

# Extract estimated marginal means of fixed effects and calculate 95% CIs
# estimated marginal means average the coefficients of selected vars over all factors
model_selected.marginals <- GetEMM(model = model_selected, 
                                   byvar = 'reported.exposure_', 
                                   label = model_selected.effectstruct)

model_selected.predictions <- GetPreds(model = model_selected, 
                                       byvar = 'reported.exposure_', 
                                       label = model_selected.effectstruct)
###################################################################################################
###################################################################################################
# Validate assumptions of selected model - Includes repeats of multicolinearity and binned residuals
# as described above

# Binned residuals (ideally >95% within SE, but >90% is satisfactory) (repeat of above)
model_selected.resid <- binned_residuals(model_selected) %>% PlotBinned()

setEPS()
postscript("./results/figure_S10.eps", width = 12, height = 12)
Sys.sleep(0.2)
model_selected.resid
dev.off()

# Multicollinearity between fixed effects
model_selected.multico <- check_collinearity(model_selected) %>% plot()

# Normality of Random Effects
model_selected.re <-check_normality(model_selected, effects = 'random') %>% plot()

# ROC curve + AUC
model_selected.roc <- performance_roc(model_selected) %>% plot()

jpeg("./results/model_assumptions.jpeg" ,width = 4500, height = 5000, res = 350 ,units = "px", pointsize = 10)
cowplot::plot_grid(model_selected.resid, model_selected.multico,  plotlist = model_selected.re, model_selected.roc,
                   ncol = 2, labels = 'AUTO', align = 'hv', axis = 'b')
dev.off()


###################################################################################################
###################################################################################################
# Sensitivity analyses on selected model
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA6. Optimisation Algorithm selected by glmerCrtl
# SA7. Inclusion of unknown sampling delay with repeated studies

# SA1. Influence of Individual Studies (LOOCV)
df_loocv <- LOOCV(df, 'publication_')[[1]]
publist_loocv <- LOOCV(df, 'publication_')[[2]]


model_selected.influence <- mclapply(df_loocv, CalcRandMetaReg,
                                     formula = model_selected.form,
                                     mc.cores = 4,
                                     mc.set.seed = FALSE) %>%
  DFInfluenceMV(., labs = publist_loocv) 


# SA2. Exclusion of small sample sizes (less than n = 10)
publist.nosmallsample <- subset(df_props , subjects > 9 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nosmallsample <- df[df$publication_ %in% publist.nosmallsample,]

model_selected.nosmallsample <- CalcRandMetaReg(df.nosmallsample, model_selected.form, opt = 'bobyqa')
model_selected.nosmallsample.out <- list(CheckModels(model_selected.nosmallsample), 
                                         GetCoefs(model_selected.nosmallsample, label = 'no_small'),
                                         GetEMM( model = model_selected.nosmallsample, 
                                                 byvar = 'reported.exposure_', 
                                                 label = 'no_small'),
                                         GetPreds(model = model_selected.nosmallsample, 
                                                  byvar = 'reported.exposure_', 
                                                  label = 'no_small')) 


# SA3. Exclusion of studies with 0 multiple founder variants 
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nozeros <- df[df$publication_ %in% publist.nozeros,]

model_selected.nozeros <- CalcRandMetaReg(df.nozeros, model_selected.form, opt = 'bobyqa')
model_selected.nozeros.out <- list(CheckModels(model_selected.nozeros), 
                                   GetCoefs(model_selected.nozeros, label = 'no_zero'),
                                   GetEMM( model = model_selected.nozeros, 
                                           byvar = 'reported.exposure_', 
                                           label = 'no_zero'),
                                   GetPreds(model = model_selected.nozeros, 
                                            byvar = 'reported.exposure_', 
                                            label = 'no_zero')) 


# SA4. Exclusion of all studies that do not use SGA
publist.sgaonly <- subset(df , sample.amplification_ == 'SGA', select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.sgaonly <- df[df$publication_ %in% publist.sgaonly,]

model_selected.sgaonly <- CalcRandMetaReg(df.sgaonly, model_selected.form, opt = 'bobyqa')
model_selected.sgaonly.out <- list(CheckModels(model_selected.sgaonly), 
                                   GetCoefs(model_selected.sgaonly, label = 'sga_only'),
                                   GetEMM(model = model_selected.sgaonly, 
                                          byvar = 'reported.exposure_', 
                                          label = 'sga_only'),
                                   GetPreds(model = model_selected.sgaonly, 
                                            byvar = 'reported.exposure_', 
                                            label = 'sga_only'))


# SA5. Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)
resampling_df <- read.csv("./data/meta_analysis_data.csv",
                          na.strings = "NA",
                          stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

model_selected.boot_participant <- BootMetaRegMV(resampling_df, 1000) #To re run

# SA6. Optimisation Algorithm selected by glmerCrtl
#opt.algo <- c('bobyqa', 'Nelder_Mead')
#algo <- lapply(opt.algo, function(x) CalcRandMetaReg(df , model_selected.form, opt = x))
#lapply(algo, check_convergence)


# SA7. Delay/Repeat permutation tests
sa7_dflist <- list()

sa7_dflist$sa7_unknown.sing <- read.csv("./data/meta_analysis_data.csv",
                                        na.strings = "NA",
                                        stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_dflist$sa7_unknown.rep <- read.csv("./data/meta_analysis_data.csv",
                                       na.strings = "NA",
                                       stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_dflist$sa7_nounknown.sing <- read.csv("./data/meta_analysis_data.csv",
                                          na.strings = "NA",
                                          stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  filter(sampling.delay_ != 'unknown') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_dflist$sa7_nounknown.rep <- read.csv("./data/meta_analysis_data.csv",
                                         na.strings = "NA",
                                         stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  filter(sampling.delay_ != 'unknown') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_mods <- lapply(sa7_dflist, CalcRandMetaReg, formula = model_selected.form)

sa7_effects <- RunParallel(GetCoefs, sa7_mods, names(sa7_dflist))
###################################################################################################
###################################################################################################
# Outputs to file
# Effect files contain intercept, fixed and random effects including CI
# Selection file includes AIC, loglikelihood, R2, logloss and LRT (as calculated)
# Sensitivity analyses to follow

ifelse(!dir.exists('./results'), dir.create(file.path('./results')), FALSE)

# Selected Model
model_selected.names <- c('multimetareg_int.csv', 'multimetareg_fe.csv', 'multimetareg_re.csv') %>% paste0('./results/', .)
mapply(write.csv, model_selected.coef , file = model_selected.names, row.names = T)

write.csv(model_selected.marginals, './results/multimetareg_emm.csv', row.names = T)

write.csv(model_selected.predictions, './results/multimetareg_preds.csv', row.names = T)

#Model Comp
t3 <- rbind.data.frame(raneff_selection, models_viable.comp )
write.csv(t3, './results/multimetareg_modelselection.csv')


# Sensitivity Analyses
s1 <- model_selected.influence
write.csv(s1, './results/multimetareg_s1.csv')

sa.coefs <- list(model_selected.nosmallsample.out[[2]],
                 model_selected.nozeros.out[[2]],
                 model_selected.sgaonly.out[[2]]) %>%
  ConcatSA()

sacoef.names <- c('multimetareg_s2-4_int.csv', 'multimetareg_s2-4_fe.csv', 'multimetareg_s2-4_re.csv')  %>% paste0('./results/', .)
mapply(write.csv, sa.coefs, file = sacoef.names, row.names = T)  

sa.preds <- list(model_selected.nosmallsample.out[[3]],
                 model_selected.nozeros.out[[3]],
                 model_selected.sgaonly.out[[3]]) %>%
  do.call(rbind.data.frame,.)

sapred.names <- c('multimetareg_s2_pred.csv', 'multimetareg_s3_pred.csv', 'multimetareg_s4_pred.csv')  %>% paste0('./results/', .)
write.csv(sa.preds, './results/multimetareg_sa_pred.csv')  

s5 <- model_selected.boot_participant 
write.csv(s5, './results/multimetareg_s5.csv')

s7 <- Effects2File(sa7_effects)
s7.names <- c('multimetareg_s7_int.csv', 'multimetareg_s7_fe.csv', 'multimetareg_s7_re.csv')  %>% paste0('./results/', .)
mapply(write.csv, s7, file = s7.names, row.names = T)


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################