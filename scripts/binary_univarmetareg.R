###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for univariate IPD meta-regression under a one-step approach
# One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies

# Sensitivity analyses conducted on all models:
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA6. Optimisation Algorithm selected by glmerCrtl
# SA7. Inclusion of unknown sampling delay with repeated studies

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
library(emmeans)
source('./scripts/generalpurpose_funcs.R')


# Extract estimates from LOOCV to create dataframe (input for influence plot)
DFInfluenceUV <- function(model,labs){
  
  nfixed <- fixef(model[[1]]) %>% length()
  #print(nfixed)
  #print(which(lapply(model, function(x) fixef(x) %>% length()) %>% unlist() !=nfixed ))
  names <- paste("Omitting" , labs %>% names(), sep = " ") %>% 
    as.factor() %>%
    rep(.,each = nfixed)
  
  beta = list()
  if (class(model[[1]]) =="glmerMod"){
    for (i in 1:length(model)){
      beta[[i]] <- summary(model[[i]])$coefficients %>% cbind.data.frame()
      
      if(nrow(beta[[i]]) != nfixed){
        #print(beta[[i]])
        fixed <- names(fixef(model[[1]]))
        col <- ncol(beta[[i]])
        missing <- fixed[which(!fixed %in% rownames(beta[[i]]))]
        #print(fixed)
        #print(missing)
        #print(rownames(beta[[i]]))
        newrow <- cbind.data.frame(rep(NA, col) %>% rbind()) %>% `rownames<-` (missing)
        names(newrow) <- names(beta[[i]])
        beta[[i]] <- rbind.data.frame(beta[[i]], newrow)
      }
      }
    }else{
      stop('no valid model detected.')
      }
  influence_out <- cbind.data.frame('trial'= names, "estimate" = do.call(rbind.data.frame,beta))
  return(influence_out)
  
  }


GetInfluence <- function(data, form){
  
  df_loocv <- LOOCV.dat(df)[[1]]
  publist_loocv <- LOOCV.dat(df)[[2]]
  
  out <- mclapply(df_loocv, CalcRandMetaReg,
                  formula = form,
                  mc.cores = 4,
                  mc.set.seed = FALSE) %>%
    DFInfluenceUV(., labs = publist_loocv)
  
  out$var <- GetName(form, effects = 'fixed')
  
  return(out)
}


# Generate resampled datasets and calculate model estimates for psuedo-bootstrap 
# sensitivity analysis of inclusion/exclusion criteria
BootMetaRegUV <- function(data, formulas, replicates){
  require(parallel)
  require(lme4)
  require(dplyr)
  
  resampled <- lapply(1:replicates, function(x,y) {y %>% group_by(participant.id_) %>% slice_sample(n=1)},
                      y = data)
  
  out <- list()
  
  for (i in (1:length(formulas))){
    cl <- detectCores() %>%
      `-` (2)
    
    start <- Sys.time()
    print(start)
    
    boot_reg <- mclapply(resampled, CalcRandMetaReg, 
                         formula = formulas[[i]],
                         opt = 'bobyqa',
                         mc.cores = cl,
                         mc.set.seed = FALSE)
    
    end <- Sys.time()
    elapsed <- end-start
    print(elapsed)
    
    remove(cl)
    
    boot_reg.coef <- lapply(boot_reg, function(mod) summary(mod)$coefficients %>% 
                             cbind.data.frame(., label = paste0(GetName(formulas[[i]], effects = 'fixed'), '.boot'))) %>%
      do.call(rbind.data.frame,.)
    
    
    boot_reg.marg <- lapply(boot_reg, 
                            GetEMM,
                            byvar = formulas[[i]],
                            label = paste0(GetName(formulas[[i]], effects = 'fixed'), '.boot')) %>%
      do.call(rbind.data.frame,.)
    
    
    
    out[[i]] <- list(boot_reg.coef, boot_reg.marg)
  }
  
  
  
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
if (!dir.exists('data')){
  Retrieve('data.zip')
}else{
  Sys.sleep(0.2)
}

df <- read.csv("./data/meta_analysis_data.csv",
               na.strings = "NA",
               colClasses=c("multiple.founders"="factor")) %>% 
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21", 'NFLG')

df <- SetBaseline(df, baseline.covar, baseline.level)


df_props <- CalcProps(df)
##################################################################################################
###################################################################################################
# STAGE 2: Univariate meta-regression of individual covariates against founder variant multiplicity
# Initial regression models with one fixed effect covariate with random effects for publication and cohort
# Equivalent to a subgroup analysis with random effects for subgroup and cohort
# Pooled heterogeneity calculation

unipooled_forms <- c(f1 = "multiple.founders_ ~ riskgroup_  + (1 | publication_)",
                  f2 = "multiple.founders_ ~ reported.exposure_  + ( 1 | publication_)",
                  f3 = "multiple.founders_ ~ grouped.method_  + ( 1 | publication_)",
                  f4 = "multiple.founders_ ~ sequencing.gene_ + (1 | publication_)",
                  f5 = "multiple.founders_ ~ alignment.bin_   + (1 | publication_)",
                  f6 = "multiple.founders_ ~ sampling.delay_  + ( 1 | publication_)")

unipooled_effectstruct <- GetName(unipooled_forms, effects = 'fixed')

# Run models
unipooled_models <- RunParallel(CalcRandMetaReg, unipooled_forms, df)

# Check model convergence and singularity
unipooled_check <- CheckModels(unipooled_models)%>% 
  `row.names<-`(unipooled_effectstruct)

# Check model convergence and singularity
unipooled_models.converged <- unipooled_models[which(unipooled_check$is.converged & !unipooled_check$is.singular)]
unipooled_forms.converged <- unipooled_forms[(which(unipooled_check$is.converged & !unipooled_check$is.singular))]
unipooled_effectstruct.converged <- unipooled_effectstruct[(which(unipooled_check$is.converged & !unipooled_check$is.singular))]


###################################################################################################
# Extract fixed and random effect coefficients and calculate bootstrapped 95% CIs
# fixed effects coefficients exponentiated to odds ratios
unipooled_models.coef <- RunParallel(GetCoefs, unipooled_models.converged, unipooled_effectstruct.converged)

# Extract marginal effects of fixed effects and calculate bootstrapped 95% CIs
unipooled_models.marginals <- mapply(GetEMM, model = unipooled_models, 
                                     byvar = as.list(unipooled_forms), 
                                     label = unipooled_effectstruct,
                                     SIMPLIFY = F) %>% do.call(rbind.data.frame,.)


###################################################################################################
# Sensitivity analyses on univariate models
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA6. Inclusion of unknown sampling delay with repeated studies
# SA7. Compare down-sampled to full dataset

# SA1. Influence of Individual Studies (LOOCV)
unipooled_models.influence <- lapply(unipooled_forms, GetInfluence, data = df) %>% 
  do.call(rbind.data.frame,.)


# SA2. Exclusion of small sample sizes (less than n = 10)
publist.nosmallsample <- subset(df_props , subjects > 9 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nosmallsample <- df[df$publication_ %in% publist.nosmallsample,]

unipooled_models.nosmallsample <- RunParallel(CalcRandMetaReg, unipooled_forms, df.nosmallsample, opt = 'bobyqa')

unipooled_models.nosmallsample.out <- list(CheckModels(unipooled_models.nosmallsample), 
                                           RunParallel(GetCoefs, 
                                                       unipooled_models.nosmallsample, 
                                                       paste0(unipooled_effectstruct.converged, '.no_small')),
                                           mapply(GetEMM,
                                                  model = unipooled_models.nosmallsample, 
                                                  byvar = as.list(unipooled_forms),
                                                  paste0(unipooled_effectstruct.converged, '.no_small'), SIMPLIFY = F) %>%
                                             do.call(rbind.data.frame,.)) 


# SA3. Exclusion of studies with 0 multiple founder variants 
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nozeros <- df[df$publication_ %in% publist.nozeros,]

unipooled_models.nozeros <- RunParallel(CalcRandMetaReg, unipooled_forms, df.nozeros, opt = 'bobyqa')

unipooled_models.nozeros.out <- list(CheckModels(unipooled_models.nozeros), 
                                           RunParallel(GetCoefs, 
                                                       unipooled_models.nozeros, 
                                                       paste0(unipooled_effectstruct.converged, '.no_zero')),
                                           mapply(GetEMM,
                                                  model = unipooled_models.nozeros, 
                                                  byvar = as.list(unipooled_forms),
                                                  label = paste0(unipooled_effectstruct.converged, '.no_zero'), SIMPLIFY = F) %>%
                                       do.call(rbind.data.frame,.)) 


# SA4. Exclusion of all studies that do not use SGA
publist.sgaonly <- subset(df , sample.amplification_ == 'SGA', select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.sgaonly <- df[df$publication_ %in% publist.sgaonly,]

unipooled_models.sgaonly <- RunParallel(CalcRandMetaReg, unipooled_forms, df.sgaonly, opt = 'bobyqa')

unipooled_models.sgaonly.out <- list(CheckModels(unipooled_models.sgaonly), 
                                           RunParallel(GetCoefs, 
                                                       unipooled_models.sgaonly, 
                                                       paste0(unipooled_effectstruct.converged, '.sga_only')),
                                           mapply(GetEMM,
                                                  model = unipooled_models.sgaonly, 
                                                  byvar = as.list(unipooled_forms),
                                                  label = paste0(unipooled_effectstruct.converged, '.sga_only'), SIMPLIFY = F) %>% 
                                       do.call(rbind.data.frame,.)) 


# SA5. Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)
resampling_df <- read.csv("./data/meta_analysis_data.csv",
                          na.strings = "NA",
                          colClasses=c("multiple.founders"="factor")) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

unipooled_models.boot_participant <- BootMetaRegUV(resampling_df, unipooled_forms, 1000) 


# SA6. Optimisation Algorithm selected by glmerCrtl - not run

# SA7. Compare down-sampled to full dataset
sa7_dflist <- list()

sa7_dflist$sing <- df
sa7_dflist$rep <- read.csv("./data/meta_analysis_data.csv",
                           na.strings = "NA",
                           colClasses=c("multiple.founders"="factor")) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  SetBaseline(baseline.covar, baseline.level) %>%
  droplevels()

sa7_modgrid <- expand.grid(data = sa7_dflist, formula = unipooled_forms, stringsAsFactors = F)
sa7_mods <- mapply(CalcRandMetaReg, data = sa7_modgrid[[1]] , formula = sa7_modgrid[[2]],SIMPLIFY = F)
sa7_coef<- RunParallel(GetCoefs, sa7_mods, names(sa7_dflist))

sa7_emm <- mapply(GetEMM,
                  model = sa7_mods, 
                  byvar = rep(unipooled_forms, each = 2),
                  label = rep(paste0(unipooled_effectstruct.converged, '.reps'),each = 2), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame,.)
sa7_emm$data <- gsub('\\..*', '\\1'  ,rownames(sa7_emm))



###################################################################################################
###################################################################################################
# Outputs to file
# Effect files contain intercept, fixed and random effects including CI
# Selection file includes AIC, loglikelihood, R2, logloss and LRT (as calculated)
# Sensitivity analyses to follow
ifelse(!dir.exists('./results'), dir.create(file.path('./results/')), FALSE)

# Model Coefficients
t1 <- Effects2File(unipooled_models.coef) # Error 
t1.names <- c('unimetareg_int.csv', 'unimetareg_fe.csv', 'unimetareg_re.csv') %>% paste0('./results/', .)
mapply(write.csv, t1, file = t1.names, row.names = T)

# Model EMM
write.csv(unipooled_models.marginals, './results/unimetareg_emm.csv')

# Sensitivity Analyses
# LOOCV (beta coefficients only)
write.csv(unipooled_models.influence, './results/unimetareg_sa1.csv')

#SA2-4 coefficients
sa2_4.coef <- list(ConcatSA(unipooled_models.nosmallsample.out[[2]]),
                   ConcatSA(unipooled_models.nozeros.out[[2]]),
                   ConcatSA(unipooled_models.sgaonly.out[[2]])) %>%
  ConcatSA() 

sa2_4.coef.names <- c('unimetareg_sa2-4_int.csv', 'unimetareg_sa2-4_fe.csv', 'unimetareg_sa2-4_re.csv')  %>% paste0('./results/', .)
mapply(write.csv, sa2_4.coef, file = sa2_4.coef.names, row.names = T)  


#SA2-4 EMM
sa2_4.emm <- rbind.data.frame(unipooled_models.nosmallsample.out[[3]],
                unipooled_models.nozeros.out[[3]],
                unipooled_models.sgaonly.out[[3]]) 

write.csv(sa2_4.emm, './results/unimetareg_sa2-4_emm.csv', row.names = T)

#SA5 - coefficients
sa5.coef <- lapply(unipooled_models.boot_participant, function(x) x[[1]]) %>% do.call(rbind.data.frame,.)
write.csv(sa5.coef, './results/unimetareg_sa5_coef.csv', row.names = T)  

#SA5 - EMM
sa5.emm <- lapply(unipooled_models.boot_participant, function(x) x[[2]]) %>% do.call(rbind.data.frame,.)
write.csv(sa5.emm, './results/unimetareg_sa5_emm.csv', row.names = T)

#SA7 - Coef
sa7_coef.out <- list(int = lapply(sa7_coef, function(x) x$int) %>% do.call(rbind.data.frame,.),
                     fe = lapply(sa7_coef, function(x) x$fe) %>% do.call(rbind.data.frame,.),
                     re = lapply(sa7_coef, function(x) x$re) %>% do.call(rbind.data.frame,.))

sa7.names <- c('unimetareg_sa7_int.csv', 'unimetareg_sa7_fe.csv', 'unimetareg_sa7_re.csv')  %>% paste0('./results/', .)
mapply(write.csv, sa7_coef.out , file = sa7.names, row.names = T)  

#SA7 - EMM
write.csv(sa7_emm , './results/unimetareg_sa7_emm.csv', row.names = T)

###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################