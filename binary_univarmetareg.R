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


###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% 
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21", 'NFLG')

df <- SetBaseline(df, baseline.covar, baseline.level)


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
models_converged.coef <- RunParallel(GetCoefs, unipooled_models.converged, unipooled_effectstruct.converged )

# Extract marginal effects of fixed effects and calculate bootstrapped 95% CIs
models_converged.marginals <- mapply(GetEMM, model = unipooled_models, 
                                     byvar = as.list(unipooled_forms), 
                                     label = unipooled_effectstruct,
                                     SIMPLIFY = F)

###################################################################################################
# Sensitivity analyses on selected model
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA6. Inclusion of unknown sampling delay with repeated studies
# SA7. Compare down-sampled to full dataset

# SA1. Influence of Individual Studies (LOOCV)
df_loocv <- LOOCV.dat(df)[[1]]
publist_loocv <- LOOCV.dat(df)[[2]]


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
                                         GetEffects(model_selected.nosmallsample, label = 'no_small')) 


# SA3. Exclusion of studies with 0 multiple founder variants 
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nozeros <- df[df$publication_ %in% publist.nozeros,]

model_selected.nozeros <- CalcRandMetaReg(df.nozeros, model_selected.form, opt = 'bobyqa')
model_selected.nozeros.out <- list(CheckModels(model_selected.nozeros), 
                                   GetEffects(model_selected.nozeros, label = 'no_zero')) 


# SA4. Exclusion of all studies that do not use SGA
publist.sgaonly <- subset(df , sample.amplification_ == 'SGA', select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.sgaonly <- df[df$publication_ %in% publist.sgaonly,]

model_selected.sgaonly <- CalcRandMetaReg(df.sgaonly, model_selected.form, opt = 'bobyqa')
model_selected.sgaonly.out <- list(CheckModels(model_selected.sgaonly), 
                                   GetEffects(model_selected.sgaonly, label = 'sga_only')) 


# SA5. Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)
resampling_df <- read.csv("data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay'), noreps = FALSE) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

model_selected.boot_participant <- BootMetaRegMV(resampling_df , 1000) #To re run

# SA6. Optimisation Algorithm selected by glmerCrtl
opt.algo <- c('bobyqa', 'Nelder_Mead')
algo <- lapply(opt.algo, function(x) CalcRandMetaReg(df , model_selected.form, opt = x))
lapply(algo, check_convergence)


# SA7. Compare down-sampled to full dataset



###################################################################################################
###################################################################################################
# Outputs to file
# Effect files contain intercept, fixed and random effects including CI
# Selection file includes AIC, loglikelihood, R2, logloss and LRT (as calculated)
# Sensitivity analyses to follow

dir.create(file.path('../results')) 

# Model Build
t1 <- Effects2File(raneff.effects) # Error 
t1.names <- c('raneff_int.csv', 'raneff_fe.csv', 'raneff_re.csv') %>% paste0('../results/', .)
mapply(write.csv, t1, file = t1.names, row.names = T)

t2 <- Effects2File(fixeff_modelbuild.nomultico.effects)
t2.names <- c('fixef_modelbuild_int.csv', 'fixef_modelbuild_fe.csv', 'fixef_modelbuild_re.csv') %>% paste0('../results/', .)
mapply(write.csv, t2, file = t2.names, row.names = T)

t3 <- rbind.data.frame(raneff.selection, fixeff_modelbuild.selection)
write.csv(t3, '../results/model_selection.csv')



###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################