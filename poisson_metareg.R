###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Interpretation and meta analysis of number of founder variants
# Raw Summary Statistics and Poisson adjusted means
# One-step Poisson GLMM allowing for clustering by study. uncorrelated random effects between studies

# Sensitivity analyses:
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
library(meta)
source('generalpurpose_funcs.R')


# Plots for each founder variant (module)
PlotNumSeqs <- function(data, plt_name){
  require(ggplot2)
  plt <- ggplot(data) +
    geom_bar(aes(minimum.number.of.founders_)) +
    theme_classic()+
    scale_y_continuous(name = "Frequency",
                       expand = c(0, 0))+
    scale_x_discrete(name = "Number of Founder Variants",
                     expand = c(0, 0))+
    labs(title = plt_name)
  
  return(plt)
}


# One-step GLMM accounting for clustering of studies using a random intercept
CalcPoissonMetaReg <- function(data, formula, opt = NULL){
  
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
                       family ="poisson",
                       nAGQ = 1,
                       control = cntrl)
  return(model)
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
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  drop_na(minimum.number.of.founders_) %>%
  droplevels()


# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "delay.sampling_")
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21")

df <- SetBaseline(df, baseline.covar, baseline.level)
df$alignment.length_ <- scale(df$alignment.length_)
df$minimum.number.of.founders_ <- as.factor(df$minimum.number.of.founders_)

###################################################################################################
###################################################################################################
# Visual inspection and raw summary statistics
# Stratified by route of exposure

# Plot
df_num.split <- split.data.frame(df,df$reported.exposure_)
plt_list <- mapply(PlotNumSeqs, data = df_num.split, plt_name = names(df_num.split), SIMPLIFY = F)

jpeg(filename = 'sa7_plot.jpeg', width = 3000, height = 4000, res = 380 ,units = "px", pointsize = 12)
cowplot::plot_grid(plotlist = plt_list, ncol = 3, align = 'hv', axis = 'b')
dev.off()


# Summary Stats - first removing all single founder infections
df_nosingles <- df %>% filter(minimum.number.of.founders_ != 1)

quant_summary <- df_nosingles %>% 
  group_by(reported.exposure_) %>%
  summarise(subjects = n(), 
            mean = mean(minimum.number.of.founders_), 
            median = median(minimum.number.of.founders_), 
            lb = min(minimum.number.of.founders_),
            ub = max(minimum.number.of.founders_)) %>%
  as.data.frame()


###################################################################################################
###################################################################################################
# Regression of number of multiple founders against reported exposure
# Random effects of study only
# Poisson log link

poisson.forms <- c(p0 = "minimum.number.of.founders_ ~  1 + (1 | publication_)",
                  p1 = "minimum.number.of.founders_ ~  reported.exposure_ -1 + (1 | publication_)")

poisson.effectstruct = GetName(poisson.forms, effects = 'random')

poisson.models <- RunParallel(CalcPoissonMetaReg, poisson.forms, df_nosingles , opt = 'bobyqa')

poisson_reg_effects <- lapply(poisson.models, GetEffects)

poisson_reg_stratified <- poisson_reg_effects[[2]]

poisson_summary <- poisson_reg_stratified$fe %>%
  select(c(poisson.mean = est, 
           poisson.cilb = ci.lb, 
           poisson.ciub = ci.ub)) %>%
  exp()

out <- cbind.data.frame(quant_summary, poisson_summary)
write.csv(out, 'numberfounders_summary.csv')


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################