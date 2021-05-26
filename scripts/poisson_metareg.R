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
source('./scripts/generalpurpose_funcs.R')


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
CalcTNBMetaReg <- function(data, formula, opt = NULL){
  
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
  model <- glmmTMB::glmmTMB(f, data = data,family = truncated_nbinom2(link = "log"))
  return(model)
}



CalcFounderProbs <- function(model){
  if (class(model) == "glmmTMB"){
    coefs <- fixef(model) %>% unlist()
    
    coefs_grid <- expand.grid(coefs[1:min(5,length(coefs))], c(1,2,3)) %>% 
      {cbind.data.frame(rep(names(coefs[1:min(5,length(coefs))]), 3),.)}
    colnames(coefs_grid) <- c('level', 'beta', 'founders')
    
  } else {
    warning('supported model not detected')
    
  }
  
  exact_probs <- mapply(function(x,y) dpois(x, lambda = exp(y)), y = coefs_grid$beta, x = coefs_grid$founders, SIMPLIFY = F) %>% 
    do.call(rbind.data.frame,.)
  
  cum_probs <- mapply(function(x,y) ppois(x, lambda = exp(y)), y = coefs_grid$beta, x = coefs_grid$founders, SIMPLIFY = F) %>% 
    do.call(rbind.data.frame,.)
  
  colnames(exact_probs) <- 'probability'
  colnames(cum_probs) <- 'probability'
  cum_probs$founders <- rep(c('<=1', "<=2", '<=3'), each = min(5,length(coefs)))
  
  out <- rbind.data.frame(cbind(coefs_grid[,-3],  cum_probs), cbind(coefs_grid,  exact_probs)) %>% arrange(., level)
  levels_order <- c('1', '<=1', '2', "<=2", '3', '<=3')
  out$founders <- factor(out$founders, levels = levels_order)
  
  return(out)
}


###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
# Note that this filters the covariates specified and removes levels where n<5
# Also removes unknown exposures
df <- read.csv("./data/data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  drop_na(minimum.number.of.founders_) %>%
  droplevels()


# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_")
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21")
df <- SetBaseline(df, baseline.covar, baseline.level)



###################################################################################################
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
# truncated negative binomial 
# Horizontal transmission only

df_nosingles <- df %>% filter(minimum.number.of.founders_ != 1) %>% filter(riskgroup_!= 'MTC') %>%  droplevels()

ggplot()+
  geom_bar(data = df_nosingles,
           aes(x = reported.exposure_, fill = as.factor(minimum.number.of.founders_),
               col =as.factor(minimum.number.of.founders_)),
           position = 'stack') +
  theme_bw()+
  scale_fill_viridis_d(name = 'Number of Founders')+
  scale_color_viridis_d(name = 'Number of Founders')+
  theme

trunc_nb2.forms <- c(nb0 = "(minimum.number.of.founders_-1) ~  1 + (1 | publication_)",
                     nb1 = "(minimum.number.of.founders_-1) ~  reported.exposure_ + (1 | publication_)",
                     nb2 = "(minimum.number.of.founders_-1) ~  reported.exposure_ + grouped.method_ + sampling.delay_ + (1 | publication_)")

trunc_nb2.effectstruct = GetName(trunc_nb2.forms, effects = 'random')

trunc_nb2.models <- RunParallel(CalcTNBMetaReg, trunc_nb2.forms, df_nosingles , opt = 'bobyqa')
lapply(trunc_nb2.models, binned_residuals)
trunc_nb2_reg_effects <- lapply(trunc_nb2.models, GetEffects)

trunc_nb2_reg_stratified <- trunc_nb2_reg_effects[[2]]

trunc_nb2_summary <- trunc_nb2_reg_stratified$fe %>%
  select(c(trunc_nb2.mean = est, 
           trunc_nb2.cilb = ci.lb, 
           trunc_nb2.ciub = ci.ub)) %>%
  exp()

out <- cbind.data.frame(quant_summary, poisson_summary)
write.csv(out, '../results/numberfounders_summary.csv')



df_pois <- df %>% filter(riskgroup_!= 'MTC') %>%  droplevels()
test_pois <- glmmTMB(minimum.number.of.founders_ ~  reported.exposure_ + grouped.method_ + sampling.delay_ + sequencing.gene_+ (1|publication_), 
                     data = df_pois ,
                     family = truncated_poisson(link = "log"))
#truncated_nbinom2(link = "log") truncated_poisson(link = "log")
test_probs <- CalcFounderProbs(test_pois)


ggplot()+
  geom_bar(data = test_probs[which(!test_probs$founders %in% c('<=1', "<=2")),], 
           aes(x = founders,
               y = ifelse(grepl('[<=]', test_probs[which(!test_probs$founders %in% c('<=1', "<=2")),]$founders), 1-probability, probability),
               fill = level,
               color = level), 
           position = position_dodge(width = 0.9),
           stat = 'identity') +
  theme_bw() + 
  scale_fill_npg(name = 'Transmission') + 
  scale_colour_npg(guide = NULL)+
  scale_x_discrete(labels = c('1', '2', '3' , '>3'), name = 'Number of Founders') +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.5), name = 'Probability')

###################################################################################################
# Visual inspection and raw summary statistics
# Stratified by route of exposure

df_nosingles$minimum.number.of.founders_ <- as.factor(df_nosingles$minimum.number.of.founders_)
# Plot
df_num.split <- split.data.frame(df_nosingles,df_nosingles$reported.exposure_)
plt_list <- mapply(PlotNumSeqs, data = df_num.split, plt_name = names(df_num.split), SIMPLIFY = F)

jpeg(filename = './results/nof_plot.jpeg', width = 3000, height = 4000, res = 380 ,units = "px", pointsize = 12)
cowplot::plot_grid(plotlist = plt_list, ncol = 3, align = 'hv', axis = 'b')
dev.off()

###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################