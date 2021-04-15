###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Interpretation and meta analysis of number of founder variants
# One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies

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

#
df_num.split <- split.data.frame(df,df$reported.exposure_)

## Plot number of founder variants, stratified by route of exposure
plt_list <- mapply(PlotNumSeqs, data = df_num.split, plt_name = names(df_num.split), SIMPLIFY = F)

cowplot::plot_grid(plotlist = plt_list, ncol = 3, align = 'hv', axis = 'b')

## Summary Stats #founders
df_nosingles <- df %>% filter(minimum.number.of.founders_ != 1)
test <- df_nosingles %>% 
  group_by(reported.exposure_) %>%
  summarise(subjects = n(), 
            mean = mean(minimum.number.of.founders_), 
            median = median(minimum.number.of.founders_), 
            lb = min(minimum.number.of.founders_),
            ub = max(minimum.number.of.founders_)) %>%
  as.data.frame()

write.csv(test, 'numberfounders_summary.csv')

# POisson Reg
poisson_reg <- glmer(minimum.number.of.founders_ ~ reported.exposure_ -1 + (1|publication_), 
                     data = df_nosingles, 
                     family ="poisson" )

exp(poisson_reg@beta)

test$adjusted.mean <- exp(poisson_reg@beta)