###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for IPD meta-regression under a one-step approach, investigating the impact of the use of 
# source sequences in phylogenetic determination of founder variant multiplicity of a acutel infected
# recipient
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
library(insight)
library(emmeans)
library(ggsci)
source('~/foundervariantsHIV_sysreview/generalpurpose_funcs.R')


###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
# Note that this filters the covariates specified and removes levels where n<5
# Also removes unknown exposures
# Role of transmitter in phylogenetic analysis (method = phylo only)

df_phylo <- read.csv("./data/data_master_11121.csv", na.strings = "NA") %>%
  filter(., grouped.method == 'phylogenetic') %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()


# Set reference levels for meta regression
# HSX:MTF, recipient only, unknown seropositivity, B, whole genome
phylo_baseline.covar <- c("reported.exposure_", "seqs.used_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_')
phylo_baseline.level <- c("HSX:MTF", "source&recipient", "B" , "whole.genome" , "<21", 'NFLG')

df_phylo <- df_phylo %>% SetBaseline(., phylo_baseline.covar, phylo_baseline.level)


###################################################################################################
# Run selected model, replacing method variable for inclusion of transmitter seqs
phylo_forms <- c(p1 = 'multiple.founders_ ~ seqs.used_ + (1 | publication_)',
                 p2 = 'multiple.founders_ ~ reported.exposure_ + seqs.used_ + sequencing.gene_ + sampling.delay_ + (1 | publication_)')

phylo_models <- RunParallel(CalcRandMetaReg, phylo_forms, df_phylo , opt = 'bobyqa') 

phylo_effectstruct <- GetName(phylo_forms, effects = 'fixed')

phylo_check <- CheckModels(phylo_models)%>% 
  `row.names<-`(phylo_effectstruct)

phylo_emm <- mapply(GetEMM, phylo_models, label = phylo_effectstruct, byvar = 'seqs.used_', SIMPLIFY = F) %>% do.call(rbind.data.frame,.)
phylo_coefs <- RunParallel(GetCoefs, phylo_models, phylo_effectstruct) %>% Effects2File()


###################################################################################################
# Plots

###################################################################################################
# Out
phylo.names <- c('./results/phylo_int.csv', './results/phylo_fe.csv', './results/phylo_re.csv') %>% paste0('../results/', .)
mapply(write.csv, phylo_coefs , file = phylo.names , row.names = T)

write.csv(phylo_emm , './results/phylo_emm.csv', row.names = T)


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################