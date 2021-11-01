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
source('./scripts/generalpurpose_funcs.R')


###################################################################################################
###################################################################################################
# Set seed
set.seed(4472)

# Import data
# Note that this filters the covariates specified and removes levels where n<5
# Also removes unknown exposures
# Role of transmitter in phylogenetic analysis (method = phylo only)
if (!dir.exists('data')){
  Retrieve('data.zip')
}else{
  Sys.sleep(0.2)
}

df_phylo <- read.csv("./data/meta_analysis_data.csv",
                     na.strings = "NA",
                     stringsAsFactors = T) %>%
  filter(., grouped.method == 'phylogenetic:R' |grouped.method == 'phylogenetic:S&R') %>%
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

lapply(phylo_models, binned_residuals)

phylo_emm <- mapply(GetEMM, phylo_models, label = phylo_effectstruct, byvar = 'seqs.used_', SIMPLIFY = F) %>% do.call(rbind.data.frame,.)
phylo_coefs <- RunParallel(GetCoefs, phylo_models, phylo_effectstruct) %>% Effects2File()


###################################################################################################
# Plots

phylo_plot.in <- rbind.data.frame(phylo_coefs$int,phylo_coefs$fe) %>% filter(., covariate == 'seqs.used' | covariate =='(Intercept)')
phylo_plot.in[c(1,2), c(3,5,6)] <- 0

phylo_plot.in[c(1,2), 2] <- 'source and recipient'

phylo_plt <- ggplot(phylo_plot.in) +
  geom_point(aes(x= exp(est), 
                 y = factor(level),
                 col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C')),
                 shape = analysis),
             size = 5,
             position = position_dodge2(width = 0.6)) +
  theme_bw() + 
  scale_shape(name = 'Model', labels = c('Univariate', 'Multivariate'))+
  geom_linerange(aes(y = factor(level),
                     xmin= exp(ci.lb), 
                     xmax= exp(ci.ub), 
                     col = ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
                 position = position_dodge2(width = 0.6))+
  scale_x_continuous(
    expand = c(0,0), 
    name = "Odds Ratio",
    #trans = 'log10'
  )+
  scale_y_discrete(name = 'Sequences Analysed',labels = c( 'Recipient Only','Source & Recipient'))+
  scale_colour_manual(values = setNames(c("#E64B35FF", "#4DBBD5FF", '#000000'), c('A',"B","C")), guide = NULL) +
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  coord_cartesian(xlim = c(0,5)) 

jpeg(filename = './results/metareg_phylo.jpeg', width = 2000, height = 2000, res = 380 ,units = "px", pointsize = 12)

phylo_plt 

dev.off() 
###################################################################################################
# Out
phylo.names <- c('phylo_int.csv', 'phylo_fe.csv', 'phylo_re.csv') %>% paste0('./results/', .)
mapply(write.csv, phylo_coefs , file = phylo.names , row.names = T)

write.csv(phylo_emm , './results/phylo_emm.csv', row.names = T)


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################