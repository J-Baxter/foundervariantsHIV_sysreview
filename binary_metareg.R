###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for IPD meta-regression under one-step and two-step approaches
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. stratified intercepts and random effects
#    for between study heterogeneity, approx ML fit
# OR. One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
#    (uncorrelated intercept and slope). approx ML fit
# OR. Two-step beta-binomial GLMM, dispersion param for study labels. Laplace approximate ML estimation

###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(mltools)
library(data.table)
library(metafor)
library(dmetar)
library(aod)
library(ggplot2)
library(influence.ME)
library(kableExtra)
source('generalpurpose_funcs.R')

# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula){
  model <- glmer(formula,
                 data = data,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  return(model)
}

###################################################################################################
###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()

set.seed(4472)

###################################################################################################
f0 <- as.formula("multiple.founders ~  1 + (1|publication) + (0 + 1|publication)")
f1 <- as.formula("multiple.founders ~  reported.exposure + (reported.exposure + 1 ||publication) -1")
f2 <- as.formula("multiple.founders ~  direction::reported.exposure + (reported.exposure + 1 ||publication) -1")
test_reg <- CalcRandMetaReg(df, f1)

props_metareg <- CalcProps(df, reported.exposure)

# Model checklist:
# 1. Convergence
# 2. 

#plot to check binned residuals
binned <- performance::binned_residuals(test_reg)
ggplot(data = binned) + 
  geom_ribbon(aes(x = xbar, ymin = -se, ymax = se), fill = "white", colour = "grey60") + 
  geom_point(aes(x = xbar, y = ybar , colour = group), shape = 4, size = 3)+
  geom_abline(intercept = 0, slope = 0, linetype = "dashed")+
  theme_classic() +
  theme(panel.background = element_rect(fill = 'gray95' )) +
  scale_color_npg() +
  scale_x_continuous(name = "Estimated Probability of Multiple Founder Variants", 
                     labels = scales::percent,
                     limits = c(0,0.73),
                     expand = c(0, 0.005)) +
  scale_y_continuous(name = "Average Residual")+
  theme(legend.position = "none")

#plot distribution of within study estimates (eg forest plot) next to modelled modifier
ran.eff <- ranef(test_reg)[[1]] %>% gather()
ran.eff$key <- gsub("reported.exposure", "", ran.eff$key)
ran.df <- ran.eff[ran.eff$key != "(Intercept)", ]
ggplot(data = ran.df,aes( value, key)) + geom_boxplot() +theme_classic()


#compare crude model to moderator model

split_transission(props_metareg ,names)
plot<- tibble(x = props_metareg$multiplefounders/props_metareg$subjects, y = props_metareg$reported.exposure)
ggplot(data = plot, aes(x = x , y = y)) +geom_violin() + theme_classic() + geom_point(data = cbind.data.frame("exposure" = gsub("reported.exposure", "", test_reg@cnms[[2]]), "estimate" = transf.ilogit(test_reg@beta)),aes(y= exposure, x = estimate))


step1 <- escalc(xi = multiplefounders , ni = subjects , data= props_metareg , add = 0.0005, measure = "PLO")
step2 <- rma.uni(yi, vi, 
                 data = step1,
                 method = "REML",
                 knha = TRUE, 
                 measure = "PLO",
                 mods = ~ reported.exposure,
                 intercept = TRUE
                 )
ggplot(data = df_props , aesx = reported.exposure, y = (multiplefounders/subjects)) +geom_density_ridges()


maybe <- cbind.data.frame(test_reg@beta, )




ggplot(cbind.data.frame("exposure" = test_reg@cnms[[2]], "estimate" = transf.ilogit(test_reg@beta)),aes(x= exposure, y = estimate)) +
  geom_point( shape = 4, size = 4) + 
  scale_y_continuous(name = "Pooled Proportion of Founder Variant Multiplicity",limits=c(0,.5),expand = c(0.01, 0.01)) +
  scale_x_discrete()+
  theme_classic() + 
  coord_flip() +
  geom_linerange(aes(ymin=estimate.lb, ymax= estimate.ub, color = analysis), position = position_dodge(0.5)) +
  scale_color_npg(name = NULL, labels = c(
    no_small = "Studies with (n<10) omitted",
    no_zeros = "Studies with (p=0) omitted",
    original = "Full analysis")) + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10.5,  family = "sans"),
        legend.text = element_text(size = 10.5,  family = "sans"),
        axis.title = element_text(size = 13,  family = "sans")
  )