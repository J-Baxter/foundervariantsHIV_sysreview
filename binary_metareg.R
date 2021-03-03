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
source('generalpurpose_funcs.R')

# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula){
  f <- as.formula(formula)
  environment(f) <- environment()
  model <- lme4::glmer(f,
                 data = data,
                 family = binomial(link = "logit"),
                 nAGQ = 1,
                 control = glmerControl(optCtrl = list(maxfun = 100000),
                                        check.nobs.vs.nlev = 'ignore',
                                        check.nobs.vs.nRE = 'ignore'))
  return(model)
}

# Execute a list of lme4 models in parallel
RunMetaReg <- function(formulas, data){
  # Set up cluster (socket)
  cl <- detectCores() %>%
    `-` (2) %>%
    makeCluster()
  
  clusterEvalQ(cl, c(library(lme4), set.seed(4472)))
  
  # Set time and run models
  start <- Sys.time()
  start
  
  metareg <- parLapply(cl = cl, formulas, CalcRandMetaReg, data = data)
  
  end <- Sys.time()
  elapsed <- end-start
  elapsed
  
  stopCluster(cl)
  remove(cl)
  
  return(metareg)
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


# Extract fixed effects from the models
GetFE <- function(model, label = "original"){
  fe <- fixef(model)
  se <- sqrt(diag(vcov(model)))
  nom <- names(fe) %>% 
    gsub("participant.seropositivity|grouped.method|riskgroup|grouped.subtype|sequencing.gene|reported.exposure" , "" , .)
  
  fix_df <- mapply(CalcCI, u=fe, se=se,threshold = 0.05) %>% 
    t() %>% 
    {cbind.data.frame(var = nom,
                      est = transf.ilogit(fe),
                      se = transf.ilogit(se),
                      transf.ilogit(.),
                      analysis = label)}
  colnames(fix_df)[4:5] <- c('fix.ub','fix.lb')
  fix_df

  return(fix_df)
}


# Note required regex for covariate names in order for this to work efficiently
# Extracts th name of the 1st covariate (as written) from lmer function syntax
GetName <- function(x) {
  require(stringr)
  
  formula <- levels(x$analysis)
  
  name <-gsub(".*[:~:] (.+?) [:+:].*", "\\1", formula) %>%
    gsub("[:.:]" , " " , .) %>%
    str_to_title()
  
  return(name)
}


FPlot <- function(data,plotname){
  
  p2 <- ggplot(data = data) + 
    geom_point(aes(x= var, y = est)) + 
    theme_classic() + 
    geom_linerange( aes(x = var, ymin=fix.lb,ymax=fix.ub))+
    scale_y_continuous(limits = c(0,0.75) , 
                       expand = c(0,0), 
                       name = "Probability of Multiple Founders")+
    scale_x_discrete(name = plotname, 
                     guide = guide_axis(angle = 50))+
    #coord_flip()+
    theme(
      axis.text = element_text(size = 10,  family = "sans"),
      legend.text = element_text(size = 10,  family = "sans"),
      axis.title.y = element_text(size = 11,  family = "sans"),
      axis.title.x = element_text(size = 12,  family = "sans"),
      plot.margin = unit(c(2,4,2,1), "lines")
    )

  return(p2)
}


###################################################################################################
###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF(.,filter = c('reported.exposure',
                                                                                     'grouped.subtype',
                                                                                     'sequencing.gene'))

# Set seed
set.seed(4472)

#Calculate IC for sequence information
# aim to be a more informative metric than gene. consists of log (sequence length x evolutionary rate x number of seqs)
# for NGS data, number of seqs is set to 1000
# TBC
evo.rate <- data.frame('env' = exp(-1.84),
                       'not.env' = exp(-2.12),
                       'wg' = )

SetIC <- function(data, rates){
  env <- c('env', 'vif+env+nef')
  not.env <- c('pol', 'gag')
  w.g <- c('whole.genome')
  
  splits <- split.data.frame(data){
    
  }
  data$sequencing.gene <- 
}


###################################################################################################
###################################################################################################
# STAGE 1: Selecting Random Effects

raneff_modelbuild.forms <- c(f0 = "multiple.founders ~  1 + (1 | publication)",
                             f1 = "multiple.founders ~  1  + (1 | publication) + (1|cohort)",
                             f2 = "multiple.founders ~  1  + (1 | publication) + (1|cohort) + (1| cohort:publication)")

raneff_modelbuild.models <- RunMetaReg(raneff_modelbuild.forms,df)
raneff_modelbuild.models 
raneff.aic <- lapply(raneff_modelbuild.models, AIC)
raneff.bic <- lapply(raneff_modelbuild.models, BIC)
raneff.confint <- lapply(raneff_modelbuild.models, CalcEstimates) %>% do.call(rbind.data.frame, .)
effectstruct = c("(1 | publication)", 
                 "(1 | publication) + (1|cohort)",
                 "(1 | publication) + (1|cohort) + (1| cohort:publication)")

raneff.fit <- rbind.data.frame(raneff.aic, raneff.bic) %>% `colnames<-`(effectstruct) %>%
  cbind.data.frame(.,criteria = c('AIC', 'BIC')) %>% reshape2::melt()

raneff.plot <- cbind.data.frame(model = effectstruct ,
                                raneff.confint[,-c(1,2)])

replot <- ggplot(raneff.plot) + 
  geom_point(aes(x = model, y = estimate))+
  geom_linerange(aes(x = model, ymin=estimate.lb, 
                     ymax= estimate.ub))+
  geom_line(aes(x = variable, y = value/2500, color = criteria, group = criteria), data = raneff.fit) +
  geom_point(aes(x = variable, y = value/2500, color = criteria,group = criteria), data = raneff.fit)+
  scale_y_continuous(name = 'Probability of Multiple Founders', expand = c(0,0.02), limits = c(0,1), sec.axis = sec_axis(~.*2500 , name = 'AIC/BIC'))+
  theme_classic() + 
  scale_color_npg()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))


# RE Selected = "(1 | publication) + (1|cohort)"

###################################################################################################
###################################################################################################
# STAGE 2: Univariate meta-regression of individual covariates against founder variant multiplicity
# Initial regression models with one fixed effect covariate with random effects for publication and cohort
# Equivalent to a subgroup analysis with random effects for subgroup and cohort

subgroup_forms <- c(as.formula("multiple.founders ~  riskgroup  + (1 | publication) + (1|cohort) - 1"),
                    as.formula("multiple.founders ~ reported.exposure + (1 | publication) + (1|cohort) - 1"),
                    as.formula("multiple.founders ~ grouped.method + (1 | publication) + (1|cohort) - 1"),
                    as.formula("multiple.founders ~ participant.seropositivity + (1 | publication) + (1|cohort) - 1"),
                    as.formula("multiple.founders ~ grouped.subtype + (1 | publication) + (1|cohort) - 1"),
                    as.formula("multiple.founders ~ sequencing.gene + (1 | publication) + (1|cohort) - 1")
                    )

subgroup_metareg <- RunMetaReg(subgroup_forms, df)
subgroup_fe <-mapply(GetFE, model = subgroup_metareg, label = as.character(subgroup_forms), SIMPLIFY = F) 


plotnames <- lapply(subgroup_fe, GetName)
  
subgroup_plotlist <- mapply(FPlot,subgroup_fe ,plotnames, SIMPLIFY = F )
subgroup_plot <- plot_grid(plotlist = subgroup_plotlist , labels = "AUTO" , align = 'hv', ncol = 2)

subgroup_plot

###################################################################################################
# STAGE 3: Selecting Fixed effects to be included in model (bottom up approach)
# Random effects as previously specified

fixeff_modelbuild.forms<- c(f0 = as.formula("multiple.founders ~  1  + (1 | publication) + (1|cohort)"),
                      
                      f1 = as.formula("multiple.founders ~ reported.exposure + 
                                      (1 | publication) + (1|cohort)"),
                      
                      f2 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                                      (1 | publication) + (1|cohort)"),
                      
                      f3 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                      participant.seropositivity + 
                                      (1 | publication) + (1|cohort)"),
                      
                      f4 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                      participant.seropositivity + sequencing.gene + 
                                      (1 | publication) + (1|cohort)"),
                      
                      f5 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                      participant.seropositivity + sequencing.length + 
                                      (1 | publication) + (1|cohort)"),
                      
                      f6 = as.formula("multiple.founders ~ reported.exposure + grouped.method +
                      participant.seropositivity + sequencing.gene  + sequencing.length + 
                                      (1 | publication) + (1|cohort) - 1"),
                      f7 = as.formula("multiple.founders ~ reported.exposure + grouped.method +
                      participant.seropositivity + sequencing.gene  + sequencing.length + grouped.subtype + 
                                      (1 | publication) + (1|cohort)"))

fixeff_modelbuild.models <- RunMetaReg(fixeff_modelbuild.forms,df)

fixeff.aic <- lapply(fixeff_modelbuild.models, AIC)
fixeff.bic <- lapply(fixeff_modelbuild.models, BIC)
fixeff.confint <- lapply(fixeff_modelbuild.models, CalcEstimates) %>% do.call(rbind.data.frame, .)
effectstruct <- c( "1",
                   "reported.exposure", 
                  "reported.exposure + grouped.method",
                  "reported.exposure + grouped.method + participant.seropositivity",
                  "reported.exposure + grouped.method + participant.seropositivity + sequencing.gene",
                  "reported.exposure + grouped.method + participant.seropositivity + sequencing.length",
                  "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene + sequencing.length",
                  "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene + sequencing.length + grouped.subtype")

fixeff.fit <- rbind.data.frame(fixeff.aic, fixeff.bic) %>% `colnames<-`(effectstruct) %>%
  cbind.data.frame(.,criteria = c('AIC', 'BIC')) %>% reshape2::melt()

fixeff.plot <- cbind.data.frame(model = effectstruct ,
                                fixeff.confint[,-c(1,2)])

replot <- ggplot(fixeff.plot) + 
  geom_point(aes(x = model, y = estimate))+
  geom_linerange(aes(x = model, ymin=estimate.lb, 
                     ymax= estimate.ub))+
  geom_line(aes(x = variable, y = value/2500, color = criteria, group = criteria), data = fixeff.fit) +
  geom_point(aes(x = variable, y = value/2500, color = criteria,group = criteria), data = fixeff.fit)+
  scale_y_continuous(name = 'Probability of Multiple Founders', expand = c(0,0.02), limits = c(0,1), sec.axis = sec_axis(~.*2500 , name = 'AIC/BIC'))+
  theme_classic() + 
  scale_x_discrete(labels = stringr::str_wrap(effectstruct , width = 26), guide = guide_axis(angle = 50))+
  scale_color_npg()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))


###################################################################################################
# STAGE 4: Evaluating the inclusion on interactions
f6 = as.formula("multiple.founders ~ reported.exposure*grouped.subtype + 
                      grouped.method + participant.seropositivity + sequencing.gene + sequencing.number + 
                                      (1 | publication) + (1|cohort) - 1")#,

# f7 = as.formula("multiple.founders ~ reported.exposure + grouped.method +
# participant.seropositivity + sequencing.ic + grouped.subtype + 
#(1 | publication) + (1|cohort) - 1"),

#f8 = as.formula("multiple.founders ~ reported.exposure*grouped.subtype + 
# grouped.method + participant.seropositivity + sequencing.ic + 
# (1 | publication) + (1|cohort) - 1")


###################################################################################################
###################################################################################################
# Evaluation of selected model
# 1. Check Convergence
# 2. Binned residuals (ideally >95% within SE)

# 1. Check Convergence
convergence <- lapply(test_reg, check_convergence, tolerance = 0.05) %>%
  {cbind.data.frame(melt(lapply(., attributes)), melt(.))} %>%
  .[,c(3,4,1)]

colnames(convergence) <- c("model" , "converged", "gradient")

print(convergence)

# 2. Check binned residuals
binned <- lapply(test_reg , binned_residuals)

binnedplots <- PlotBinned(binned)


###################################################################################################
###################################################################################################
# Sensitivity Analyses
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)
# SA5. Optimisation Algorithm selected by glmerCrtl



