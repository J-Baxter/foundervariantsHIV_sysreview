###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for IPD meta-regression under a one-step and two-step approaches
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
#    (uncorrelated intercept and slope). approx ML fit

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
source('generalpurpose_funcs.R')

# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula){
  model <- lme4::glmer(formula,
                 data = data,
                 family = binomial(link = "logit"),
                 nAGQ = 1,
                 control = glmerControl(optCtrl = list(maxfun = 100000),
                                        check.nobs.vs.nlev = 'ignore',
                                        check.nobs.vs.nRE = 'ignore'))
  return(model)
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


###################################################################################################
###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()

# Set seed
set.seed(4472)

###################################################################################################

# Outline formulas for meta-regression (Hierarchical model fitting)
# Models 'A' use broadly defined risk groups (MSM, HSX etc) whereas models 'B', subgroup by direction
# and/or timing

forms <- c(f0 = as.formula("multiple.founders ~  1 + 
                           (1 | publication)"),
           
           f1a = as.formula("multiple.founders ~  riskgroup + 
                           (1 | cohort) + (1 | publication) - 1"),
           f1b = as.formula("multiple.founders ~  reported.exposure + 
                           (1 | cohort) + (1 | publication) - 1"),
           
           f2a = as.formula("multiple.founders ~  riskgroup + grouped.method + 
                           (1| cohort) + (1 | publication) - 1"),
           f2b = as.formula("multiple.founders ~  reported.exposure + grouped.method + 
                           (1| cohort) + (1 | publication) - 1"),
           
           f3a = as.formula("multiple.founders ~  riskgroup + grouped.method + sequencing.region + 
                           (1 | cohort) + (1 | publication) - 1"),
           f3b = as.formula("multiple.founders ~ reported.exposure + grouped.method + sequencing.region + 
                           (1 | cohort) + (1 | publication) - 1")#,
           
           #f5 = as.formula("multiple.founders ~  riskgroup*grouped.subtype + grouped.method + sequencing.region  + seropositivity + (1 | cohort) + (1 | publication) - 1"),
           #f6 = as.formula("multiple.founders ~  reported.exposure*grouped.subtype + grouped.method + sequencing.region  + seropositivity + (1 | cohort) + (1 | publication) - 1")
           )


###################################################################################################

# Set up cluster (socket)
cl <- detectCores() %>%
  `-` (2) %>%
  makeCluster()

clusterEvalQ(cl, c(library(lme4), set.seed(4472)))

# Set time and run models
start <- Sys.time()
start

test_reg <- parLapply(cl = cl, forms, CalcRandMetaReg, data = df_splittrans)

end <- Sys.time()
elapsed <- end-start
elapsed

stopCluster(cl)
remove(cl)

###################################################################################################
###################################################################################################

# Model Evaluation
# 1. Check Convergence
# 2. Binned residuals (ideally >95% within SE)
# 3. Compare AIC

# 1. Check Convergence
convergence <- lapply(test_reg, check_convergence, tolerance = 0.05) %>%
  {cbind.data.frame(melt(lapply(., attributes)), melt(.))} %>%
  .[,c(3,4,1)]

colnames(convergence) <- c("model" , "converged", "gradient")

print(convergence)

# 2. Check binned residuals
#plot to check binned residuals
binned <- lapply(test_reg , binned_residuals)

binnedplots <- PlotBinned(binned)

#plot distribution of within study estimates (eg forest plot) next to modelled modifier
ran.eff <- ranef(test_reg_2)[[1]] %>% gather()
ran.eff$key <- gsub("reported.exposure", "", ran.eff$key)
ran.df <- ran.eff[ran.eff$key != "(Intercept)", ]
ggplot(data = ran.df,aes( value, key)) + geom_boxplot() +theme_classic() 
fix.eff <- cbind.data.frame(gsub("reported.exposure", "", test_reg@cnms[[2]]),fixef(test_reg))

ggplot(data = fix.eff) + geom_point(aes(y= exposure, x = estimate)) +theme_classic
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

