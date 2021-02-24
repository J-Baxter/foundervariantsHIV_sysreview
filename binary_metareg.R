###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# Framework for IPD meta-regression under a one-step and two-step approaches
# One-step binomial GLMM allowing for clustering by study. uncorrelated random effects between studies
# 1. Initial regression models with one fixed effect covariate with random effects for publication and cohort
# 2. Hierachichal Model Building with additional parameters

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

###################################################################################################
# Initial regression models with one fixed effect covariate with random effects for publication and cohort
# Equivalent to a subgroup analysis

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

# Outline formulas for meta-regression (Hierarchical model fitting)
# Models 'A' use broadly defined risk groups (MSM, HSX etc) whereas models 'B', subgroup by direction
# and/or timing

modelbuild_forms <- c(f0 = as.formula("multiple.founders ~  1  + (1 | publication) + (1|cohort) - 1"),
                      
                      f1 = as.formula("multiple.founders ~ reported.exposure + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f2 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f3 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                      participant.seropositivity + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f4 = as.formula("multiple.founders ~ reported.exposure + grouped.method + 
                      participant.seropositivity + sequencing.gene + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f5 = as.formula("multiple.founders ~ reported.exposure + grouped.method +
                      participant.seropositivity + sequencing.gene + sequencing.number + grouped.subtype + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f6 = as.formula("multiple.founders ~ reported.exposure*grouped.subtype + 
                      grouped.method + participant.seropositivity + sequencing.gene + sequencing.number + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f7 = as.formula("multiple.founders ~ reported.exposure + grouped.method +
                      participant.seropositivity + sequencing.ic + grouped.subtype + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
                      f8 = as.formula("multiple.founders ~ reported.exposure*grouped.subtype + 
                      grouped.method + participant.seropositivity + sequencing.ic + 
                                      (1 | publication) + (1|cohort) - 1"),
                      
           )


###################################################################################################


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
binned <- lapply(test_reg , binned_residuals)

binnedplots <- PlotBinned(binned)

# 3. Extract AIC and compare


#plot distribution of within study estimates (eg forest plot) next to modelled modifier
ran.eff <- ranef(test_reg_2)[[1]] %>% gather()
ran.eff$key <- gsub("reported.exposure", "", ran.eff$key)
ran.df <- ran.eff[ran.eff$key != "(Intercept)", ]
ggplot(data = ran.df,aes( value, key)) + geom_boxplot() +theme_classic() 



simple <-c(f1a = as.formula("multiple.founders ~  riskgroup  + (1 | publication) - 1"),
           f1a = as.formula("multiple.founders ~  riskgroup  + (1 | publication) + (1|cohort) - 1"),
           f1a = as.formula("multiple.founders ~  riskgroup  + grouped.method + (1 | publication) + (1|cohort) - 1"),
           f1a = as.formula("multiple.founders ~  riskgroup  + grouped.method + participant.seropositivity + (1 | publication) + (1|cohort) - 1")
          )

simple2 <-c(as.formula("multiple.founders ~  riskgroup  + (1 | publication) + (1|cohort) - 1"),
            as.formula("multiple.founders ~ reported.exposure + (1 | publication) + (1|cohort) - 1"),
            as.formula("multiple.founders ~ grouped.method + (1 | publication) + (1|cohort) - 1"),
            as.formula("multiple.founders ~ participant.seropositivity + (1 | publication) + (1|cohort) - 1"),
            as.formula("multiple.founders ~ grouped.subtype + (1 | publication) + (1|cohort) - 1"),
            as.formula("multiple.founders ~ sequencing.region + (1 | publication) + (1|cohort) - 1")
)

test_reg2 <- lapply(simple2, CalcRandMetaReg, data = df)



fe2 <- mapply(GetFE, model = test_reg2, label = c("~  riskgroup  + (1 | publication) + (1|cohort) - 1",
                                                  "~  reported.exposure + (1 | publication) + (1|cohort) - 1",
                                                '~ grouped.method + (1 | publication) + (1|cohort) - 1',
                                               '~ participant.seropositivity + (1 | publication) + (1|cohort) - 1',
                                               "~ grouped.subtype + (1 | publication) + (1|cohort) - 1",
                                               "~ sequencing.region + (1 | publication) + (1|cohort) - 1"), SIMPLIFY = F) %>% do.call(rbind.data.frame, .)

p1 <- ggplot(data = fe) + 
  geom_point(aes(x= var, y = est , colour = analysis) ,position = position_dodge(0.5)) + 
  theme_classic() + 
  geom_linerange( aes(x = var, ymin=fix.lb,ymax=fix.ub, color = analysis),position = position_dodge(0.5))+
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )+guides(col = guide_legend(nrow=2))+
  coord_flip()



