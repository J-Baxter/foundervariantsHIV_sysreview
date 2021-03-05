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
library(data.table)
source('generalpurpose_funcs.R')

# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula){
  f <- as.formula(formula)
  environment(f) <- environment()
  model <- lme4::glmer(f,
                 data = data,
                 family = binomial(link = "logit"),
                 nAGQ = 1,
                 control = glmerControl(optCtrl = list(maxfun = 1000000),
                                        check.nobs.vs.nlev = 'ignore',
                                        check.nobs.vs.nRE = 'ignore',
                                        optimizer = "bobyqa"))
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
GetEffects <- function(model, label = "original"){
  
  ci <- confint.merMod(model, 
                       method = 'boot', 
                       .progress="txt", 
                       PBargs=list(style=3), 
                       nsim = 100)
  modeldf <- model@frame


  
  #Fixed Effects
  fe <- fixef(model)
  sd <- sqrt(diag(vcov(model)))
  ci.fe <- ci[-c(1,2),]
  
  nom <- names(fe) %>% 
    strsplit("_")
  
  fix_df <- cbind.data.frame(level = nom,
                             est = fe,
                             sd = sd,
                             ci.lb = ci.fe[,1],
                             ci.ub = ci.fe[,2],
                             analysis = label, 
                             make.row.names =FALSE)
 
    ref <- list()
    for (i in 1:ncol(fix_df)){
    ref[[i]] = which(param$level %in% modeldf[,i]) %>% cbind.data.frame
    }
    names(ref) <- colnames(fix_df)
    
    covar <- rbindlist(ref,idcol = "names") %>%
    
    dataframe <- data.frame()
    ############# UNFINISHED
    for (i in 1:length(covar[,2])){
      coord <- covar[i,2]
      data.frame[coord,1] <- covar[1,i]
    }
    
    fix_df$covar <- rbindlist(idcol = "names")


  

  #Random Effects
  #re <- ranef(model)
  #re.mean <- lapply(re, mean) %>% do.call(rbind.data.frame,.)
  #re.sd <- ?
  #ci.re <- ci[c(1,2),]
  #random_df <- cbind.data.frame(var = names(re),
                                #est = re.mean,
                                #sd = sd,
                                #ci.lb = ci.re[,1],
                                #ci.ub = ci.re[,2],
                                #analysis = label, 
                                #make.row.names =FALSE)
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

raneff_modelbuild.forms <- c(r0 = "multiple.founders ~  1 + (1 | publication)",
                             r1 = "multiple.founders ~  1  + (1 | publication) + (1|cohort)",
                             r2 = "multiple.founders ~  1  + (1 | publication) + (1|cohort) + (1| cohort:publication)")

raneff_modelbuild.models <- RunMetaReg(raneff_modelbuild.forms,df)
raneff_modelbuild.models 
raneff.aic <- lapply(raneff_modelbuild.models, AIC)
raneff.bic <- lapply(raneff_modelbuild.models, BIC)
raneff.confint <- lapply(raneff_modelbuild.models, CalcEstimates) %>% do.call(rbind.data.frame, .)
effectstruct = c("(1 | publication)", 
                 "(1 | publication) + (1|cohort)",
                 "(1 | publication) + (1|cohort) + (1| cohort:publication)")

raneff_selection <- rbind.data.frame(raneff.aic, raneff.bic) %>% 
  `colnames<-`(effectstruct) %>% 
  t() %>%
  `colnames<-`(c('AIC', 'BIC')) %>%
  cbind.data.frame(raneff.confint[,-c(1,2)])

# RE Selected = "(1 | publication) + (1|cohort)"

###################################################################################################
###################################################################################################
# STAGE 2: Univariate meta-regression of individual covariates against founder variant multiplicity
# Initial regression models with one fixed effect covariate with random effects for publication and cohort
# Equivalent to a subgroup analysis with random effects for subgroup and cohort

fixeff_uni.forms <- c(f0 = "multiple.founders ~  1 + (1 | publication)",
                      f1 = "multiple.founders ~  riskgroup  + (1 | publication) + (1|cohort) - 1",
                      f2 = "multiple.founders ~ reported.exposure + (1 | publication) + (1|cohort) - 1",
                      f3 = "multiple.founders ~ grouped.method + (1 | publication) + (1|cohort) - 1",
                      f5 = "multiple.founders ~ grouped.subtype + (1 | publication) + (1|cohort) - 1",
                      f6 = "multiple.founders ~ sequencing.gene + (1 | publication) + (1|cohort) - 1",
                      f7 = "multiple.founders ~ participant.seropositivity + (1 | publication) + (1|cohort) - 1",
                      f8 = "multiple.founders ~ alignment.length + (1 | publication) + (1|cohort) - 1")

fixeff_uni.models <- RunMetaReg(fixeff_uni.forms, df)
fixeff_uni.fe <-mapply(GetFE, model = fixeff_uni.models, label = as.character(fixeff_uni.forms), SIMPLIFY = F) 
plotnames <- lapply(fixeff_uni.fe, GetName)
names(fixeff_uni.fe) <- plotnames 
fixeff_uni.fe_df <- rbindlist(fixeff_uni.fe, idcol = 'names')


###################################################################################################
# STAGE 3: Selecting Fixed effects to be included in model (bottom up approach)
# Random effects as previously specified

fixeff_modelbuild.forms<- c(f0 = "multiple.founders ~  1  + (1 | publication) + (1|cohort)",
                            f1 = "multiple.founders ~ reported.exposure + (1 | publication) + (1|cohort)-1",
                            f2 = "multiple.founders ~ reported.exposure + grouped.method + (1 | publication) + (1|cohort)-1",
                            f3 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + (1 | publication) + (1|cohort)-1",
                            f4 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene + (1 | publication) + (1|cohort)-1",
                            f5 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + alignment.length + (1 | publication) + (1|cohort)-1",
                            f6 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + grouped.subtype + (1 | publication)-1",
                            f7 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene  + alignment.length + (1 | publication) + (1|cohort)-1",
                            f8 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene  +  grouped.subtype + (1 | publication) + (1|cohort)-1",
                            f9 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene  +  alignment.length + grouped.subtype + (1 | publication) + (1|cohort)-1")

fixeff_modelbuild.models <- RunMetaReg(fixeff_modelbuild.forms,df)

varcov_mat <- cov2cor(get_varcov(fixeff_modelbuild.models[[6]]))%>%
  reshape2::melt() %>%
  `colnames<-` (c('X','Y','Correlation'))
heatmap <- ggplot(varcov_mat)+geom_tile(aes(x=X,y=Y,fill=Correlation))+
  scale_fill_viridis_c()+
  theme_classic()+
  scale_x_discrete(guide = guide_axis(angle = 50))

fe <- GetFE(fixeff_modelbuild.models[[7]])
fe.var <- get_varcov(fixeff_modelbuild.models[[7]])
fixeff.aic <- lapply(fixeff_modelbuild.models, AIC)
fixeff.bic <- lapply(fixeff_modelbuild.models, BIC)
fixeff.confint <- lapply(fixeff_modelbuild.models, confint, method = 'boot' , nsim = ) %>% do.call(rbind.data.frame, .)

ggplot(t) +
  geom_point(aes(x = var, y = transf.ilogit(est))) +
  geom_linerange(aes(x = var,
                     ymin=transf.ilogit(ci.lb),
                     ymax= transf.ilogit(ci.ub)))+
  
  coord_flip()+
  theme_classic()


effectstruct <- c( "1",
                   "reported.exposure", 
                  "reported.exposure + grouped.method",
                  "reported.exposure + grouped.method + participant.seropositivity",
                  "reported.exposure + grouped.method + participant.seropositivity + sequencing.gene",
                  "reported.exposure + grouped.method + participant.seropositivity + alignment.length",
                  "reported.exposure + grouped.method + participant.seropositivity + grouped.subtype",
                  "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene + alignment.length",
                  "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene + grouped.subtype",
                  "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene + alignment.length + grouped.subtype")

fixeff.fit <- rbind.data.frame(fixeff.aic, fixeff.bic) %>% `colnames<-`(effectstruct) %>%
  cbind.data.frame(.,criteria = c('AIC', 'BIC')) %>% reshape2::melt()

fixeff.plot <- cbind.data.frame(model = effectstruct ,
                                fixeff.confint[,-c(1,2)])


###################################################################################################
# STAGE 4: Evaluating the inclusion on interactions
interaction_modelbuild.forms <- c(
  f0 = "multiple.founders ~ reported.exposure + grouped.method + participant.seropositivity + sequencing.gene  +
  (1 | publication) + (1|cohort) - 1",
  f1 = "multiple.founders ~ reported.exposure*grouped.subtype + grouped.method + participant.seropositivity + 
  sequencing.gene + (1 | publication) + (1|cohort) - 1")
  
interaction_modelbuild..models <- RunMetaReg(interaction_modelbuild.forms ,df)


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


###################################################################################################
###################################################################################################
# Outputs to file
write.csv(raneff_selection, file = 'raneff_selection.csv', row.names = T)
write.csv(fixeff_uni.fe_df, file = 'fixeff_uni.csv', row.names = T)
###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################