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

# Set baseline contasts for GLMM
SetBaseline <- function(data,covar,baseline){
  dataframe <- data
  
  stopifnot(length(covar) == length(baseline))
  
  for (i in 1:length(covar)){
    covar. = covar[i]
    baseline. = paste0("(?<!\\S)", baseline[i], "(?!\\S)")
    if(class(dataframe[,covar.]) == 'factor'){
      int <- grep(baseline. , levels(dataframe[,covar.]), perl = T)
      dataframe[,covar.] <- relevel(dataframe[,covar.],  int)
      
      print(levels(dataframe[,covar.])[1])
      
    }else{
      warning('requires factor as input.')
    }
  }
  
  return(dataframe)
}

# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula){
  options(warn = 1)
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
RunParallel <- function(formulas, func, data){
  options(warn = 1)
  
  # Set up cluster (forking)
  cl <- detectCores() %>% `-` (2) 
  
  # Set time and run models
  start <- Sys.time()
  start
  
  metareg <- mclapply(formulas,
                      func, 
                      data = data,
                      mc.cores = cl,
                      mc.set.seed = FALSE) #child process has the same initial random number generator (RNG) state as the current R session
  
  end <- Sys.time()
  elapsed <- end-start
  print(elapsed)
  
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
  options(warn = 1)
  ci <- confint.merMod(model, 
                       method = 'boot', 
                       .progress="txt", 
                       PBargs=list(style=3), 
                       nsim = 10)
  
  if (length(fixef(model)) > 1){
    
    fe <- fixef(model) 
    sd <- sqrt(diag(vcov(model)))
    ci.fe <- ci[-c(1,2),]
    nom <- names(fe)
    
    
    fix_df <- cbind.data.frame(nom = nom,
                               est = fe,
                               sd = sd,
                               ci.lb = ci.fe[,1],
                               ci.ub = ci.fe[,2],
                               analysis = label) %>% 
      `row.names<-` (NULL) %>%
      separate(nom , c('covariate' , 'level') , '_')
    
  }else{
    fe <- fixef(model) 
    ci.fe <- ci[2,]
    nom <- names(fe)
    sd <-  NA
    
    fix_df <- cbind.data.frame(nom = nom,
                               est = fe,
                               sd = sd,
                               ci.lb = ci.fe[1],
                               ci.ub = ci.fe[2],
                               analysis = label) %>% 
      `row.names<-` (NULL) %>%
      separate(nom , c('covariate' , 'level') , '_')
  }
  
  
  return(fix_df)
}


#Random Effects
re <- ranef(model)
re.mean <- lapply(re, function(x) mean(x$`(Intercept)`)) %>%
  do.call(rbind.data.frame,.)
re.sd <- VarCorr(fixeff_modelbuild.converged[[5]])
ci.re <- ci[c(1,2),]
#random_df <- cbind.data.frame(var = names(re),
#est = re.mean,
#sd = sd,
#ci.lb = ci.re[,1],
#ci.ub = ci.re[,2],
#analysis = label, 
#make.row.names =FALSE)
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
# Set seed
set.seed(4472)

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% 
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene'))
  

# Set reference levels for meta regression
# HSX:MTF, phylogenetic, unknown seropositivity, B, env
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "participant.seropositivity_")
baseline.level <- c("HSX:MTF", "phylogenetic", "B" , "env" , "positive")

#df <- SetBaseline(df, baseline.covar, baseline.level)
#df$alignment.length_ <- scale(df$alignment.length_)

###################################################################################################
###################################################################################################
# STAGE 1: Selecting Random Effects

raneff_modelbuild.forms <- c(r0 = "multiple.founders_ ~  1 + (1 | publication_)",
                             r1 = "multiple.founders_ ~  1  + (1 | publication_) + (1|cohort_)",
                             r2 = "multiple.founders_ ~  1  + (1 | publication_) + (1|cohort_) + (1| cohort_:publication_)")

raneff_modelbuild.models <- RunParallel(raneff_modelbuild.forms, CalcRandMetaReg, df)
raneff_modelbuild.models 
raneff.aic <- lapply(raneff_modelbuild.models, AIC)
raneff.bic <- lapply(raneff_modelbuild.models, BIC)
raneff.confint <- lapply(raneff_modelbuild.models, CalcEstimates) %>% do.call(rbind.data.frame, .)
effectstruct = c("(1 | publication)", 
                 "(1 | publication) + (1 | cohort)",
                 "(1 | publication) + (1 | cohort) + (1 | cohort:publication)")

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

fixeff_uni.forms <- c(f0 = "multiple.founders_ ~  1 + (1 | publication_)",
                      f1 = "multiple.founders_ ~  riskgroup_  + (1 | publication_) + (1| cohort_) - 1",
                      f2 = "multiple.founders_ ~ reported.exposure_ + (1 | publication_) + (1| cohort_) - 1",
                      f3 = "multiple.founders_ ~ grouped.method_ + (1 | publication_) + (1| cohort_) - 1",
                      f5 = "multiple.founders_ ~ grouped.subtype_ + (1 | publication_) + (1| cohort_) - 1",
                      f6 = "multiple.founders_ ~ sequencing.gene_ + (1 | publication_) + (1| cohort_) - 1",
                      f7 = "multiple.founders_ ~ participant.seropositivity_ + (1 | publication_) + (1| cohort_) - 1",
                      f8 = "multiple.founders_ ~ alignment.length_ + (1 | publication_) + (1| cohort_) - 1")

fixeff_uni.models <- RunParallel(fixeff_uni.forms, CalcRandMetaReg, df)
fixeff_uni.fe <-mapply(GetFE, model = fixeff_uni.models, label = as.character(fixeff_uni.forms), SIMPLIFY = F) 
plotnames <- lapply(fixeff_uni.fe, GetName)
names(fixeff_uni.fe) <- plotnames 
fixeff_uni.fe_df <- rbindlist(fixeff_uni.fe, idcol = 'names')


###################################################################################################
# STAGE 3: Selecting Fixed effects to be included in model (bottom up approach)
# Random effects as previously specified
# Baseline covariates: HSX:MTF, phylogenetic, unknown seropositivity, B, env

fixeff_modelbuild.forms<- c(f0 = "multiple.founders_ ~  1  + (1 | publication_) + (1 | cohort_)",
                            f1 = "multiple.founders_ ~ reported.exposure_ + (1 | publication_) + (1 | cohort_)",
                            f2 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + (1 | publication_) + (1 | cohort_)",
                            f3 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + (1 | publication_) + (1 | cohort_)",
                            f4 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + sequencing.gene_ + (1 | publication_) + (1 | cohort_)",
                            f5 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + alignment.length_ + (1 | publication_) + (1 | cohort_)",
                            f6 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + grouped.subtype_ + (1 | publication_)",
                            f7 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + sequencing.gene_  + alignment.length_ + (1 | publication_) + (1 | cohort_)",
                            f8 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + sequencing.gene_  +  grouped.subtype_ + (1 | publication_) + (1 | cohort_)",
                            f9 = "multiple.founders_ ~ reported.exposure_ + grouped.method_ + participant.seropositivity_ + sequencing.gene_  +  alignment.length_ + grouped.subtype_ + (1 | publication_) + (1 | cohort_)")

fixeff_modelbuild.models <- RunParallel(fixeff_modelbuild.forms, CalcRandMetaReg, df) 
 
# Model diagnostics prior to selection of fixed effects structure
# 1. Identify models that satisfy convergence threshold
# 2. Check for multicollinearity between fixed effects
# 3. Binned residuals (ideally >95% within SE, but >90% is satisfactory)

# 1. 
fe_convergence <- lapply(fixeff_modelbuild.models, check_convergence) %>% 
  do.call(rbind,.)

fixeff_modelbuild.converged <- fixeff_modelbuild.models[which(fe_convergence)]

# 2.
fe_multico <- lapply(fixeff_modelbuild.converged, check_collinearity)

# 3.
binned <- lapply(fixeff_modelbuild.converged, binned_residuals)
binnedplots <- PlotBinned(binned)

# Extract fixed and random effects estimates from models
FE__modelbuild.converged <- lapply(fixeff_modelbuild.converged, GetEffects)


varcov_mat <- cov2cor(get_varcov(fixeff_modelbuild.models[[6]]))%>%
  reshape2::melt() %>%
  `colnames<-` (c('X','Y','Correlation'))
heatmap <- ggplot(varcov_mat)+geom_tile(aes(x=X,y=Y,fill=Correlation))+
  scale_fill_viridis_c()+
  theme_classic()+
  scale_x_discrete(guide = guide_axis(angle = 50))

fe <- GetEffects(fixeff_modelbuild.models[[7]])
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

binned_residuals()

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