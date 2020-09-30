#HIV-1 Founder Variants Meta-Analysis: Principle Dataset
#J Baxter
#created 27/09/2020

#dependencies
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)
library(tidyr)

formatDF <- function(df, covar){
  if (is.null(covar)){
    df_grouped <- df %>% 
      group_by(publication) %>%
      summarise(subjects = n(), multiplefounders = sum(df$founderclass))
  }else{
    df_grouped <- df %>% 
      group_by(publication,covar) %>%
      summarise(subjects = n(), multiplefounders = sum(df$founderclass))}
  
  return(df_grouped)
}

MetaAll <- function(agregated_df) {
  meta::metaprop(multiplefounders,
           subjects,
           studlab = publication,
           data = agregated_df,
           subset = NULL,
           exclude = NULL,
           method = 'Inverse', #Inverse variance method to pool
           sm = "PLN", #log transform proportions
           incr = gs("incr"),#A numeric which is added to event number and sample size of studies with zero or all events, i.e., studies with an event probability of either 0 or 1
           allincr = gs("allincr"),#A logical indicating if incr is considered for all studies if at least one study has either zero or all events. If FALSE (default), incr is considered only in studies with zero or all events
           addincr = gs("addincr"),#A logical indicating if incr is used for all studies irrespective of number of events
           method.ci = "NAsm",
           level = gs("level"),#The level used to calculate confidence intervals for individual studies
           level.comb = gs("level.comb"),#he level used to calculate confidence intervals for pooled estimates
           comb.fixed = FALSE,
           comb.random = TRUE,
           overall = comb.fixed | comb.random,
           overall.hetstat = comb.fixed | comb.random,
           hakn = gs("hakn"), 
           adhoc.hakn = "se",
           method.tau = 'ML', #Maximum likelihood estimation of Tau
           method.tau.ci = "QP", #Q-Profile method (Viechtbauer, 2007)
           tau.preset = NULL,
           TE.tau = NULL,#Overall treatment effect used to estimate the between-study variance tau-squared
           prediction = gs("prediction"),#A logical indicating whether a prediction interval should be printed
           level.predict = gs("level.predict"),#The level used to calculate prediction interval for a new study.
           null.effect = NA,#A numeric value specifying the effect under the null hypothesis.
           method.bias = gs("method.bias"), #A character string indicating which test is to be used. Either "rank", "linreg", or "mm", can be abbreviated. 
  
           #presentation
           backtransf = TRUE,
           pscale = 1,
           title = gs("title"),
           complab = gs("complab"),
           outclab = "",
           print.byvar = gs("print.byvar"),
           byseparator = gs("byseparator"),
           keepdata = TRUE,
           warn = TRUE, 
           control = NULL)
  }

##START##  
#import dataset
principle <- read.csv('sysreview_indivi')

#define subgroups
subgroups <- list(NULL, 'subtype' , 'exposure' , 'seroconversion' , 'method')

#format dataframes for meta-analysis
df_labelled <- unite(df, "publication", c(df$FAU , df$YOP), sep = '_')

data <- lapply(subgroups , function(x) formatDF(df_labelled,x))
  
names(data) <- c('df_master' , 'df_subtype',  'df_seroconversion' , 'df_method')

#pool effect size
meta_init <- lapply(data , MetaAll)

#presentation of initial 'main' meta-analysis (no subgrouping)
#visualisation ideas: forest plots, posterior probaility distributions, l'abbe, GOSH, funnel


#ascertainment of publication bias


#subgroup analysis of principle dataset
metabins4subanalysis <- meta_init[[2:5]]

#presenting subgroup analyses

