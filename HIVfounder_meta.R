#HIV-1 Founder Variants Meta-Analysis
#J Baxter
#created 27/09/2020

#dependencies
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)

#import datasets
principle <- read.csv('sysreview_indivi')

multimethods <- read.csv('sysreview_indivi')

enumeration<- read.csv('sysreview_indivi')

#calculate log odds ratio of multiple founder variants initiating infection for categorical data

CalcOddsRatio <- function(df, covar){
  require(metafor)
  library(dplyr)
  #group by publication, split by covar
  df_comp <- 
  #scale
  df_OR <- escalc(measure="OR", ... , data = df_comp)
  return(df_OR)
}

principle_OR <- CalcOddsRatio(principle)

#