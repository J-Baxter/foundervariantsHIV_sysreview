#HIV-1 Founder Variants Meta-Analysis: Principle Dataset
#J Baxter
#created 27/09/2020

#dependencies
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)
library(tidyr)

#import dataset
principle <- read.csv('sysreview_indivi')

#format dataframe for meta-analysis
df_labelled <- unite(df, "publication", c(df$FAU , df$YOP), sep = '_')

df_agg <- df_labelled %>% 
  group_by(publication) %>%
  summarise(N = n(), prop.multiplefounders = mean(df$founderclass, na.rm = TRUE), sd = sd(df$founderclass)) %>%
  mutate(se = sd / sqrt(N),
         lower.ci = prop.multiplefounders - qt(1 - (0.05 / 2), N - 1) * se,
         upper.ci = prop.multiplefounders + qt(1 - (0.05 / 2), N - 1) * se,)

#calculate log odds
log_odds <- mapply(qlogis , df_agg$prop.multiplefounders)

df_lgodds <- cbind.data.frame(df_agg , log_odds)

#pool effect size
metagen(TE = log_odds , seTE = se , studlab = publication , method.tau = 'REML' , sm = '' , data = df_logodds , hakn = TRUE , )

 wrapperfunction <- function(df, covar){
  #dependencies
  require(metafor)
  library(dplyr)
  library(tidyr)

  #group by publication, split by covar (if covar is present...)
  if (is.null(covar)){
    
    #group by publication ony
    df_comp <- df_labelled %>% group_by(publication) %>%
      summarise(mean_size = mean(df$founderclass, na.rm = TRUE), n = n())
    
    }else{
    
    #group by publication, subgroup covar
      df_comp <- df_labelled %>% group_by(publication,covar) %>%
        summarise(mean_size = mean(df$founderclass, na.rm = TRUE), n = n())
  }
  
  #scale
  df_OR <- escalc(measure="OR", ... , data = df_comp)
  
  #output dataframe
  return(df_OR)
}

principle_OR <- CalcOddsRatio(principle)

#