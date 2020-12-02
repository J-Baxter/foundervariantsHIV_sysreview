#evaluate distribution of data to assumptions of meta analysis

#dependencies
library(fitdistrplus)
library(tidyr)

#functions

#sanity check dataframe, created author/year labels and remove NA from founder multiplicty column
formatDF <- function(df){
  if ('multiple.founders' %in% colnames(df)){
    no_na <- df[!is.na(df$multiple.founders),]
    labelled <-  unite(no_na, "publication", c(author ,year), sep = '_')
  }else{
    labelled <- NULL
    stop('no founder multiplicity column found')
  }
  formatted <- labelled
  return(formatted)
}

#group counts of founder multiplicty by study and selected covariates (NB does not count proportions)
groupDF <- function(df, covar){
  if (class(df$multiple.founders)=="factor"){
    df$multiple.founders = 1 - (as.numeric(df$multiple.founders)-1)
  }else{
    stop('founder multiplicity is not a factor')
  }
  
  #summarise
  if (is.null(covar)){
    df_grouped <- df %>% 
      group_by(publication) %>%
      summarise(subjects = n(), multiple.founders = sum(multiple.founders))
  }else{
    df_grouped <- df %>% 
      group_by(publication, sym(covar)) %>%
      summarise(subjects = n(), multiple.founders = sum(multiple.founders))}
  
  return(df_grouped)
}

#continuity correction for log transformations where proportions = 0 or 1. correction set at 0.3, the expected proportion of founder variants from previous studies
contCorrection <- function(x){
  if(x==0){
    x=0.3
  }else if(x==1){
      x=0.3
  }else{
        x=x}
  return(x)
}

#calculate proportions
calcProps <- function(df, logtransfom = TRUE){
  if ('multiple.founders' %in% colnames(df)){
    props.init <- df$multiple.founders/df_grouped$subjects
  }else{
    stop('no founder multiplicity column found')
  }
  if (logtransfom == TRUE){
    print('calculating log transformed proportions')
    props.corrected <- lapply(props , contCorrection) %>% unlist()
    props <- log(props.corrected)
  }else if (logtransform == FALSE){
    print('calculating raw proportions')
    props <- props.init
  }else{
    stop('specify whether a log transform should be applied by stating logtransfom = TRUE or logtransfom = FALSE')
  }
  return(props)
}


#import dataset
data_master<- read.csv("data_master.csv", na.strings = "NA")



props <- df_grouped$multiplefounders/df_grouped$subjects
props.corrected.log =lapply(props , contcorrection) %>% unlist() %>% log()
fit.norm.cont = fitdist(props.corrected.log , distr = 'norm')
plot(fit.norm.cont)
gofstat(fit.norm.cont)