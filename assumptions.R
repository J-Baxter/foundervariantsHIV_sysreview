#evaluate distribution of data to assumptions of meta analysis

#dependencies
setwd("~/foundervariantsHIV_sysreview")
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
groupDF <- function(df, covar = NULL){
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
      group_by(publication, get(covar)) %>%
      summarise(subjects = n(), multiple.founders = sum(multiple.founders))
    }
  
  return(df_grouped)
}


#continuity correction for log transformations where proportions = 0 or 1. correction set at 0.3, the expected proportion of founder variants from previous studies
contCorrection <- function(x){
  if(x == 0 | x == 1){
    x <- 0.3
    
  }else{
    x <- x}
  
  return(x)
}


#calculate proportions
calcProps <- function(df, logtransform = TRUE){
  if ('multiple.founders' %in% colnames(df)){
    props.init <- df$multiple.founders/df$subjects
    
  }else{
    stop('no founder multiplicity column found')
  }
  
  if (logtransform == TRUE){
    cat('\ncalculating log transformed proportions \n')
    props.corrected <- lapply(props.init , contCorrection) %>% unlist()
    props <- log(props.corrected)
    
  }else if (logtransform == FALSE){
    cat('\n calculating raw proportions \n')
    props <- props.init
    
  }else{
    stop('specify whether a log transform should be applied by stating logtransfom = TRUE or logtransfom = FALSE')
  }
  
  return(props)
}


#exploring assumption 1: assessing distirbution of proportions
#input is study based proportions
assumption1 <- function(num){
  fitted.norm = fitdist(num , distr = 'norm')
  summary(fitted.norm)
  plot(fitted.norm)
  gofstat(fitted.norm)
}

      
#exploring assumption 2: evaluating independence of covariates
#input main dataframe
assumption2 <- function(df, combinations, simulate.p.value = FALSE){
  nvar <- length(combinations)/2
  chiseq.val <- list()
  
  for (i in 1:nvar){
    
    a = paste0("df", '$' , combinations[1,i]) %>% parse(text=.) %>% eval()
    b = paste0("df", '$' , combinations[2,i])%>% parse(text=.) %>% eval()
    cont.tab <- table(a,b) %>% 
    chiseq.val[[i]] <- chisq.test(cont.tab , simulate.p.value = simulate.p.value)
  }

  return(chiseq.val)

}


#main
main <- function(data, splitby = NULL, logtransformprops = TRUE, covars, simulate.p.value = FALSE){
  if (class(data) != "data.frame"){
    stop('input is not a data frame')
    
  }else{
  data = data
  }
  
  #calculate proportions
  props <- formatDF(data) %>%
    groupDF(., covar = splitby) %>%
    calcProps(., logtransform = logtransformprops)

  #test assumption 1: is our measure variable parametric?
  assumption1(props) %>% print()
  
  #get covar combinations
  combos <- combn(covars , m = 2)

  #testing assumption 2: independence of covariates
  assumption2(data, combinations = combos , simulate.p.value = simulate.p.value)
  cat("\nassumption evaluation complete \nDONE")

}  


#import dataset
data_master<- read.csv("data/data_master.csv", na.strings = "NA")

##RUN##
#potential covars: "participant.seropositivity","grouped.method","reported.exposure","grouped.subtype" 
covars <- c("participant.seropositivity","grouped.method","reported.exposure")
main(data_master , covars = "reported.exposure")

##END##