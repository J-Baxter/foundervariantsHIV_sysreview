###################################################################################################
###################################################################################################
# General purpose functions for Inferring the multiplicity of founder variants initiating HIV-1 
# infection: a systematic review and IPD meta-analysis
###################################################################################################
###################################################################################################

# Formats data spreadsheet for analysis. Removes duplicates and NAs.
# noreps removes repeated measurements (reason for false would be to conduct sensitivity analyses)
# filter removes factors containing less than 5 individuals from specified df
formatDF <-  function(df, noreps = TRUE, filter = NULL){
  require(tidyr)
  require(reshape2)
  
  #create dummy variables in founder multiplicity col
  if (class(df$multiple.founders)=="factor"){
    df$multiple.founders = 1 - (as.numeric(df$multiple.founders)-1)
    
  }else{
    stop('founder multiplicity is not a factor')
  }
  
  df_nona <- df[!is.na(df$multiple.founders),]
  
  if (noreps == TRUE){
    df_nodups <- df_nona[(df_nona$include.main == '') & (df_nona$exclude.repeatstudy == ''),]
    df_labelled <- unite(df_nodups, "publication", c(author ,year), sep = '_', remove = FALSE) %>%
      select(-author)
    
  }else{
    df_labelled <- unite(df_nona, "publication", c(author ,year), sep = '_', remove = FALSE) %>%
      select(-author)
  }
  
  if(is.character(filter)){
    looped_df <- df_labelled
    for (i in 1:length(filter)){
      fac <- filter[i]
      tbl <- table(looped_df[fac])
      looped_df <- looped_df[!looped_df[[fac]] %in% names(tbl)[tbl < 6],]
    }
    
    
    df_splittrans <- colsplit(looped_df$reported.exposure, ":" , c("riskgroup" , "direction")) %>%
      type.convert() %>%
      cbind.data.frame(.,looped_df)
    #df_splittrans <- df_splittrans[!(df_splittrans$reported.exposure %in% "unknown.exposure"), ]
    
  }else{
    df_splittrans <- colsplit(df_labelled$reported.exposure, ":" , c("riskgroup" , "direction")) %>%
      cbind.data.frame(.,df_labelled)
  }

  colnames(df_splittrans) <- paste0(colnames(df_splittrans), "_")

  return(df_splittrans)
}


###################################################################################################
# Sums number of patients within each study and number of infections initiated by multiple founders. 
# Can be stratified with additional covariates
CalcProps <- function(.data, ...){
  if (nargs()>1){
    .data %>% 
      group_by(publication_, !!sym(...)) %>%
      summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>% 
      as.data.frame()
  }else{
    .data %>% 
      group_by(publication_) %>%
      summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>% 
      as.data.frame()
  }

}


###################################################################################################
# Does what is says on the tin
# Assumes sampling distribution of parameters is multivariate normal
CalcCI <- function(u,se,threshold){
  value <- 1-(threshold/2)
  upper <- u + se*qnorm(value)
  lower <- u - se*qnorm(value)
  ci <- c(lower,upper)
  return(ci)
}


###################################################################################################
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
      
      print('Baseline Covariates:')
      print(levels(dataframe[,covar.])[1])
      
    }else{
      warning('requires factor as input.')
    }
  }
  
  return(dataframe)
}


###################################################################################################
# Execute a list of lme4 models in parallel
RunParallel <- function(func, v1, v2, ...){
  options(warn = 1)
  
  # Set up cluster (fork)
  cl <- detectCores() %>% `-` (2) 
  
  if (class(v2) == 'data.frame'){
    
    start <- Sys.time()
    start
    
    para <- mclapply(v1,
                     func, 
                     data = v2,
                     ...,
                     mc.cores = cl,
                     mc.set.seed = FALSE) #child process has the same initial random number generator (RNG) state as the current R session
    
    end <- Sys.time()
    elapsed <- end-start
    print(elapsed)
  }else if(class(v2) != 'data.frame'){
    
    start <- Sys.time()
    start
    para <- mcmapply(func,
                     v1,
                     v2,
                     ...,
                     mc.cores = cl,
                     mc.set.seed = FALSE,
                     SIMPLIFY = F) 
    
    end <- Sys.time()
    elapsed <- end-start
    print(elapsed)
  }
  
  
  return(para)
}


###################################################################################################
###################################################################################################