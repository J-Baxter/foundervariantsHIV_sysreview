###################################################################################################
###################################################################################################
# General purpose functions for Inferring the multiplicity of founder variants initiating HIV-1 
# infection: a systematic review and IPD meta-analysis
###################################################################################################
###################################################################################################

# Formats data spreadsheet for analysis. Removes duplicates and NAs.
formatDF <-  function(df){
  require(tidyr)
  require(reshape2)
  #create dummy variables in founder multiplicity col
  if (class(df$multiple.founders)=="factor"){
    df$multiple.founders = 1 - (as.numeric(df$multiple.founders)-1)
  }else{
    stop('founder multiplicity is not a factor')
  }
  df_nona <- df[!is.na(df$multiple.founders),]
  df_nodups <- df_nona[(df_nona$include.main == '') & (df_nona$exclude.repeatstudy == ''),]
  df_labelled <- unite(df_nodups, "publication", c(author ,year), sep = '_')
  df_splittrans <- colsplit(df_labelled$reported.exposure, ":" , c("riskgroup" , "direction")) %>% cbind.data.frame(.,df_labelled)
  return(df_splittrans)
}


###################################################################################################
# Sums number of patients within each study and number of infections initiated by multiple founders. 
# Can be stratified with additional covariates
CalcProps <- function(.data, ...){
  .data %>% 
    group_by(publication, ...) %>%
    summarise(subjects = n(), multiplefounders = sum(multiple.founders)) %>% 
    as.data.frame()
  
}


###################################################################################################
###################################################################################################