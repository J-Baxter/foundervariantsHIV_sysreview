#test script for meta analysis models using keele, abrahams, haaland and li papers

formatDF <-  function(df){
  df_nona <- df[!is.na(df$multiple.founders),]
  df_nodups <- df_nona[(df_nona$include.main == '') & (df_nona$exclude.repeatstudy == ''),]
  df_labelled <- unite(df_nodups, "publication", c(author ,year), sep = '_')
  return(df_labelled)
}


CalcProps <- function(df, covar = NULL){
  #create dummy variables in founder multiplicity col
  if (class(df_labelled$multiple.founders)=="factor"){
    df_labelled$multiple.founders = 1 - (as.numeric(df_labelled$multiple.founders)-1)
  }else{
    stop('founder multiplicity is not a factor')
  }
  
  #summarise
  if (is.null(covar)){
    df_grouped <- df_labelled %>% 
      group_by(publication) %>%
      summarise(subjects = n(), multiplefounders = sum(multiple.founders))
  }else{
    df_grouped <- df_labelled %>% 
      group_by(publication, sym(covar)) %>%
      summarise(subjects = n(), multiplefounders = sum(multiple.founders))}
  
  return(df_grouped)
}


step1 <- function(data, study_id){
  df_subset <- subset(data, publication==study_id)
  events <- nrow(df_subset)
  success <-  sum(df_subset$multiple.founders)
  prop.multifounders <- success/events
  incr = 0.0005 
  incr.event = incr/events
  
  if (prop.multifounders == 0){
    lr <- glm(multiple.founders+incr.event ~ 1, df_subset , family = binomial(link = "logit"))
    log_or <- summary(lr)$coefficients[,1]
    se <- summary(lr)$coefficients[, 2]
    
  }else{
    lr <- glm(multiple.founders ~ 1, df_subset , family = binomial(link = "logit"))
    log_or <- summary(lr)$coefficients[,1]
    se <- summary(lr)$coefficients[, 2]
  }
  
  agg.results <- cbind.data.frame(study_id, log_or, se)
  
  return(agg.results)
}

#define main
main <- function(){
  #import data
  df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()
  
  #set test data
  testlist <- c('Keele_2008' , "Abrahams_2009", "Haaland_2009", "Li_2010")
  testset_df <- lapply(testlist, function(x,y) subset(x, publication == y), x = df) %>% do.call(rbind.data.frame,.)
  
 
  #pooling (no evaluation of covariates)
  #twostep - binomial/normal model coded manually
  #step1 pooling within studies
  summary = lapply(testlist, step1, data = testset_df) %>% do.call(rbind.data.frame,.)
  
  #step 2a pooling across studies using random effects (normal model). 
  
  
}