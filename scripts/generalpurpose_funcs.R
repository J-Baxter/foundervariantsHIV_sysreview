###################################################################################################
###################################################################################################
# General purpose functions for Inferring the multiplicity of founder variants initiating HIV-1 
# infection: a systematic review and IPD meta-analysis
###################################################################################################
###################################################################################################

# unzip file or directory containing data csv
Retrieve <- function(zipfile){
  if(file.exists(zipfile) & grepl('zip', zipfile)){
    newname <- gsub('.zip', '', zipfile)
    
    if (!dir.exists(newname)){
      dir.create(newname)
    }
    else{
      stop('data directory already exists')
    }
    
    
    unzip(zipfile, exdir= newname)
    filepath <- file.path(newname)
    filelist <- list.files(path = newname, recursive = T)
    file.copy(paste0('data/',filelist), filepath)
    dirlist <- list.dirs(path = filepath, recursive = T)
    unlink(dirlist[-1], recursive = T)
  }
  else{
    stop('not a .zip file')
  }
}


# Formats data spreadsheet for analysis. Removes duplicates and NAs.
# noreps removes repeated measurements (reason for false would be to conduct sensitivity analyses)
# filter removes factors containing less than 5 individuals from specified df
formatDF <-  function(df, noreps = TRUE, filter = NULL){
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
    df_labelled <- unite(df_nodups, "publication", c(author, year), sep = '_', remove = FALSE) %>%
      dplyr::select(-author)  
    
  }else{
    df_labelled <- unite(df_nona, "publication", c(author, year), sep = '_', remove = FALSE) %>%
      dplyr::select(-author)
  }
  
  if(is.character(filter)){
    looped_df <- df_labelled
    for (i in 1:length(filter)){
      fac <- filter[i]
      tbl <- table(looped_df[fac])
      looped_df <- looped_df[!looped_df[[fac]] %in% names(tbl)[tbl < 6],]
    }
    
    
    df_splittrans <- reshape2::colsplit(looped_df$reported.exposure, ":" , c("riskgroup" , "direction")) %>%
      type.convert(., as.is = FALSE) %>%
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
      dplyr::summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>% 
      as.data.frame()
  }else{
    .data %>% 
      group_by(publication_) %>%
      dplyr::summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>% 
      as.data.frame()
  }

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
# One-step GLMM accounting for clustering of studies using a random intercept
CalcRandMetaReg <- function(data, formula, opt = NULL){
  
  options(warn = 1)
  f <- as.formula(formula)
  environment(f) <- environment()
  
  if(is.character(opt)){
    model <- lme4::glmer(f,
                         data = data,
                         family = binomial(link = "logit"),
                         nAGQ = 1,
                         control = glmerControl(optCtrl = list(maxfun = 50000000),
                                                check.nobs.vs.nlev = 'ignore',
                                                check.nobs.vs.nRE = 'ignore',
                                                optimizer = 'bobyqa'))

  }else{
    model <- lme4::glmer(f,
                         data = data,
                         family = binomial(link = "logit"),
                         nAGQ = 1,
                         control = glmerControl(optCtrl = list(maxfun = 50000000),
                                                check.nobs.vs.nlev = 'ignore',
                                                check.nobs.vs.nRE = 'ignore'))
  }
  
 
  return(model)
}


###################################################################################################
# Execute a list of lme4 models in parallel
RunParallel <- function(func, v1, v2, ...){9
  options(warn = 1)
  
  # Set up cluster (fork)
  cl <- detectCores() %>%
    `-` (2) 
  
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

    remove(cl)
    
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
# Extracts covariate names from lmer function syntax
GetName <- function(x, effects = NULL) {
  require(stringr)
  
  if(is.null(effects)){
    print('requires user to specify whether fixed or random effects are required')
  }
  
  if(effects == 'fixed'){
    name <-gsub(".*[:~:] (.+?) [:(:].*", "\\1", x) %>%
      gsub("[:+:]([:^+:]*)$","",.) %>%
      gsub("[:.:]" , " " , .) %>% 
      gsub("[:_:]" , "" , .) %>%
      str_to_title()%>%
      str_trim()
    
  }else if( effects == 'random'){
    name <-gsub("^[^\\(]+", "\\1", x, perl = T) %>%
      gsub("_" , "" , .) %>%
      gsub(":" , " : " , .) %>%
      str_to_title()%>%
      str_trim()
  }
  
  return(name)
}


###################################################################################################
# Extracts covariate names from lmer function syntax
# Pipeline for check_singularity, check_convergence and logloss functions
CheckModels <- function(modellist){
  require(performance)
  
  if (class(modellist) == 'list'){
    is.sing <- lapply(modellist, check_singularity, tolerance = 1e-05) %>% do.call(rbind.data.frame, .)
    is.con <- lapply(modellist, check_convergence) %>% do.call(rbind.data.frame, .)
    
  }else{
    is.sing <- check_singularity(modellist, tolerance = 1e-05) 
    is.con <- check_convergence(modellist)
  }
  
  
  out <- cbind.data.frame(is.sing, is.con) %>% 
    `colnames<-` (c('is.singular', 'is.converged'))
  
  rownames(out) <- names(modellist)
  
  return(out)
}


###################################################################################################
# Extract intercept, fixed effects and random effects from the models
# Output is a list of dataframes
GetCoefs <- function(model, label = "original"){
  # Calculate CIs
  options(warn = 1)
  
  n_cores <- detectCores() %>% `-` (2)
  
  ci <- lme4::confint.merMod(model, 
                             method = 'boot',
                             FUN = NULL,
                             #.progress="txt", 
                             PBargs=list(style=3), 
                             nsim = 500,
                             parallel= "multicore",
                             ncpus = n_cores
  )
  
  re.num <- ranef(model) %>% length()
  
  # Extract Fixed Effects
  if (length(fixef(model)) > 1){
    
    fe <- fixef(model) 
    sd <- sqrt(diag(vcov(model)))
    ci.fe <- ci[-c(1,re.num),]
    nom <- names(fe)
    sum <- summary(model)$coefficients
    
    fix_df <- cbind.data.frame(nom = nom,
                               est = fe,
                               sd = sd,
                               ci.lb = ci.fe[,1],
                               ci.ub = ci.fe[,2],
                               z.val = sum[,3],
                               p.val = sum[,4],
                               analysis = label) %>% 
      `row.names<-` (NULL) %>%
      separate(nom , c('covariate' , 'level') , '_')
    
  }else{
    fe <- fixef(model) 
    ci.fe <- ci[re.num,]
    nom <- names(fe)
    sd <-  NA
    sum <- summary(model)$coefficients
    
    fix_df <- cbind.data.frame(nom = nom,
                               est = fe,
                               sd = sd,
                               ci.lb = ci.fe[1],
                               ci.ub = ci.fe[2],
                               z.val = sum[,3],
                               p.val = sum[,4],
                               analysis = label) %>% 
      `row.names<-` (NULL) %>%
      separate(nom , c('covariate' , 'level') , '_')
  }
  
  int_df <- fix_df[which(is.na(fix_df$level)),]
  fix_df <- fix_df[which(!is.na(fix_df$level)),]
  
  # Extract Random Effects
  if (length(ranef(model)) > 1){
    re <- ranef(model)
    re.mean <- lapply(re, function(x) mean(x$`(Intercept)`)) %>%
      do.call(rbind.data.frame,.) %>% `colnames<-` ('mean')
    
    re.sd <- VarCorr(model) %>% as.data.frame()
    
    ci.re <- ci[1:re.num,]
    
    re_df <- cbind.data.frame(groups = gsub('_', '', re.sd[,1]),
                              mean = re.mean,
                              vcov = re.sd[,4],
                              sd = re.sd[,5],
                              ci.lb = ci.re[,1],
                              ci.ub = ci.re[,2],
                              analysis = label) %>% 
      `row.names<-` (NULL)
  }else{
    re <- ranef(model)
    re.mean <- lapply(re, function(x) mean(x$`(Intercept)`)) %>%
      do.call(rbind.data.frame,.) %>% `colnames<-` ('mean')
    
    re.sd <- VarCorr(model) %>% as.data.frame()
    
    ci.re <- ci[1:re.num,]
    
    re_df <- cbind.data.frame(groups = gsub('_', '', re.sd[,1]),
                              mean = re.mean,
                              vcov = re.sd[,4],
                              sd = re.sd[,5],
                              ci.lb = ci.re[1],
                              ci.ub = ci.re[2],
                              analysis = label) %>% 
      `row.names<-` (NULL)
  }
  
  
  out <- list(int_df, fix_df, re_df) %>% `names<-` (c('int' , 'fe', 're'))
  
  return(out)
}


###################################################################################################
# Extract Marginal effects and prediction/confidence intervals
GetEMM <- function(model, byvar, label = "original"){
  #if byvar is formula:
  if (grepl('~', byvar)){
    specs <- gsub(".*[:~:] (.+?) [:(:].*", "\\1", byvar) %>%
      gsub("[:+:]([:^+:]*)$","",.) %>%
      str_trim()%>%
      paste('~',., sep = ' ') %>%
      as.formula()
  }
  
  else{
    specs <- paste('~', byvar, sep = ' ') %>%
      str_trim()%>%
      as.formula()
    
  }
  
  environment(specs) <- environment()
  
  out <- emmeans(model, specs = specs, type = 'source') %>%
    {cbind.data.frame(.,label = label)}
  colnames(out)[1] <- 'covariate_level'
  return(out)
}

# Extract predictions and calculate bootstraped prediction/confidence intervals
GetPreds <- function(model, byvar, label = "original"){
  require(ggeffects)
  #if byvar is formula:
  if (grepl('~', byvar)){
    terms <- gsub(".*[:~:] (.+?) [:(:].*", "\\1", byvar) %>%
      gsub("[:+:]([:^+:]*)$","",.) %>%
      str_trim()
  }
  
  else{
    terms <- str_trim(byvar)
    
  }
  
  out <- ggpredict(model, terms = terms, type = 'fe') %>%
    {cbind.data.frame(.,label = label)}
  colnames(out)[1] <- 'covariate_level'
  return(out)
}

###################################################################################################
Effects2File <- function(effectslist){
  int <- lapply(effectslist, function(x) x$int) %>% do.call(rbind.data.frame, .)
  fe <- lapply(effectslist, function(x) x$fe) %>% do.call(rbind.data.frame, .)
  re <- lapply(effectslist, function(x) x$re) %>% do.call(rbind.data.frame, .)
  
  out <- list(int, fe, re)
  names(out) <- c('int', 'fe', 're')
  return(out)
}


###################################################################################################
###################################################################################################