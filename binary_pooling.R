###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity
# calculates summary proportion of infections initiated by multiple founder variants
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 3. One-step binomial GLMM allowing for clustering by study. random effects between studies
#    (random intercepts). approx ML fit
# 4. Two-step beta-binomial GLMM, dispersion param for study labels. Laplace approximate ML estimation

# estimations of mean effect size, confidence intervals and heterogeneity are presented 
# Sensitivity analyses are conducted as follows:
# SA1: leave-one-out cross validation of all models
# SA2: exclude studies with fewer than 10 participants
# SA3: excluding studies that only report single founder infectioins (ie prop multifounder =0)
# SA4: exclude studies published before 2008 (ie studies that predate SGA)
# SA6: bootstrap/resampling of participants for whom multiple measurements are available

###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(mltools)
library(data.table)
library(metafor)
library(aod)
library(ggplot2)
library(influence.ME)
library(kableExtra)
library(reshape2)
source('~/foundervariantsHIV_sysreview/generalpurpose_funcs.R')


# Add increment to studies with only single founder infections (for one-step models)
AddIncr <- function(df, incr){
  split_df <- split.data.frame(df , df$publication)
  
  list_incr = list()
  for (i in 1:length(split_df)){
    data = split_df[[i]]
    if(!(1 %in% data$multiple.founders)){
      data$multiple.founders = incr
    }else{
      data = data}
    list_incr[[i]] = data
  }
  df_incr = do.call(rbind.data.frame,list_incr)
  stopifnot(nrow(df) == nrow(df_incr))
  
  return(df_incr)
}


# Two-step binomial/normal model, pooling studies using Inverse Variance method,
# random effects, REML estimator of tau2.
CalcTwostepBiNorm <- function(data){
  step1 <- escalc(xi = multiplefounders ,
                  ni = subjects , 
                  data= data ,
                  drop00 = FALSE,
                  add = 0.0005, 
                  to = "only0",
                  measure = "PLO")
  
  step2 <- rma.uni(yi, vi, 
                   data = step1,
                   method = "REML",
                   knha = TRUE,
                   measure = "PLO")
  
  list(step1, step2) %>%
    return()
}


# One-step GLMM accounting for clustering of studies using a random intercept (
# random slope, random & uncorrelated intercept)
CalcOnestepBiRand <- function(data){
  model <- rma.glmm(xi = multiplefounders,
                    ni = subjects, 
                    data = data , 
                    drop00 = FALSE, 
                    add = 0.0005, 
                    measure= 'PLO', 
                    nAGQ = 7)
  return(model)
}


# Extracts estimates of summary effect from models
CalcEstimates <- function(model , analysis = "original", mermod.method = 'Profile'){
  if (class(model) == "rma" || class(model) == "rma.uni" || class(model) == "rma.glmm"){
    beta <- model$beta
    ci.lb <- model$ci.lb
    ci.ub <- model$ci.ub
  }
  else if (class(model) =="glmerMod"){
    beta <- summary(model)$coefficients[1,1]
    ci <- confint.merMod(model, method = mermod.method)
    ci.lb <- ci[nrow(ci),1]
    ci.ub <- ci[nrow(ci),2]
  }
  else if (class(model) =="metaprop"){
    beta <- model$TE.random
    ci.lb <- model$lower.random
    ci.ub <- model$upper.random
  }
  else if (class(model) =="glimML"){
    beta <- model@param[1]
    ci <- varbin(subjects,multiplefounders, data = model@data)@tab[5,c(3,4)]
    ci.lb <- as.numeric(ci[1]) %>% transf.logit()
    ci.ub <- as.numeric(ci[2]) %>% transf.logit()
  }

  results <- cbind.data.frame("estimate" = beta,
                              "estimate.lb" = ci.lb,
                              "estimate.ub" = ci.ub) %>%
    mapply(transf.ilogit, ., SIMPLIFY = FALSE)
  
  results_lab <- cbind.data.frame("model" = substitute(model) %>% 
                                    deparse() %>% 
                                    gsub("[:.:].*", "", .),
                                  "analysis" = analysis,
                                  results)
  
  return(results_lab)
}


# Extracts estimates of heterogeneity (tau2) from models
CalcHet <- function(model , analysis = "original"){
  if (class(model)[1] == "rma.uni"){
    tau2 <- model$tau2 
    q <- model$QE
    i2 <- model$I2
    phi <- NA
  }
  else   if (class(model)[1] == "rma.glmm"){
    tau2 <- model$tau2 
    q <- model$QE.Wld
    i2 <- model$I2
    phi <- NA
  }else   if (class(model)[1] == "metaprop"){
    tau2 <- model$tau2$TE
    q <- model$Q
    i2 <- model$I2
    phi <- NA
  }
  else if (class(model) =="glmerMod"){
    tau2 <- VarCorr(model)[[1]][1]
    q <- NA
    i2 <- NA
    phi <- NA
  }
  else if (class(model) =="glimML"){
    tau2 <- NA
    q <- NA
    i2 <- NA
    phi <- model@param[2] %>% as.numeric()
  }
  
  results <- cbind.data.frame("model" = substitute(model) %>% 
                                deparse() %>% 
                                gsub("[:.:].*", "", .),
                              "tau2" = tau2,
                              "Q" = q,
                              "I2" = i2,
                              "phi" = phi,
                              "analysis" = analysis)
  
  return(results)
}


# Create list of dataframes for leave-one-out cross validation
LOOCV.dat <- function(data){
  pubs <- unique(data$publication)
  loo <- list()
  loo.pubs <- list()
  
  for (i in pubs){
    loo[[i]] <- data[data$publication != i, ]
    loo.pubs[[i]] <- pubs[pubs != i]
  }
  out <- list(loo, loo.pubs)
  stopifnot(length(loo) == length(loo.pubs))
  return(out)
}


# Extract estimates from LOOCV to create dataframe (input for influence plot)
DFInfluence <- function(model,labs){
  
  names <- paste("Omitting" , labs %>% names(), sep = " ") %>% as.factor()
  beta = numeric()
  ci.lb = numeric()
  ci.ub = numeric() 
  
  if (class(model[[1]]) == "rma" || class(model[[1]]) == "rma.uni" || class(model[[1]]) == "rma.glmm"){
    for (i in 1:length(model)){
      beta[i] <- model[[i]]$beta
      ci.lb[i] <-model[[i]]$ci.lb
      ci.ub[i] <- model[[i]]$ci.ub
    }
  }else if (class(model[[1]]) =="glmerMod"){
    for (i in 1:length(model)){
      ci <- confint(model[[i]])
      beta[i] <- summary(model[[i]])$coefficients[1,1]
      ci.lb[i] <-ci[nrow(ci),1]
      ci.ub[i] <- ci[nrow(ci),2]
    }
    }else if (class(model) =="metaprop"){
      beta <- model$TE.random
      ci.lb <- model$lower.random
      ci.ub <- model$upper.random
    
  }else if (class(model[[1]]) =="glimML"){
    #Bootstrapped Binomial CIs-check Chuang-Stein 1993
    for (i in 1:length(model)){
      binom.ci <- varbin(subjects,multiplefounders, data = model[[i]]@data)@tab[5,c(3,4)] %>%
        as.numeric() %>%
        transf.logit()
      beta[i] <- model[[i]]@param[1]
      ci.lb[i] <- binom.ci[1]
      ci.ub[i] <- binom.ci[2]
    }
  }else{
    stop('no valid model detected.')
  }
  
  influence_out <- cbind.data.frame('trial'= names,
                                    "estimate" = transf.ilogit(beta),
                                    "ci.lb" = transf.ilogit(ci.lb),
                                    "ci.ub"= transf.ilogit(ci.ub)) 
  

  return(influence_out)
}

# Generate resampled datasets and calculate model estimates for psuedo-bootstrap 
# sensitivity analysis of inclusion/exclusion criteria
BootParticipant <- function(data, replicates){
  require(parallel)
  require(lme4)
  require(metafor)
  require(dplyr)
  
  resampled <- lapply(1:replicates, function(x,y) {y %>% group_by(participant.id_) %>% slice_sample(n=1)},
                      y = data)
  
  resampled_props <- lapply(resampled , CalcProps)
  
  cl <- detectCores() %>%
    `-` (2) %>%
    makeCluster()
  
  clusterEvalQ(cl, c(library(lme4), library(metafor), library(aod), library(dplyr),
                     set.seed(4472),'CalcOnestepBiRand', 'CalcTwostepBiNorm'))
  
  start <- Sys.time()
  print(start)
  
  twostep_boot <- parLapply(cl = cl, resampled_props, CalcTwostepBiNorm)
  rand_boot <- parLapply(cl = cl, resampled_props, CalcOnestepBiRand)
  
  end <- Sys.time()
  elapsed <- end-start
  print(elapsed)
  
  stopCluster(cl)
  remove(cl)
  
  twostep_boot.est <- lapply(twostep_boot, function(mod) mod[[2]]$beta) %>%
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  twostep_boot.het <- lapply(twostep_boot, function(mod) CalcHet(mod[[2]], analysis = "twostep_binorm")) %>%
    do.call(rbind.data.frame,.)
  
  rand_boot.est <- lapply(rand_boot, function(mod) mod$beta) %>%
    do.call(rbind.data.frame,.) %>%
    {cbind.data.frame("estimate"=transf.ilogit(.[,1]))}
  
  rand_boot.het <- lapply(rand_boot, function(mod) CalcHet(mod, analysis = "onestep_bi_rand")) %>%
    do.call(rbind.data.frame,.)
  
  boot_estimates <- rbind.data.frame(twostep_boot.est,
                                     rand_boot.est)
  
  boot_het <- rbind.data.frame(twostep_boot.het,
                               rand_boot.het)
  
  out <- cbind.data.frame(boot_estimates, boot_het) %>% .[,-2]

  
  return(out)
}


###################################################################################################
###################################################################################################

# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF(., noreps = TRUE)
df_props <- CalcProps(df)  
publist <- df %>%
  pull(.,var=publication_) %>%
  unique()

# Set seed
set.seed(4472)


###################################################################################################
# 1. Two-step binomial-normal model
# step 1: pooling within studies using binomial model/logit 
# step 2: pooling across studies using random effects (normal model). Inverse variance method used to pool. 
# REML estimator of Tau. 

twostep_binorm.all <- CalcTwostepBiNorm(df_props)
twostep_binorm.step1 <- twostep_binorm.all[[1]]
twostep_binorm <- twostep_binorm.all[[2]]

twostep_binorm.sum <- summary(twostep_binorm)
twostep_binorm.sum

twostep_binorm.est <- CalcEstimates(twostep_binorm)
twostep_binorm.het <- CalcHet(twostep_binorm)
###################################################################################################
###################################################################################################

# 2. One-step binomial GLMM allowing for clustering by study. 
# uncorrelated random intercept and random slope within group 
# approx ML fit
# assumes conditional independence and follow binomial distribution

onestep_bi_rand <- CalcOnestepBiRand(df_props) 
onestep_bi_rand.sum <- summary(onestep_bi_rand)
onestep_bi_rand.sum

onestep_bi_rand.est <- CalcEstimates(onestep_bi_rand)
onestep_bi_rand.het <- CalcHet(onestep_bi_rand)


###################################################################################################
###################################################################################################
# Model comparison: Estimated sumary effects (prop, CI), between study variance (tau, I^)
# Tau2 = Var(theta_i), theta_i = E[theta_i]
  estimates <- rbind.data.frame(twostep_binorm.est,
                                onestep_bi_rand.est)


heterogeneity <- rbind.data.frame(twostep_binorm.het,
                                  onestep_bi_rand.het)

  
###################################################################################################
###################################################################################################
# Sensitivity Analyses
# SA1. Influence of Individual Studies
# SA2. Exclusion of small sample sizes (less than n = 10)
# SA3. Exclusion of studies with 0 multiple founder variants
# SA4. Exclusion of all studies that do not use SGA
# SA5. Resampling of participants for which we have multiple measurments (takes pre-formatted DF)


# SA1. Influence of Individual Studies (LOOCV)
df_loocv <- LOOCV.dat(df)[[1]]
publist_loocv <- LOOCV.dat(df)[[2]]
dfp_loocv <- LOOCV.dat(df_props)[[1]]

twostep_binorm.influence <- lapply(dfp_loocv, function(x) CalcTwostepBiNorm(data = x)[[2]]) %>%
  DFInfluence(., labs = publist_loocv) %>% {cbind.data.frame(.,'model' = 'twostep_binorm')}

onestep_bi_rand.influence <-  lapply(dfp_loocv ,CalcOnestepBiRand) %>%
  DFInfluence(., labs =publist_loocv) %>% {cbind.data.frame(.,'model' = 'onestep_bi_rand')}

influence_df <- rbind.data.frame(twostep_binorm.influence,
                                 onestep_bi_rand.influence)


# SA2. Exclusion of small sample sizes (less than n = 10)
publist.nosmallsample <- subset(df_props , subjects > 9 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nosmallsample <- df[df$publication_ %in% publist.nosmallsample,]

df_props.nosmallsample <- df_props[df_props$publication %in% publist.nosmallsample,]

twostep_binorm.nosmallsample <- CalcTwostepBiNorm(df_props.nosmallsample)[[2]]
twostep_binorm.nosmallsample.out <- list(CalcEstimates(twostep_binorm.nosmallsample, analysis = "no_small"),
                                         CalcHet(twostep_binorm.nosmallsample)) %>%
  cbind.data.frame(.)

onestep_bi_rand.nosmallsample <- CalcOnestepBiRand(df_props.nosmallsample)
onestep_bi_rand.nosmallsample.out <- list(CalcEstimates(onestep_bi_rand.nosmallsample, analysis = "no_small"),
                                          CalcHet(onestep_bi_rand.nosmallsample)) %>%
  cbind.data.frame(.)


SA2_results <-  rbind.data.frame(twostep_binorm.nosmallsample.out,
                                 onestep_bi_rand.nosmallsample.out)


# SA3. Exclusion of studies with 0 multiple founder variants 
publist.nozeros <- subset(df_props , multiplefounders != 0 , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.nozeros <- df[df$publication_ %in% publist.nozeros,]

df_props.nozeros <- df_props[df_props$publication %in% publist.nozeros,]

twostep_binorm.nozeros <- CalcTwostepBiNorm(df_props.nozeros)[[2]]
twostep_binorm.nozeros.out <- list(CalcEstimates(twostep_binorm.nozeros, analysis = "no_zeros"),
                                   CalcHet(twostep_binorm.nozeros)) %>%
  cbind.data.frame(.)

onestep_bi_rand.nozeros <- CalcOnestepBiRand(df_props.nozeros)
onestep_bi_rand.nozeros.out <- list(CalcEstimates(onestep_bi_rand.nozeros, analysis = "no_zeros"),
                                    CalcHet(onestep_bi_rand.nozeros)) %>%
  cbind.data.frame(.)


SA3_results <- rbind.data.frame(twostep_binorm.nozeros.out, 
                                onestep_bi_rand.nozeros.out)

# SA4. Exclusion of all studies that do not use SGA
publist.sgaonly <- subset(df , sample.amplification_ == 'SGA' , select = publication_) %>%
  pull(.,var=publication_) %>%
  unique()

df.sgaonly <- df[df$publication_ %in% publist.sgaonly,]

df_props.sgaonly <- df_props[df_props$publication_ %in% publist.sgaonly,]

twostep_binorm.sgaonly <- CalcTwostepBiNorm(df_props.sgaonly)[[2]]
twostep_binorm.sgaonly.out <- list(CalcEstimates(twostep_binorm.sgaonly, analysis = "sga_only"),
                                   CalcHet(twostep_binorm.sgaonly)) %>%
  cbind.data.frame(.)

onestep_bi_rand.sgaonly <- CalcOnestepBiRand(df_props.sgaonly)
onestep_bi_rand.sgaonly.out <- list(CalcEstimates(onestep_bi_rand.sgaonly, analysis = "sga_only"),
                                    CalcHet(onestep_bi_rand.sgaonly)) %>%
  cbind.data.frame(.)


SA4_results <- rbind.data.frame(twostep_binorm.sgaonly.out, 
                                onestep_bi_rand.sgaonly.out)

# SA5. Resampling of participants for which we have multiple measurments (aim is to generate a distribution of possible answers)
resampling_df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF(., noreps = FALSE)

boot_participant <- BootParticipant(resampling_df , 1000)


###################################################################################################
###################################################################################################
# Outputs
# CSV of pooling model results (estimates only) with SAs 2 + 3
originals <- cbind.data.frame(estimates, heterogeneity)
pooled_est <- rbind.data.frame(originals,
                               SA2_results,
                               SA3_results,
                               SA4_results) %>% .[,-c(6,11)]

write.csv(pooled_est , file = 'bp_estsa2sa3sa4.csv', row.names = F)


# CSV study influence
write.csv(influence_df, file = 'bp_sa1.csv',row.names = F)

# CSV resampling (pseudo bootstrapping)
write.csv(boot_participant, file = 'bp_resampl.csv',row.names = F)

# Forest Plot 1-step BN
jpeg("testplot.jpeg" ,width = 5250, height = 6500, res = 380 ,units = "px", pointsize = 12)
meta::metaprop(data = df_props,
               n = subjects,
               event = multiplefounders,
               studlab = sort(publist) %>% gsub("_" , " " , .),
               method = 'GLMM',
               sm = 'PLOGIT',
               incr = 0.0005,
               allincr = F,
               comb.fixed = F,
               method.tau = 'ML') %>% meta::forest(hetstat = 'random')


dev.off()


# Funnel Plot
funnel_data <- cbind.data.frame('se' = sqrt(onestep_bi_rand.nozeros$vi), 'b' =  onestep_bi_rand.nozeros$yi)
lb <- onestep_bi_rand.nozeros$ci.lb 
ub <- onestep_bi_rand.nozeros$ci.ub 
u <- mean(funnel_data$b)
se <- onestep_bi_rand.nozeros$se

poldgpn <- data.frame(x=c(-3.8,u,1.9), y = c(1.5,0,1.5))
plt <- ggplot(funnel_data ) +
  geom_polygon(aes(x=x, y = y), data =  poldgpn ,fill = 'white', linetype = 'dashed' , color = 'black')  +
  geom_point( aes(y = se, x = b, colour = ), shape = 4, size = 3)+
  theme_classic() +
  scale_x_continuous(limits = c(-5 , 3), expand = c(0,0), name = 'Log Odds of Multiple Founders')+
  scale_y_reverse(limit=c(1.5,0),  expand = c(0,0), name = 'Standard Error') +
  
  geom_segment(aes(x=u, y =1.5, xend = u, yend=0)) +
  theme(panel.background = element_rect(fill = 'gray97' )) +
  scale_color_npg()
###################################################################################################
###################################################################################################
# END #
###################################################################################################
###################################################################################################