###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity: Comparing placebo vs vaccine escape participants
# in HIV-1 vaccine trials

# calculates summary proportion of infections initiated by multiple founder variants
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. random effects between studies
#    (random intercepts). approx ML fit

# Then, a univaribale metaregression, segregating by Vaccine status
###################################################################################################
###################################################################################################
# RUN FROM HERE #
# Dependencies
source('./scripts/load_packages.R')
source('./scripts/generalpurpose_funcs.R')

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
                    slab = gsub('[_]', ', ' , publication_),
                    nAGQ = 1)
  return(model)
}


# Extracts estimates of summary effect from models
CalcEstimates <- function(model , analysis = "original"){
  
  # One and two step binomial models
  if (class(model) == "rma" || class(model) == "rma.uni" || class(model) == "rma.glmm"){
    beta <- model$beta
    ci.lb <- model$ci.lb
    ci.ub <- model$ci.ub
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

###################################################################################################
###################################################################################################
# Import data
if (!dir.exists('data')){
  Retrieve('data.zip')
}else{
  Sys.sleep(0.2)
}

df <- read.csv("./data/meta_analysis_data.csv",
               na.strings = "NA",
               stringsAsFactors = T) %>%
  formatDF(., filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  filter(!is.na(df$vaccine.status_ )) %>%
  droplevels()

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
# STAGE 2: Univariate meta-regression of individual covariates against founder variant multiplicity
# Set reference levels for meta regression
###################################################################################################
# Equivalent to a subgroup analysis with random effects for subgroup and cohort
form <- "multiple.founders_ ~ vaccine.status_  + (1 | publication_)"

unipooled_effectstruct <- GetName(form, effects = 'fixed')

# Run models
unipooled_model <- CalcRandMetaReg(df, form, opt = '1')

# Check model convergence and singularity
unipooled_check <- CheckModels(unipooled_models)%>% 
  `row.names<-`(unipooled_effectstruct)


###################################################################################################
# Extract fixed and random effect coefficients and calculate bootstrapped 95% CIs
# fixed effects coefficients exponentiated to odds ratios
unipooled_models.coef <- GetCoefs(unipooled_model, unipooled_effectstruct)


###################################################################################################
###################################################################################################
# Outputs
pooled_models <-  read.csv('./results/pooling_estsa2sa3sa4sa6sa7.csv')

models <- c('Two-Step Binomial-Normal',
            'One-Step Binomial GLMM')

estimates$analysis <- 'vaccine'

pooled <- rbind("model" = pooled_models[1:2,1:5], estimates) %>% arrange(., model)


# Export csv with pooling and univariable metaregression to file
originals <- cbind.data.frame(estimates, heterogeneity)

# Figure S9a - Pooled Original vs Pooled Vaccine Only 
# Figure S9b - Vaccine Subgroups comparison

#
figureS9a <- ggplot(pooled,
                  aes(x= forcats::fct_rev(model), y = estimate, color = analysis)) +
  
  geom_point( shape = 4, 
              size = 4,
              position = position_dodge(0.5)) + 
  
  scale_y_continuous(name = "Probability of Multiple Founders",
                     #limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  coord_cartesian(ylim = c(0,.5))+
  
  scale_x_discrete(name = "Model", 
                   labels = c(
                     onestep_bi_rand = "GLMM",
                     twostep_binorm = "B-N"
                   )) +
  theme_bw() + 
  
  coord_flip() +
  
  guides(colour = guide_legend(reverse=T))+
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_colour_npg(name = 'Analysis', labels = c(
    original = "Full analysis",
    vaccine = "Vaccine trial participants only")) + 
  
  theme(legend.position = c(0.8,0.86),
        axis.text = element_text(size = 9.5),
        legend.text = element_text(size = 9.5),
        axis.title = element_text(size = 11),
        legend.background = element_blank()#,
        #plot.margin = unit(c(2,4,2,1), "lines")
  )

vaccine_ref <- cbind.data.frame(level = 'placebo', est = 0, ci.lb = NA, ci.ub = NA)
vaccine_subgroup <- rbind.data.frame(unipooled_models.coef$fe[, c(2,3,5,6)], vaccine_ref)

figureS9b <- ggplot(vaccine_subgroup,
                    aes(x = level , y = exp(est))) +
  
  geom_point( shape = 18, 
              size = 4) + 
  
  scale_y_continuous(name = "Odds Ratio",
                     #limits=c(0,.5),
                     expand = c(0.01, 0.01),
                     limits = c(0, 4)) +
  
  scale_x_discrete(name = "Participant Subgroup", 
                   labels = c(
                     vaccine = "Vaccine",
                     placebo = "Placebo"
                   )) +
  theme_bw() + 
  
  coord_flip() +
  
  guides(colour = guide_legend(reverse=T))+
  
  geom_linerange(aes(ymin=exp(ci.lb), 
                     ymax=exp(ci.ub))) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  
  theme(legend.position = c(0.8,0.86),
        axis.text = element_text(size = 9.5),
        legend.text = element_text(size = 9.5),
        axis.title = element_text(size = 11),
        legend.background = element_blank()#,
        #plot.margin = unit(c(2,4,2,1), "lines")
  )

figureS9 <- cowplot::plot_grid(figureS9a, 
                               figureS9b, 
                               ncol = 2,  rel_widths  = c(1,1) ,labels = "AUTO", align = 'h', axis = 'b', greedy = F)



# Save to file (ggsave rather than setEPS() to preseve transparencies)
ggsave("./results/figureS9.eps", device=cairo_ps, width = 16, height = 10, units= 'in')
Sys.sleep(0.5)
figureS9
dev.off()