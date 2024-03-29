###################################################################################################
###################################################################################################
# IPD meta analysis of HIV founder variant multiplicity: Comparing sequencing technologies

# Calculates summary proportion of infections initiated by multiple founder variants
# models implemented:
# 1. Two-step binomial-normal model (Random effects, inverse variance pooling, reml estimator of tau)
# 2. One-step binomial GLMM allowing for clustering by study. random effects between studies
#    (random intercepts). approx ML fit

# Then, a univaribale metaregression, segregating by sequencing technology
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
  filter(!is.na(sequencing.method_)) %>%
  filter(sequencing.method_ != 'unknown') %>%
  mutate(sequencing.method_ = gsub('_', '', sequencing.method_) %>% as.factor())%>%
  droplevels()

table(df$sequencing.method_)

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
form <- "multiple.founders_ ~ sequencing.method_  + (1 | publication_)"

unipooled_effectstruct <- GetName(form, effects = 'fixed')

# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_",
                    "sampling.delay_",'alignment.bin_', 'sequencing.method_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21", 'NFLG', 'sangerSGA')

df <- SetBaseline(df, baseline.covar, baseline.level)
df$alignment.length_ <- scale(df$alignment.length_)

df_props <- CalcProps(df)



# Run models
unipooled_model <- CalcRandMetaReg(df, form, opt = '1')

# Check model convergence and singularity
unipooled_check <- CheckModels(unipooled_model)%>% 
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

estimates$analysis <- 'seq_tech'

pooled <- rbind("model" = pooled_models[1:2,1:5], estimates) %>% arrange(., model)


# Export csv with pooling and univariable metaregression to file
originals <- cbind.data.frame(estimates, heterogeneity)

# Figure S9a - Pooled Original vs Pooled Vaccine Only 
# Figure S9b - Vaccine Subgroups comparison

# Set colour palettes 
mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)] #c("#E64B35FF", "#4DBBD5FF")

# Set Labels
labs <- c('Multiple','Single')

figureS7a <- ggplot(df, aes(x = sequencing.method_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)))+
  scale_fill_brewer( palette = 'YlGn', labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  xlab('Sequencing Technology')+
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") + 
  scale_x_discrete(labels = c("Sanger (with SGA)",
                              "Illumina MiSeq",
                              "Roche 454",
                              "PacBio HiFi",
                              "Sanger (no SGA)",
                              "Sanger (with SGA precursor)")%>%
                     str_wrap(width = 10))+
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
         legend.position = c(0.75,0.86))


figureS7b <- ggplot(pooled,
                    aes(x= forcats::fct_rev(model), y = estimate, color = analysis)) +
  
  geom_point( shape = 4, 
              size = 4,
              position = position_dodge(0.5)) + 
  
  scale_y_continuous(name = "Probability of Multiple Founders",
                     limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  #coord_cartesian(ylim = c(0,.5))+
  
  scale_x_discrete(name = "Model", 
                   labels = c(
                     onestep_bi_rand = "GLMM",
                     twostep_binorm = "B-N"
                   )) +
  
  coord_flip() +
  
  guides(colour = guide_legend(reverse=T))+
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_colour_brewer(palette = 'YlGn', name = 'Analysis', labels = c(
    original = "Full analysis",
    seq_tech = "Sequence methods only")) + 
  theme_bw(base_family = "LM Sans 10") + 
  
  theme(legend.position = c(0.8,0.8),
        legend.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
        axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
        strip.text  = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) 
    

seqtech_ref <- cbind.data.frame(level = 'sangersga', est = 0, ci.lb = NA, ci.ub = NA)
seqtech_subgroup <- rbind.data.frame(unipooled_models.coef$fe[, c(2,3,5,6)], seqtech_ref)

levs <- c("sangersga", "2G:illuminamiseq", "2G:iontorrent","2G:roche454", "3G:PacBiohifi","sanger" ,   "sangerprecSGA")
seqtech_subgroup$level <- factor(x = seqtech_subgroup$level, levels = levs)

figureS7c <- ggplot(seqtech_subgroup,
                    aes(x = level , y = exp(est))) +
  
  geom_point( shape = 18, 
              size = 4) + 
  
  scale_y_continuous(name = "Odds Ratio",
                     #limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  coord_cartesian(ylim = c(0,6))+
  
  scale_x_discrete(name = "Sequencing Technology", 
                  labels = c( "sangerprecSGA" = "Sanger (with SGA precursor)" ,
                             '2G:roche454'= "Roche 454",
                             '3G:PacBiohifi'= "PacBio HiFi" ,
                             'sanger'= "Sanger (no SGA)",
                             '2G:illuminamiseq'= "Illumina MiSeq",
                             'sangersga' = "Sanger (with SGA)" ) 
                   ) +
  
  coord_flip() +
  
  guides(colour = guide_legend(reverse=T))+
  
  geom_linerange(aes(ymin=exp(ci.lb), 
                     ymax=exp(ci.ub))) +
  geom_hline(yintercept = 1, linetype = 'dashed')  + 
  theme_bw(base_family = "LM Sans 10") + 
  
  theme(legend.position = c(0.8,0.86),
        legend.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
        axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
        strip.text  = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) 


figureS7 <- cowplot::plot_grid(figureS7a,
                               figureS7b, 
                               figureS7c, 
                               ncol = 1,  rel_heights = c(1,1,1) ,labels = "AUTO", align = 'hv', axis = 'b', greedy = F)



# Save to file (ggsave rather than setEPS() to preseve transparencies)
ggsave("figureS6.eps", device=cairo_ps,  height = 240, width = 150, units = 'mm')
Sys.sleep(0.5)
figureS7
dev.off()