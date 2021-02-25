###################################################################################################
###################################################################################################
# Visualisation for pooling and sensitivity analyses
# Requires some model inputs from binary_pooling.R
###################################################################################################
# Dependencies
library(ggplot2)
library(ggsci)
library(kableExtra)
library(metafor)
library(dplyr)


# Plot for pseudo-bootstrap replicates of participant selection
PltBoot <- function(data, intercept, ci.lb, ci.ub){
  plt <- ggplot(data) + 
    geom_histogram(aes(x = estimate),
                   color="black", 
                   fill="grey96",
                   binwidth = 0.0005) + 
    
    #geom_vline(aes(xintercept=mean(estimate)),color="#DC0000B2", linetype="dashed", size=1) +
    geom_vline(aes(xintercept= intercept),
               color="#00A087B2", 
               linetype="dashed", 
               size=1) +
    
    theme_classic() +
    
    scale_y_continuous(name = "Frequency" ,
                       limits=c(0,75),
                       expand = c(0, 0)) +
    
    scale_x_continuous(name = "Probability of Multiple Founders",
                       limits=c(0,0.35), 
                       expand = c(0.01, 0.01)) +
    
    annotate( "rect",
              ymin=0, ymax=Inf,
              xmin=ci.lb, xmax=ci.ub
              ,alpha = .1,
              fill = "#00A087B2") +
    
    theme(axis.text = element_text(size = 10.5,  family = "sans"),
          legend.text = element_text(size = 10.5,  family = "sans"),
          axis.title = element_text(size = 13,  family = "sans"))
  
  return(plt)
}


###################################################################################################
###################################################################################################
# Import data
setwd("./data")

influence_df <- read.csv("bp_sa1.csv") %>% arrange(., model)

pooled_models <-  read.csv('bp_estsa2sa3.csv')

resampled_models <- read.csv('bp_resampl.csv')

models <- c('Two-Step Binomial Normal',
            'One-Step Binomial',
            "Two-Step Beta-Binomial")

og_models <- cbind("model" = pooled_models[1:3,1], pooled_models[1:3,3:9] %>% round(digits = 3)) %>% arrange(., model)

###################################################################################################
###################################################################################################
# Panel: Table of estimate and tau for pooling

og_models_formatted <- cbind.data.frame("Estimate" = paste0(og_models$estimate, ' ' ,
                                                               '[' , og_models$estimate.lb ,' - ',
                                                               og_models$estimate.ub, ']'),
                                        og_models[,5:8],
                                        row.names =models ) 
og_models_formatted[is.na(og_models_formatted )] <- "-"

# Not hapy with formatting
tbl <- kbl(og_models_formatted , 
          booktabs = T,
           col.names = c('Estimate', '$$\\hat{\\tau}^2$$', "$$\\text{Q}$$", "$$\\text{I}^2$$", "$$\\phi$$"),
           escape = FALSE,
           align = 'c', 
           ) %>%
  add_header_above(c(" " = 1, "Probability of Multiple Founders" = 1, "Heterogeneity" = 4), bold = T, line = T) %>%
 kable_classic(full_width = F, html_font = 'arial')

  
###################################################################################################
# Plot log odds of individual studies. should be normally distiributed to satisfy binomial-normal model.
ggplot(twostep_binorm.step1, aes(x=log_or)) + geom_histogram(binwidth = 0.25,color="black", fill="white")+
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))


###################################################################################################
###################################################################################################
# Panel: Sensitivity analyses (exclusion criteria and resampling)


# Exclusion critera dot and whisker plot
senseplot <- ggplot(pooled_models,
                    aes(x= model, y = estimate, color = analysis)) +
  
  geom_point( shape = 4, 
              size = 4,
              position = position_dodge(0.5)) + 
  
  scale_y_continuous(name = "Probability of Multiple Founders",
                     limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  
  scale_x_discrete(name = "Model", 
                   labels = c(
                     twostep_binorm = "BN",
                     twostep_betabi = "BB",
                     onestep_bi_rand= "Random")) +
  theme_classic() + 
  
  coord_flip() +
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_color_npg(name = NULL, labels = c(
    no_small = "Studies with (n<10) omitted",
    no_zeros = "Studies with (p=0) omitted",
    original = "Full analysis")) + 
  
  theme(legend.position = "top",
        axis.text = element_text(size = 10.5,  family = "sans"),
        legend.text = element_text(size = 10.5,  family = "sans"),
        axis.title = element_text(size = 13,  family = "sans"),
        plot.margin = unit(c(2,4,2,1), "lines")
  )

# Resampling histograms - requires model input from binary_pooling.R
resampled_models.list <- split.data.frame(resampled_models, resampled_models$analysis) %>%

plt_boot.list <- mapply(PltBoot, data = resampled_models.list, 
                        intercept = og_models[,2], 
                        ci.lb = og_models[,3],
                        ci.ub = og_models[,4],
                        SIMPLIFY = FALSE)

plt_boot.grid <- cowplot::plot_grid(plotlist = plt_boot.list , align = "hv" , nrow = 2, ncol = 2 , labels = c("B" , 'C' , "D" , "E"))

cowplot::plot_grid(senseplot, 
                   plt_boot.grid , nrow = 2,  rel_heights = c(1, 1) ,labels = "A")


###################################################################################################
# Panel: Sensitivity analyses (Faceted Influence Plots)
influence_df$trial <- gsub("_" , " ", influence_df$trial)

plt <- ggplot(influence_df,aes(x = trial , y = estimate) ) +
  geom_point() + 
  
  scale_y_continuous(limits=c(0,0.4),
                     expand = c(0, 0),
                     name = "Probability of Multiple Founders") +
  
  scale_x_discrete(labels = influence_df$trial) +
  
  theme_classic() + 
  
  coord_flip() +
  
  facet_grid(cols = vars(model),
             labeller = labeller(model = c(
               twostep_binorm = "Binomial-Normal",
               twostep_betabi = "Beta-Binomial",
               onestep_bi_rand = "GLMM Random"))) +
  
  geom_errorbar(data =influence_df, 
                aes(x = trial,
                    ymin=ci.lb,
                    ymax=ci.ub)) +
  
  geom_hline(data = og_models, 
             aes(yintercept = estimate),
             linetype="dashed", 
             color = "#DC0000FF",
             size=0.75) +
  
  geom_rect(data = og_models, 
            inherit.aes = FALSE,
            aes(ymin=estimate.lb, ymax=estimate.ub, 
                xmin = -Inf, xmax = Inf),
            fill = "#DC000033") +
  
  theme(axis.line.y =element_blank(),
    axis.title.y =element_blank(),
    axis.ticks.y=element_blank(),
    strip.text.x = element_text(face = "bold" , colour = 'black' , size = 10.5),
    panel.spacing.x = unit(0.6 , 'cm'),
    strip.background = element_rect(fill = NA, colour= NA)
  )

jpeg(filename = 'influenceplot.jpeg', res = 350, width=2750, height=4000 , units = 'px', pointsize = 10)

plt

dev.off()


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################