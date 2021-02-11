###################################################################################################
###################################################################################################
# Visualisation for pooling and sensitivity analyses

###################################################################################################
# Dependencies
library(ggplot2)
library(ggsci)
library(kableExtra)
library(metafor)
library(readxl)

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
                       limits=c(0,50),
                       expand = c(0, 0)) +
    
    scale_x_continuous(name = "Probability of Multiple Founders",
                       limits=c(0,0.4), 
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
pooled_models <-  readxl::read_xlsx("pooled_models.xlsx")
models <- c('Two-Step Binomial Normal',
            'One-Step Binomial (random slope) and correlated intercept',
            'One-Step Binomial (uncorrelated random intercept and slope)',
            "Two-Step Beta-Binomial")
og_models <- cbind(pooled_models[1:4,1], mapply(transf.ilogit, pooled_models[1:4,3:5]) %>% round(digits = 3))

# resampling data

###################################################################################################
###################################################################################################
# Panel: Table of estimate and tau for pooling

og_models_formatted <- cbind.data.frame("probability" = paste0(og_models$estimate, ' ' ,
                                                               '[' , og_models$estimate.lb ,' - ',
                                                               og_models$estimate.ub, ']'),
                                        row.names =models )

tbl <- kbl(og_models_formatted, longtable = T, booktabs = T, col.names = c("Probability of Multiple Founders")) %>%
  kable_classic(html_font = "Arial")


###################################################################################################
# Plot log odds of individual studies. should be normally distiributed to satisfy binomial-normal model.
ggplot(twostep_binorm.step1, aes(x=log_or)) + geom_histogram(binwidth = 0.25,color="black", fill="white")+
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))


###################################################################################################
###################################################################################################
# Panel: Sensitivity analyses (exclusion criteria and resampling)


# Exclusion critera dot and whisker plot
senseplot <- ggplot(cbind.data.frame(pooled_models[,1:2],mapply(metafor::transf.ilogit, pooled_models[,3:5])),
                    aes(x= model, y = estimate, color = analysis)) +
  
  geom_point( shape = 4, 
              size = 4,
              position = position_dodge(0.5)) + 
  
  scale_y_continuous(name = "Probability of Multiple Founders",
                     limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  
  scale_x_discrete(name = "Model", 
                   labels = c(
                     twostep.bn = "BN",
                     twostep.bb = "BB",
                     onestep.strat = "Stratified",
                     onestep.rand = "Random")) +
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

# Resampling histograms
plt_boot.list <- mapply(PltBoot, data = t, 
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
  
  scale_x_discrete(labels = influence.labs) +
  
  theme_classic() + 
  
  coord_flip() +
  
  facet_grid(cols = vars(model),
             labeller = labeller(model = c(
               twostep_binorm = "Binomial-Normal",
               twostep_betabi = "Beta-Binomial",
               onestep_strat = "GLMM Stratified",
               onestep_rand = "GLMM Random"))) +
  
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

tiff(filename = 'influenceplot.tif', res = 180, width=2750, height=2650 , units = 'cm', pointsize = 10)

plt

dev.off()


