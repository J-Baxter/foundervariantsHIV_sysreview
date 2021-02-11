###################################################################################################
###################################################################################################
# Visualisation for pooling and sensitivity analyses

###################################################################################################
# Dependencies
library(ggplot2)
library(ggsci)
library(kableExtra)

# Plot for pseudo-bootstrap replicates of participant selection
PltBoot <- function(data, intercept, ci.lb, ci.ub){
  plt <- ggplot(data) + 
    geom_histogram(aes(x = estimate),color="black", fill="grey96", binwidth = 0.0005) + 
    #geom_vline(aes(xintercept=mean(estimate)),color="#DC0000B2", linetype="dashed", size=1) +
    geom_vline(aes(xintercept= intercept),color="#00A087B2", linetype="dashed", size=1) +
    theme_classic() +
    scale_y_continuous(name = "Frequency" ,limits=c(0,50),expand = c(0, 0))+
    scale_x_continuous(name = "Probability of Multiple Founders", limits=c(0,0.4), expand = c(0.01, 0.01))+
    annotate( "rect", ymin=0, ymax=Inf, xmin=ci.lb, xmax=ci.ub ,alpha = .1, fill = "#00A087B2") +
    theme(axis.text = element_text(size = 10.5,  family = "sans"),
          legend.text = element_text(size = 10.5,  family = "sans"),
          axis.title = element_text(size = 13,  family = "sans"))
  
  return(plt)
}


###################################################################################################
###################################################################################################
# Import data
pooled_models <-  readxl::read_xlsx("pooled_models.xlsx")
og_models <- pooled_models[1:4,]
# resampling data

###################################################################################################
###################################################################################################
# Panel: Table of estimate and tau for pooling

tbl <- kbl(sensitivity_df, digits = 3) %>%
  kable_classic(html_font = "Arial") %>%
  add_header_above(c(" " = 1, 
                     'Summary Estimate' = 3,
                     'Exclusion of Small (<10) Studies' = 3,
                     'Exclusion of Studies reporting 0 MF' = 3))

###################################################################################################
# Plot log odds of individual studies. should be normally distiributed to satisfy binomial-normal model.
ggplot(twostep_binorm.step1, aes(x=log_or)) + geom_histogram(binwidth = 0.25,color="black", fill="white")+
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))


###################################################################################################
###################################################################################################
# Panel: Sensitivity analyses (exclusion criteria and resampling)


# Exclusion critera dot and whisker plot
senseplot <- ggplot(cbind.data.frame(pooled_models[,1:2],mapply(metafor::transf.ilogit, pooled_models[,3:5])),aes(x= model, y = estimate, color = analysis)) +
  geom_point( shape = 4, size = 4, position = position_dodge(0.5)) + 
  scale_y_continuous(name = "Probability of Multiple Founders",limits=c(0,.5),expand = c(0.01, 0.01)) +
  scale_x_discrete(name = "Model", labels = c(
    twostep.bn = "BN",
    twostep.bb = "BB",
    onestep.strat = "Stratified",
    onestep.rand = "Random "))+
  theme_classic() + 
  coord_flip()+
  geom_linerange(aes(ymin=estimate.lb, ymax= estimate.ub, color = analysis), position = position_dodge(0.5)) +
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
                        intercept = mapply(transf.ilogit, og_models[,3]), 
                        ci.lb = mapply(transf.ilogit, og_models[,4]),
                        ci.ub = mapply(transf.ilogit, og_models[,5]),
                        SIMPLIFY = FALSE)

plt_boot.grid <- cowplot::plot_grid(plotlist = plt_boot.list , align = "hv" , nrow = 2, ncol = 2 , labels = c("B" , 'C' , "D" , "E"))

cowplot::plot_grid(senseplot, 
                   plt_boot.grid , nrow = 2,  rel_heights = c(1, 1) ,labels = "A")

###################################################################################################
# Panel: Sensitivity analyses (Faceted Influence Plots)


