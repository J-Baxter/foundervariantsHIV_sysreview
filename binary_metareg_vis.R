###################################################################################################
###################################################################################################
# Visualisation for meta-regression
###################################################################################################
###################################################################################################
# Dependencies
library(ggplot2)
library(ggsci)
library(kableExtra)
library(metafor)
library(dplyr)

###################################################################################################
###################################################################################################
# Import data
setwd("./data")

influence_df <- read.csv("bp_sa1.csv") %>% arrange(., model)

###################################################################################################
# Plot and tabulate random effects structure selection

raneff_selection <- read.csv('raneff_selection.csv')
ranneff_ic <- raneff_selection[,c(1,2,3)] %>%
  reshape2::melt() %>%
  `colnames<-`(c('X' , 'criteria' , 'value'))

replot <- ggplot(raneff_selection) + 
  geom_point(aes(x = X, y = estimate))+
  geom_linerange(aes(x = X, ymin=estimate.lb, 
                     ymax= estimate.ub))+
  geom_line(aes(x = X, y = value/2500, color = criteria, group = criteria), data = ranneff_ic) +
  geom_point(aes(x = X, y = value/2500, color = criteria,group = criteria), data = ranneff_ic)+
  scale_y_continuous(name = 'Probability of Multiple Founders', expand = c(0,0.02), limits = c(0,1), sec.axis = sec_axis(~.*2500 , name = 'AIC/BIC'))+
  theme_classic() + 
  scale_color_npg()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))



###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################