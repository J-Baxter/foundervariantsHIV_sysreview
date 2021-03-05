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
library(magick)

###################################################################################################
###################################################################################################
# Import data
setwd("./data")

###################################################################################################
# Plot and tabulate random effects structure selection

raneff_selection <- read.csv('raneff_selection.csv') 
raneff_ic <- raneff_selection[,c(1,2,3)] %>%
  reshape2::melt() %>%
  `colnames<-`(c('X' , 'criteria' , 'value'))

raneff_plot <- ggplot(raneff_selection) + 
  geom_point(aes(x = X, y = estimate))+
  geom_linerange(aes(x = X, ymin=estimate.lb, 
                     ymax= estimate.ub))+
  geom_line(aes(x = X,
                y = value/2500,
                color = criteria, 
                group = criteria), data = raneff_ic) +
  geom_point(aes(x = X,
                 y = value/2500,
                 color = criteria,
                 group = criteria), data = raneff_ic)+
  scale_y_continuous(name = 'Probability of Multiple Founders',
                     expand = c(0,0.02),
                     limits = c(0,1), 
                     sec.axis = sec_axis(~.*2500 , name = 'AIC/BIC'))+
  scale_x_discrete(name = 'Random Effects Structure',
                   labels = stringr::str_wrap(as.factor(raneff_selection$X), width = 30))+
  theme_classic() + 
  scale_color_npg()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Format dataframe for table
raneff_comb <- cbind.data.frame("Estimate" = paste0(round(raneff_selection$estimate, digits = 3), ' ' ,
                                                    '[' , round(raneff_selection$estimate.lb, digits = 3) ,' - ',
                                                    round(raneff_selection$estimate.ub, digits = 3), ']'),
                                row.names = raneff_selection$X )

raneff_tbl_df <- cbind.data.frame(raneff_comb, raneff_selection[,c(2,3)])

knitr::kable(raneff_tbl_df, 
                   format = 'latex',
                   booktabs = T,
                   col.names = c('Estimate', 'AIC', 'BIC'),
                   escape = FALSE,
                   align = 'c', 
                   linesep = c("\\addlinespace")) %>% 
  kable_styling(latex_options = "striped") %>%as_image(filename = 'something.png')
###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################