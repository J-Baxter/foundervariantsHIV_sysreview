#TEMP file for code for plots to compare senstivity analyses with original models for pooling

library(ggplot2)
library(ggsci)

pooled_models <-  read_xlsx("pooled_models.xlsx")
  
ggplot(cbind.data.frame(pooled_models[,1:2],mapply(transf.ilogit, pooled_models[,3:5])),aes(x= model, y = estimate, color = analysis)) +
  geom_point( shape = 4, size = 4, position = position_dodge(0.5)) + 
  scale_y_continuous(name = "Pooled Proportion of Founder Variant Multiplicity",limits=c(0,.5),expand = c(0.01, 0.01)) +
  scale_x_discrete(name = "Model", labels = c(
    twostep.bn = "Binomial-Normal",
    twostep.bb = "Beta-Binomial",
    onestep.strat = "GLMM with Stratified Intercept",
    onestep.rand = "GLMM with Random Intercept"))+
  theme_classic() + 
  coord_flip()+
  geom_linerange(aes(ymin=estimate.lb, ymax= estimate.ub, color = analysis), position = position_dodge(0.5)) +
  scale_color_npg(name = NULL, labels = c(
    no_small = "Studies with (n<10) omitted",
    no_zeros = "Studies with (p=0) omitted",
    original = "Full analysis")) + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10.5,  family = "sans"),
        legend.text = element_text(size = 10.5,  family = "sans"),
        axis.title = element_text(size = 13,  family = "sans")
  )