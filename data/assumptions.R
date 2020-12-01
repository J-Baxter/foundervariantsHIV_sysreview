#evaluate distribution of data to assumptions of meta analysis
#import dataset
data_master<- read.csv("data_master.csv", na.strings = "NA")
data_master.nona <- df[!is.na(df$multiple.founders),]
data_master.labelled <- unite(data_master.nona, "publication", c(author ,year), sep = '_')

data_master.labelled$multiple.founders = 1 - (as.numeric(data_master.labelled$multiple.founders)-1)

df_grouped <- data_master.labelled %>% 
  group_by(publication) %>%
  summarise(subjects = n(), multiplefounders = sum(multiple.founders))
contcorrection <- function(x){if(x==0){x=0.3}else if(x==1){x=0.3} else{x=x}}
props <- df_grouped$multiplefounders/df_grouped$subjects
props.corrected.log =lapply(props , contcorrection) %>% unlist() %>% log()
fit.norm.cont = fitdist(props.corrected.log , distr = 'norm')
plot(fit.norm.cont)
gofstat(fit.norm.cont)