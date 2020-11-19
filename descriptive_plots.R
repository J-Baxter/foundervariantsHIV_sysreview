#descriptive plots for selected studies
#J Baxter
#14 July 2020

library(ggplot2)
library(ggsci)
library(dplyr)
library(forcats)
library(tidyr)
library(RColorBrewer)
#set wd
setwd("./data/")

#import data
df <- read.csv("data_master.csv", na.strings = "NA")
dfnona <- df[!is.na(df$multiple.founders),]

#function for stacking categories and calculating summary frequencies
stacked_categories <- function(x, catnames){
  characters <- as.character(x)
  split_list <- strsplit(characters , '[:]')
  split_df <- do.call(rbind.data.frame, split_list)
  
  freq <- split_df %>%
    group_by(split_df[,1] , split_df[,2]) %>%
    summarise(frequency = n())
  
  colnames(freq) <- catnames
  return(freq)
}


#my colours
mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)]
nb.cols <- 12
mycols_method <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols)
#remove NAs from founder.multiplicity column


##basic plots for panel 1 (four bar plots to panel displaying frequencies of grouped method, seropositivity, subtype and number of sequences)
#method
p1 <- ggplot(dfnona, aes(grouped.method))+
  geom_bar(aes(fill = multiple.founders, colour = multiple.founders))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  theme_classic()+
  xlab('Grouped Method')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

p1

#Seroconversion #stack infant and NA together
p3 <- ggplot(dfnona, aes(participant.seropositivity))+ 
  geom_bar(aes(fill = multiple.founders, colour = multiple.founders))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  theme_classic()+
  xlab('Seropositivity')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

p3

#subtype
p4 <- ggplot(dfnona, aes(grouped.subtype))+
  geom_bar(aes(fill = multiple.founders, colour = multiple.founders))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  theme_classic()+
  xlab('Subtype')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

p4

#number of sequences
#must drop NA and NGS first 

numseqs_df <- dfnona[!dfnona$sequencing.number %in% c('NGS', NA),]
numseqs_df$sequencing.number <- as.numeric(numseqs_df$sequencing.number)

p5 <- ggplot(numseqs_df , aes(sequencing.number))+
  geom_histogram(binwidth=5, aes(fill = multiple.founders, colour = multiple.founders))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  theme_classic()+
  xlab('Number of Consensus Genomes Analysed')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') +
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")
  
  
library(ggpubr)
ggarrange(p1,p3,p4, p5,
          ncol = 2 , nrow = 2 , labels = "AUTO" , common.legend = TRUE , legend = 'bottom' , align = 'hv') 


#bar plot for exposure
#Risk
names <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure, names)


p2 <- ggplot(exposures_df, aes(x = reported.exposure , y = frequency))+
  geom_bar(stat = 'identity' , aes(fill = sub.exposure, colour = sub.exposure), position = 'stack')+
  scale_color_manual(values = mycols_method)+
  scale_fill_manual(values = mycols_method)+
  theme_classic()+
  xlab('Exposure')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') + labs(fill = "Reported Exposure", colour = "Reported Exposure")

p2


#Year of Publication

#Minimum numcitqtber of founders plot





###others
df <- read.csv("yop.csv")
dates <- as.Date(ISOdate(df$Year.of.Publication, 1, 1))
df.yrs <- as.data.frame(dates)
start <- as.Date('1991-6-30')
end <- as.Date('2020-6-30')
colnames(df.yrs) <- 'Year.of.Publication'
p1 <- ggplot(df.yrs, aes(Year.of.Publication ))+
  geom_bar(fill = '#000066')+
  geom_vline(xintercept = as.numeric(dates[58]),linetype=2, colour="red" )+
  scale_color_brewer(palette = 'RdBu')+
  scale_fill_brewer(palette = 'RdBu')+
  theme_classic()+
  ylab('')+
  xlab('Year of Publication')+scale_x_date(date_labels = "%Y", date_breaks = '1 year' , limits = c(start,end))+
  theme(axis.text.x=element_text(angle=45 , hjust = 1))

p1


df <- read.csv("subtype.csv")
dates <- as.Date(ISOdate(df$Year.of.Publication, 1, 1))
df.yrs <- as.data.frame(dates)
start <- as.Date('1991-6-30')
end <- as.Date('2020-6-30')
colnames(df.yrs) <- 'Year.of.Publication'
p1 <- ggplot(df.yrs, aes(Year.of.Publication ))+
  geom_bar(fill = '#000066')+
  geom_vline(xintercept = as.numeric(dates[58]),linetype=2, colour="red" )+
  scale_color_npg()+
  scale_fill_npg()+
  theme_classic()+
  ylab('')+
  xlab('Year of Publication')+scale_x_date(date_labels = "%Y", date_breaks = '1 year' , limits = c(start,end))+
  theme(axis.text.x=element_text(angle=45 , hjust = 1))

p1




