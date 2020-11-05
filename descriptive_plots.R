#descriptive plots for selected studies
#J Baxter
#14 July 2020

library(ggplot2)
library(ggsci)
library(dplyr)
library(forcats)
library(tidyr
        )
#set wd
setwd("./data/")

#import data
df <- read.csv("data_master.csv", na.strings = "NA")

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

#remove NAs from founder.multiplicity column

#method
p1 <- ggplot(df, aes(grouped.method))+
  geom_bar()+
  scale_color_npg()+
  scale_fill_npg()+
  theme_classic()+
  xlab('Grouped Method')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Count')

p1

#Risk
names <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure, names)


p2 <- ggplot(exposures_df, aes(x = reported.exposure , y = frequency))+
  geom_bar(stat = 'identity' , aes(fill = sub.exposure, colour = sub.exposure), position = 'stack')+
  scale_color_npg()+
  scale_fill_npg()+
  theme_classic()+
  xlab('Exposure')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Count')

p2

#Seroconversion
p3 <- ggplot(df, aes(participant.seropositivity))+ 
  geom_bar()+
  scale_color_npg()+
  scale_fill_npg()+
  theme_classic()+
  xlab('Seropositivity')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Count')

p3

#multiplicity
p4 <- ggplot(df, aes(grouped.subtype))+
  geom_bar()+
  scale_color_npg()+
  scale_fill_npg()+
  theme_classic()+
  xlab('Subtype')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Count')

p4
library(ggpubr)
ggarrange(p1,p2,p3,p4,
          ncol = 2 , nrow = 2 , labels = "AUTO")

#poisson plot
x<- 0:15
y <-  dpois(x, lambda=0.90)
df <- cbind.data.frame(x,y )
p5 <- ggplot(df , aes(x,y))+
  geom_line(colour = 'red' , size = 2)+
  geom_bar(stat = 'identity', colour = 'black',fill = NA)+
  theme_classic()+
  xlab('Hamming Distance')+
  ylab('')

x <- 0:15
y <-  runif(16,max = 0.4, min = 0)
fail.df <- cbind.data.frame(x,y)
p6 <- ggplot(fail.df , aes(x,y))+
  geom_bar(stat = 'identity', colour = 'black',fill = NA)+
  geom_line(data = df, color = "red", size = 2)+
  theme_classic()+
  xlab('Hamming Distance')+
  ylab('')

ggarrange(p5,p6,
          ncol = 1 , nrow = 2 , labels = "AUTO")




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
  scale_color_npg()+
  scale_fill_npg()+
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




