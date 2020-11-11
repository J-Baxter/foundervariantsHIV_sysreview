#expand pooled method data and group according to single/multiple founder thresholds as stipulated in source paper
#Keele et al. PNAS 2008
#102 participants included
#methods include model, distance, phylogenetic and haplotype(ish) categories
#requires csv with column names in the format as follows: participant.var, method1.results, method2.results

#dependencies
library(dpylr)
library(stats)

#define functions required in script
#split column names about '.'
getsplits <- function(cols) {
  splitted <- strsplit(cols, '[.]')
  first_element <- unlist(splitted)[1]
  
  return(first_element)
  }


#split methods into independent dataframes
groupbycols <- function(df){
  stopifnot(is(df , 'data.frame'))
  col_names = names(df)
  split <- lapply(col_names, getsplits) %>% unique()
  selected <- lapply(split , function(x) df[ , grepl(x, col_names) ])
  
  names(selected) <- split
  return(selected)
}


#create labelled dataframes with individual methods
labeldfs <- function(listofdfs){
  covar = listofdfs[[1]]
  measures = listofdfs[-1]
  varnames = names(measures)
  
  withparticipant <- lapply(measures, function(x) cbind.data.frame(covar, x))
  
  output <- list()
  for (i in 1:length(withparticipant)){
    df <- cbind.data.frame(withparticipant[[i]] , varnames[i])
    names(df)[ncol(df)] <- 'method'
    output[[i]] <- df
  }
  
  names(output) <- varnames
  return(output)
}


#classification according to generalised threshold
classifyfixed <- function(df, criteria){
  require(dplyr)
  require(stats)
  
  multiple <- do.call(dplyr::filter_, list(df, criteria[[1]])) %>% cbind(founder.multiplicity = 'multiple') 
  single <- do.call(dplyr::filter_, list(df, criteria[[2]]))  %>% cbind(founder.multiplicity = 'single') 
  nomeasure <- do.call(dplyr::filter_, list(df, criteria[[3]])) 
  
  if (nrow(nomeasure) > 0){
     nomeasure_na <- cbind.data.frame(nomeasure, founder.multiplicity = 'NA')
    classified <- rbind.data.frame(single, multiple, nomeasure_na)
  }else{
    classified <- rbind.data.frame(single, multiple)
    }

  stopifnot(nrow(df)==nrow(classified))
  return(classified)
}


#keele et al classification for BEAST lower bound tmrca (relative to particpant.fiebig)
classifyrelative <- function(df){
  EDI <- data.frame(fiebig.stage = c('I', 'II', 'III', 'IV', 'V'), EDI.upperbound = c(8, 34, 37, 43, 154)) 
  
  coord <- which(EDI$fiebig.stage %in% df$participant.feibig)
  
  criteria_b <- c(~ beast.lower > EDI[coord,2], ~ beast.lower <= EDI[coord,2], ~ is.na(beast.lower)) #is "VI" force NA

  classified <- classifyfixed(fiebig_split$III, criteria_b)
  
  return(classified)
}


#assign multiple/single founder classification according to a constant threshold value
#currently written for poisson and distance only
#must edit threshold, df name and criteria to adapt to other vars (not generalisable at present)

assignclassification<- function(listofdfs, threshold){
  stopifnot(length(threshold) == length(listofdfs))

  #Distance
  distance_df <-  listofdfs$distance
  DISTANCE <- threshold[[1]]
  criteria_d <- c(~ distance.meanpercent >= DISTANCE, ~ distance.meanpercent < DISTANCE, ~ is.na(distance.meanpercent))
  distance_classified <- classifyfixed(distance_df, criteria_d)
  stopifnot(nrow(distance_df) == nrow(distance_classified))
  
  #Poisson
  poisson_df <- listofdfs$poisson
  POISSON <- threshold[[2]]
  criteria_p <- c(~ poisson.GOF >= POISSON, ~ poisson.GOF < POISSON, ~ is.na(poisson.GOF))
  poisson_classified <- classifyfixed(poisson_df, criteria_p)
  stopifnot(nrow(poisson_df) == nrow(poisson_classified))
  
  #BEAST TMRCA
  beast_df <- listofdfs$beast
  fiebig_split <- split.data.frame(beast_df , beast_df$participant.feibig)
  beast_classified <- lapply(fiebig_split , classifyrelative) %>% do.call(what =rbind.data.frame)
  stopifnot(nrow(beast_df) == nrow(beast_classified))
  
  #out
  stopifnot(nrow(poisson_classified) == nrow(distance_classified) == nrow(beast_classified))
  output <- list(distance_classified, poisson_classified, beast_classified)
  names(output) <- names(listofdfs)
  return(output)
}


#wrapper function for script.
stratifypooledmethods <- function(data, thresholds, doi){

  labelled_dfs <- groupbycols(keele_combined) %>%
    labeldfs()
  
  #classification
  classified_dfs <- assignclassification(labelled_dfs, THRESHOLDS, EDI) #note exclusion
  
  #write output csv(s) to file
  yymmdd <- format(Sys.Date(), '%Y-%b-%d')
  filenames <- lapply(names(classified_dfs) , function(x) paste0('keele', '_', x , '_', yymmdd, '.csv'))
  names(classified_dfs) <- filenames
  
  for (filename in filenames){
    write.csv(filename, file = filename)
  }
}


###############################
##START##
#import data and set groups
keele_combined <- read_csv("keele_combined.csv")

#define thresholds as stipulated in keele et al 2008
THRESHOLDS <- c(0.86,0.05)
names(THRESHOLDS) <- c('DISTANCE','POISSON') 

#stage VI is open ended so cannot place an upper bound of time to infection with any confidence
#estimates obtained from Lee et al Theoretical Bio 2009

#run script
stratifypooledmethods(keele_combined, THRESHOLDS)

#END#
###############################