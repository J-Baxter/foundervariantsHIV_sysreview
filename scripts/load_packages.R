###################################################################################################
###################################################################################################
# Loads dependencies for Inferring the multiplicity of founder variants initiating HIV-1 
# infection: a systematic review and IPD meta-analysis

# Packages should be installed using renv::restore(), which will load the required packages and 
# versions from renv.lock
###################################################################################################
###################################################################################################
require(tidyverse)

pkgs <- c( "cowplot", "tidyverse", "emmeans", "lme4", "parallel", "performance", "scales",
          "ggsci", "grDevices", "aod", "metafor", "RColorBrewer", "ggeffects", "reshape2", 
          "data.table", "glmmTMB", "influence.ME", "meta", "mltools", "actuar")

for (pkg in pkgs){
  
  if (!(pkg %in% installed.packages())){
    
    install.packages(pkg)
    
    library(pkg, character.only = T)
    
  }else if (pkg %in% installed.packages()){
    
    library(pkg, character.only = T)
  }
}


###################################################################################################
# Create directories for results and figures

ddmonthyy <- format(Sys.Date(), '%d%b%y')
check_dirs <- paste(c('./results', './figures'), ddmonthyy, sep = '/')
dirs <- list.dirs()

for (check_dir in check_dirs){
  
  if (!(check_dir %in% dirs)){
    dir.create(check_dir)
  }

}

results_dir <- check_dirs[[1]]
figs_dir <- check_dirs[[2]]


###################################################################################################
# END #
###################################################################################################