###################################################################################################
###################################################################################################
# Installs dependencies for Inferring the multiplicity of founder variants initiating HIV-1 
# infection: a systematic review and IPD meta-analysis
###################################################################################################
###################################################################################################

required_packages <- c("tidyr", "lme4", "dplyr", "mltools", "data.table", "metafor", "aod", "ggplot2", 
                       "influence.ME", "reshape2", "ggsci", "forcats", "RColorBrewer", "cowplot", 
                       "stringr", "parallel", "performance",  "emmeans", "insight", "magick", "ggeffects")


for (package in required_packages){
  print(package)
  require(package, character.only = T)
}


installed_packages <- (.packages())

if(all(required_packages %in% installed_packages)){
  print('All packages installed.')
}else{
  missing_packages <- required_packages[!(required_packages %in% installed_packages)]
  error_message <- sapply(missing_packages, function(x) paste0(x,'\n'))
  stop('The following packages are not installed: \n', error_message)
}
     
###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################