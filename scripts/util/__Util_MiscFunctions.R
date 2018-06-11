##################################################
#
# Misc. functions
#
##################################################

####################
# Install missing packages
####################
install_missing_packages <- function(list_of_packages) {
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[ , "Package"])]
  if(length(new_packages)) {
    install.package(new_packages)
  } 
}