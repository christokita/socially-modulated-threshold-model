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
    install.packages(new_packages, repo = 'https://cloud.r-project.org/')
  } 
}

####################
# Tag runs for parallel computing
####################
label_parallel_runs <- function(matrix, n, simulation, chunk) {
  rows <- nrow(matrix)
  col_names <- colnames(matrix)
  matrix <- cbind(matrix, rep(n, rows))
  matrix <- cbind(matrix, rep(simulation, rows))
  matrix <- cbind(matrix, rep(chunk, rows))
  colnames(matrix) <- c(col_names, 'n', 'sim', 'chunk')
  return(matrix)
}

####################
# Save parallel computing data
####################
save_parallel_data <- function(data, path, sub_directory, n, chunk) {
  binded_data <- do.call('rbind', data)
  write_path <- paste0(path, "/", sub_directory, "/", n, "-", chunk, ".Rdata")
  save(binded_data, file = write_path)
}