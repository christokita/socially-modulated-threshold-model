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

label_parallel_runs_beta <- function(matrix, beta, simulation, chunk) {
  rows <- nrow(matrix)
  col_names <- colnames(matrix)
  matrix <- cbind(matrix, rep(beta, rows))
  matrix <- cbind(matrix, rep(simulation, rows))
  matrix <- cbind(matrix, rep(chunk, rows))
  colnames(matrix) <- c(col_names, 'beta', 'sim', 'chunk')
  return(matrix)
}

####################
# Save parallel computing data
####################
save_parallel_data <- function(data, path, sub_directory, n, chunk) {
  n <- str_pad(string = n, width = 3, pad = "0")
  chunk <- str_pad(string = chunk, width = 2, pad = "0")
  write_path <- paste0(path, "/", sub_directory, "/", n, "-", chunk, ".Rdata")
  save(data, file = write_path)
}

save_parallel_data_beta <- function(data, path, sub_directory, beta, chunk) {
  chunk <- str_pad(string = chunk, width = 2, pad = "0")
  write_path <- paste0(path, "/", sub_directory, "/", beta, "-", chunk, ".Rdata")
  save(data, file = write_path)
}
