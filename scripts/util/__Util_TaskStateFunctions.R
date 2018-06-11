##################################################
#
# Task performance/state functions
#
##################################################


####################
# Choose task probabilistically 
####################
update_task_performance <- function(task_probs, state_matrix, quit_prob) {
  # Create possible task space
  tasks <- seq(1:ncol(task_probs))
  # Loop through individuals
  for(row in 1:nrow(state_matrix)) {
    # Inactive workers randomly sample one stimulus
    if (sum(state_matrix[row, ]) == 0) {
      # Sample task probability
      tasks_order <- sample(x = tasks, size = length(tasks), replace = FALSE)
      # Loop through tasks and go with first one that results in activity
      for (task in tasks_order) {
        prob <- task_probs[row, task]
        activity <- sample(x = c(0, 1), size = 1, prob = c(1 - prob, prob))
        if (activity == 1) {
          state_matrix[row, task] <- activity
          break
        }
      }
    } 
    else { #active workers quit with certain probability
      quit_now <- sample(x = c("yes", "no"), size = 1, prob = c(quit_prob, (1 - quit_prob)))
      if (quit_now == "yes") {
        state_matrix[row, ] <- 0
      }
    }
  }
  # Return
  colnames(state_matrix) <- paste0("Task", 1:ncol(task_probs))
  rownames(state_matrix) <- paste0("v-", 1:nrow(task_probs))
  return(state_matrix) 
}
