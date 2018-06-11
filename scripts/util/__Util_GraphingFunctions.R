##################################################
#
# Graphing functions
#
##################################################


####################
# Graph probabilities
####################
graphTaskDistr <- function(TaskMat) {
  # Check for number of tasks
  cols <- ncol(TaskMat)
  if (cols == 5) {
    # Plot task distributions for 2 tasks
    plot_TaskMat <- as.data.frame(TaskMat)
    gg_dist <- ggplot(data = plot_TaskMat, aes(y = Task1, x = set)) +
      geom_point(aes(colour = n)) +
      theme_classic() +
      ylab("Frequency Task 1") 
      
      
  } else if (cols == 6) {
    # Plot task distributions for 3 tasks
    plot_TaskMat <- as.data.frame(TaskMat)
    gg_simplex <- ggtern(data = plot_TaskMat, aes(Task1, Task2, Task3)) +
      stat_density_tern(
        geom = 'polygon',
        aes(fill = ..level..),
        alpha = 0.1,
        bins = 100,
        show.legend = FALSE) +
      # geom_density_tern(aes(color = ..level..), 
      #                   show.legend = FALSE, 
      #                   alpha = 0.7) +
      geom_point(color = "black",
                 fill = "black", 
                 stroke = 0,
                 size = 1.5, 
                 alpha = 1) +
      theme_bw() +
      # scale_colour_distiller(palette = "RdYlBu") +
      scale_fill_distiller(palette = "RdYlBu") +
      theme(axis.title = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 12, face = "italic", hjust = 0.5)) +
      labs(x = "Task 1",
           y = "Task 2",
           z = "Task 3", 
           title = paste0("Simulation ", sim))
    return(gg_simplex)
  } else {
    message("ERROR: Must be 2 or 3 tasks for graphing to work")
  }
}

####################
# Turn matrix data and turn into igraph object
####################
matrixToGraphObject <- function(AdjacencyMat, TaskStateMat) {
  # Turn into graph from adjacency matrix
  graphObj <- graph_from_adjacency_matrix(adjmatrix = AdjacencyMat, mode = c("undirected"))
  # Get task state of each individual
  tasks <- lapply(1:nrow(X_g), function(i) {
    task <- which(X_g[i, ] == 1)
    label <- paste("Task", task)
  })
  tasks <- unlist(tasks)
  # Add attribute and return
  graphObj <- set_vertex_attr(graphObj, "Task", value = tasks)
  return(graphObj)
}

####################
# Empty theme for plotting graphs
####################
theme_invisible <- function() {
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank())
}

####################
# Multiplot function
####################
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

####################
# Graph network by tasks
####################
graphNetwork <- function(graph) {
  # Get task atrributes
  tasks <- vertex_attr(graph)[[2]]
  # Make graph
  g <- ggraph(graph = graph, layout = "fr") +
    geom_edge_link(alpha = 0.8) +
    geom_node_point(size = 2, aes(colour = tasks)) +
    theme_invisible() +
    scale_colour_brewer(palette = "Set1") +
    labs(colour = "Task State") +
    theme(legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8))
  return(g)
}



