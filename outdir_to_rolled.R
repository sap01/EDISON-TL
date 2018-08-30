## Goal: Given an unrolled net, roll it into a rolled net.
##
#' @title Convert a given unrolled Dynamic Bayesian Network (DBN) into a rolled DBN using different rolling methods.
#' @param unrolled.net.adj.matrix.list Given time-varying network adjacency list. Its length = number of time-interval-specific nets. 
#' @param roll.method Which rolling method to use from {'any', 'all', or some real number in (0, 1), like - 0.5}.
#' @param allow.self.loop Boolean to decide whehter to allow self loop or not in the rolled DBN.
#' @return rolled.net.adj.matrix Return the rolled DBN adjacency matrix.
##
OutdirToRolled <- function(out.dir.name, edge.post.prob.threshold, 
                           roll.method, allow.self.loop, series.combn.method) {
  
  ##------------------------------------------------------------
  ## Begin: Load the Required Packages
  ##------------------------------------------------------------
  ##
  ##------------------------------------------------------------
  ## End: Load the Required Packages
  ##------------------------------------------------------------
  
  ##------------------------------------------------------------
  ## Begin: Load the Required External Functions
  ##------------------------------------------------------------
  init.path <- getwd()
  
  source(paste(init.path, 'unrolled_to_rolled.R', sep = '/'))
  ##------------------------------------------------------------
  ## End: Load the Required External Functions
  ##------------------------------------------------------------
  
  
  ## Absolute path of the output dir
  out.dir.abs.path <- paste(init.path, 'asset', out.dir.name, sep = '/')
  
  ## Names of the files inside the output dir
  out.dir.file.names <- list.files(out.dir.abs.path)
  
  ## Names of the files that contain unrolled nets
  unrolled.net.file.names <- grep('unrolled.net.adj.matrices.[0-9]+.RData', out.dir.file.names, value = TRUE)
  rm(out.dir.file.names)

  ## Initialize rolled net adjacency matrix
  rolled.net.adj.matrix <- NULL
    
  for (unrolled.net.idx in 1:length(unrolled.net.file.names)) {
    
    if (unrolled.net.idx == 1) {
      
      ## Loads obj 'unrolled.net.adj.matrices' 
      load(paste(out.dir.abs.path, unrolled.net.file.names[1], sep = '/'))
      
      ## source(paste(init.path, 'unrolled_to_rolled.R', sep = '/'))
      rolled.net.adj.matrix <- UnrolledToRolled(unrolled.net.adj.matrices, 
                                                edge.post.prob.threshold,  
                                                roll.method, 
                                                allow.self.loop)
      rm(unrolled.net.adj.matrices)
    } else {
      
      ## Loads obj 'unrolled.net.adj.matrices' 
      load(paste(out.dir.abs.path, unrolled.net.file.names[unrolled.net.idx], sep = '/'))
      
      ## source(paste(init.path, 'unrolled_to_rolled.R', sep = '/'))
      rolled.net.adj.matrix.to.combine <- UnrolledToRolled(unrolled.net.adj.matrices, 
                                                           edge.post.prob.threshold, 
                                                           roll.method, 
                                                           allow.self.loop)
      
      rm(unrolled.net.adj.matrices)
      
      rolled.net.adj.matrix <- rolled.net.adj.matrix + rolled.net.adj.matrix.to.combine
    }
  }
  rm(unrolled.net.idx)
  
  ##------------------------------------------------------------
  ## Begin: Evaluate the time-series combination threshold
  ##------------------------------------------------------------
  series.combn.threshold <- NULL
  
  if(is.character(series.combn.method)) {
    if (series.combn.method == 'any') {
      ## Insert an edge in the final rolled net iff 
      ## it is present in at least one interim rolled net
      series.combn.threshold <- 1
    } else if (series.combn.method == 'all') {
      ## Insert an edge in the final rolled net iff 
      ## it is present in all interim rolled nets
      series.combn.threshold <- length(unrolled.net.file.names)
    }
  } else if (is.numeric(series.combn.method)) {
    ## Insert an edge in the final lrolled net iff it is present in at least 
    ## (series.combn.method * length(unrolled.net.file.names)) number of 
    ## interim rolled nets
    
    if ((series.combn.method > 0) & (series.combn.method < 1)) {
      series.combn.threshold <- series.combn.method * length(unrolled.net.file.names)
    } else {
      # print('\'series.combn.method\' accepts numeric values in the interval (0,1)')
      stop('\'series.combn.method\' accepts numeric values in the interval (0,1)')
    }
  }
  ##------------------------------------------------------------
  ## End: Evaluate the time-series combination threshold
  ##------------------------------------------------------------

  ## Apply the time-series combination threshold on the rolled net adjacency matrix
  for (tgt.node.idx in 1:ncol(rolled.net.adj.matrix)) {
    for (src.node.idx in 1:nrow(rolled.net.adj.matrix)) {
      if (rolled.net.adj.matrix[src.node.idx, tgt.node.idx] >= series.combn.threshold) {
        rolled.net.adj.matrix[src.node.idx, tgt.node.idx] <- 1
      } else {
        rolled.net.adj.matrix[src.node.idx, tgt.node.idx] <- 0
      }
    }
    rm(src.node.idx)
  }
  rm(tgt.node.idx)
  
  return(rolled.net.adj.matrix)
}