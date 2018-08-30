## Goal: Given an EDISON output directory, produce an unrolled net 
##
###########################################################################################
## Goal: Given multiple sets of time-varying nets, aggregate them into a 
## single set of time-varying nets. Usually, each input set of time-varying nets
## is inferred from a distinct time series.
##
#' @title Convert a given unrolled Dynamic Bayesian Network (DBN) into a rolled DBN using different rolling methods.
#' @param unrolled.net.adj.matrix.list Given time-varying network adjacency list. Its length = number of time-interval-specific nets. 
#' @param allow.self.loop Boolean to decide whehter to allow self loop or not in the rolled DBN.
#' @param series.combn.method Which method to use for combining multiple time series from {'any', 'all', or some real number in (0, 1), like - 0.5}.
#' @return output.unrolled.net Return the adjacency matrices of the unrolled DBN.
###########################################################################################
OutdirToUnrolled <- function(out.dir.name, edge.post.prob.threshold, 
                             allow.self.loop, series.combn.method) {
  
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
  
  num.time.series <- length(unrolled.net.file.names)
  
  ## Initialize adjacency matrices of the unrolled net
  output.unrolled.net <- NULL
  
  for (unrolled.net.idx in 1:length(unrolled.net.file.names)) {
    
    ## Loads obj 'unrolled.net.adj.matrices' 
    load(paste(out.dir.abs.path, unrolled.net.file.names[1], sep = '/'))
    
    ## A list of length = num of time pts.
    ## Each element represnts the net adjacency matrix of the
    ## corresponding time pt.
    unrolled.net.adj.matrices <- unrolled.net.adj.matrices$probs.all
    
    for (time.pt.idx in 1:length(unrolled.net.adj.matrices)) {
      
      time.pt.spec.adj.matrix <- unrolled.net.adj.matrices[[time.pt.idx]]
      
      ## Remove edges with marginal posterior probability less than
      ## 'edge.post.prob.threshold'
      time.pt.spec.adj.matrix[time.pt.spec.adj.matrix < edge.post.prob.threshold] <- 0
      
      ## Convert 'time.pt.spec.adj.matrix' to a binary matrix of zeros and ones
      time.pt.spec.adj.matrix[time.pt.spec.adj.matrix != 0] <- 1
      
      ## Remove self loops if 'allow.self.loop' = FALSE
      if (!allow.self.loop) {
        diag(time.pt.spec.adj.matrix) <- 0
      }
      
      unrolled.net.adj.matrices[[time.pt.idx]] <- time.pt.spec.adj.matrix
    }
    rm(time.pt.idx)
    
    if (unrolled.net.idx == 1) {
      
      output.unrolled.net <- unrolled.net.adj.matrices
      
    } else {
      
      for (time.pt.idx in 1:length(output.unrolled.net)) {
        output.unrolled.net[[time.pt.idx]] <- (output.unrolled.net[[time.pt.idx]] + 
                                                 unrolled.net.adj.matrices[[time.pt.idx]])
      }
      rm(time.pt.idx)
      
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
      series.combn.threshold <- num.time.series
    }
  } else if (is.numeric(series.combn.method)) {
    ## Insert an edge in the final lrolled net iff it is present in at least 
    ## (series.combn.method * length(unrolled.net.file.names)) number of 
    ## interim rolled nets
    
    if ((series.combn.method > 0) & (series.combn.method < 1)) {
      series.combn.threshold <- (series.combn.method * num.time.series)
    } else {
      # print('\'series.combn.method\' accepts numeric values in the interval (0,1)')
      stop('\'series.combn.method\' accepts numeric values in the interval (0,1)')
    }
  }
  ##------------------------------------------------------------
  ## End: Evaluate the time-series combination threshold
  ##------------------------------------------------------------
  
  ## Apply the time-series combination threshold on the unrolled net adjacency matrices
  for (time.pt.idx in 1:length(output.unrolled.net)) {
    time.pt.spec.adj.matrix <- output.unrolled.net[[time.pt.idx]]
    time.pt.spec.adj.matrix[time.pt.spec.adj.matrix < series.combn.threshold] <- 0
    time.pt.spec.adj.matrix[time.pt.spec.adj.matrix != 0] <- 1
    
    output.unrolled.net[[time.pt.idx]] <- time.pt.spec.adj.matrix
    
  }
  rm(time.pt.idx)
  
  return(output.unrolled.net)
}
###########################################################################################