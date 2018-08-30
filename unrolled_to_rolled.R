## Goal: Given an unrolled net, roll it into a rolled net.
##
#' @title Convert a given unrolled Dynamic Bayesian Network (DBN) into a rolled DBN using different rolling methods.
#' @param unrolled.net.adj.matrix.list Given time-varying network adjacency list. Its length = number of time-interval-specific nets. 
#' @param roll.method Which rolling method to use from {'any', 'all', or some real number in (0, 1), like - 0.5}.
#' @param allow.self.loop Boolean to decide whehter to allow self loop or not in the rolled DBN.
#' @return rolled.net.adj.matrix Return the rolled DBN adjacency matrix.
##
## Loading the Packages
## End(Loading the Packages)
##
UnrolledToRolled <- function(unrolled.net.adj.matrix.list, edge.post.prob.threshold, 
                             roll.method, allow.self.loop) {
  
  ## Number of nets or
  ## number of adjacency matrices
  num.nets <- length(unrolled.net.adj.matrix.list$probs.all) 
  
  ## Number of nodes per net
  num.nodes <- nrow(unrolled.net.adj.matrix.list$probs.all[[1]])
  
  ## Names of the ndoes
  node.names <- rownames(unrolled.net.adj.matrix.list$probs.all[[1]])
  
  ## Initialize rolled DBN adj matrix as a zero matrix
  rolled.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                                  dimnames = c(list(node.names), list(node.names)))
  
  for (time.pt.idx in 1:num.nets) {

    ## Time point specific adjacency matrix.
    ## Adjacency matrix at time pt 'time.pt.idx'.
    time.pt.spec.adj.matrix <- unrolled.net.adj.matrix.list$probs.all[[time.pt.idx]]

    ## Remove edges with marginal posterior probability less than
    ## 'edge.post.prob.threshold'
    time.pt.spec.adj.matrix[time.pt.spec.adj.matrix < edge.post.prob.threshold] <- 0

    ## Convert 'time.pt.spec.adj.matrix' to a binary matrix of zeros and ones
    time.pt.spec.adj.matrix[time.pt.spec.adj.matrix != 0] <- 1
    
    rolled.net.adj.matrix <- (rolled.net.adj.matrix + time.pt.spec.adj.matrix)
  }
  rm(time.pt.idx, edge.post.prob.threshold)
  
  ##------------------------------------------------------------
  ## Begin: Evaluate the rolling threshold
  ##------------------------------------------------------------
  roll.threshold <- NULL
  
  if(is.character(roll.method)) {
    if (roll.method == 'any') {
      ## Insert an edge in the rolled net iff it is present in at least one time-interval-specific net
      roll.threshold <- 1
    } else if (roll.method == 'all') {
      ## Insert an edge in the rolled net iff it is present in all time-interval-specific nets
      roll.threshold <- num.nets
    }
  } else if (is.numeric(roll.method)) {
    ## Insert an edge in the rolled net if it is present in at least 
    ## (roll.method * num.nets) number of time-interval-specific nets
    
    if ((roll.method > 0) & (roll.method < 1)) {
      roll.threshold <- num.nets * roll.method
    } else {
      # print('\'roll.method\' accepts numeric values in the interval (0,1)')
      stop('\'roll.method\' accepts numeric values in the interval (0,1)')
    }
  }
  ##------------------------------------------------------------
  ## End: Evaluate the rolling threshold
  ##------------------------------------------------------------
  
  
  # writeLines('\n rolled.net.adj.matrix = \n')
  # print(rolled.net.adj.matrix)
  
  ## Apply the rolling threshold on the rolled net adjacency matrix
  for (tgt.node.idx in 1:ncol(rolled.net.adj.matrix)) {
    for (src.node.idx in 1:nrow(rolled.net.adj.matrix)) {
      if (rolled.net.adj.matrix[src.node.idx, tgt.node.idx] >= roll.threshold) {
        rolled.net.adj.matrix[src.node.idx, tgt.node.idx] <- 1
      } else {
        rolled.net.adj.matrix[src.node.idx, tgt.node.idx] <- 0
      }
    }
    rm(src.node.idx)
  }
  rm(tgt.node.idx)
  
  ## Remove self loops if 'allow.self.loop' = FALSE
  if (!allow.self.loop) {
    diag(rolled.net.adj.matrix) <- 0
  }
  
  return(rolled.net.adj.matrix)
}