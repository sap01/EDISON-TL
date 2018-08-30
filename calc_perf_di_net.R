## Goal: Calculate performance metrics of predicted rolled or unrolled 
## Dynamic Bayesian Nets (DBNs).
##
###########################################################################
## Goal: Calculate performance metrics of the predicted directed net 
## adjacency matrix 'predicted.net.adj.matrix' w.r.t. 
## the true directed net adjacency matrix 'true.net.adj.matrix'.
###########################################################################
CalcPerfDiNet <-function(predicted.net.adj.matrix, true.net.adj.matrix)
{
  if ((nrow(predicted.net.adj.matrix) != nrow(true.net.adj.matrix)) | 
    (ncol(predicted.net.adj.matrix) != ncol(true.net.adj.matrix))) {
    stop('Dimension of the predicted net adj matrix and that of the 
         true net adj matrix do not match.')
  }
  
  num.nodes <- nrow(predicted.net.adj.matrix)
  
  TrPos = 0
  TrNeg = 0
  FlPos = 0
  FlNeg = 0
  
  for(i in 1:num.nodes)
  {
    for(j in 1:num.nodes)
    {
      if(predicted.net.adj.matrix[i,j] == 1 )
      {
        if(predicted.net.adj.matrix[i,j] == true.net.adj.matrix[i,j])
        {
          #True Positive
          TrPos = TrPos + 1	 	
        }
        else
        {
          FlPos = FlPos +1
        }
      }	   
      if(predicted.net.adj.matrix[i,j] == 0 )
      {
        if (predicted.net.adj.matrix[i,j] == true.net.adj.matrix[i,j])
        {
          # True Negative
          TrNeg = TrNeg + 1
        }
        else
        {
          FlNeg = FlNeg +1
        }
        
      }
    }	
  }
  
  #------------------------------------------------------------
  ## Begin: Create the format for result
  #------------------------------------------------------------
  Result <- matrix(0, nrow = 1, ncol = 11)
  colnames(Result) <- list('TP', 'TN', 'FP', 'FN', 'TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F1')
  #------------------------------------------------------------
  ## End: Create the format for result
  #------------------------------------------------------------
  
  
  #------------------------------------------------------------
  # Begin: Calculate Performance Metrics
  #------------------------------------------------------------
  TPR <- TrPos/(TrPos + FlNeg)
  FPR <- FlPos/(FlPos + TrNeg)
  FDR <- FlPos/(FlPos + TrPos)
  
  PPV <- TrPos/(TrPos + FlPos)
  ## If there is no positive i.e. 
  ## TrPos = FlPos = 0
  if (is.nan(PPV)) {
    PPV <- 0
  }
  
  ACC <- (TrPos + TrNeg)/(TrPos + FlPos + TrNeg + FlNeg)
  MCC <- ((TrPos * TrNeg) - (FlNeg * FlPos)) / sqrt((TrPos + FlPos) * (TrPos + FlNeg) * (TrNeg + FlPos) * (TrNeg+FlNeg))
  
  ## Begin: Evaluate F1 score
  F1 <- NULL
  if ((PPV + TPR) != 0) {
    F1 <- (2 * PPV * TPR) / (PPV + TPR)
  } else {
    ## Ref: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
    F1 <- 0
  }
  ## End: Evaluate F1 score
  
  ## Calculate AUC under ROC
  # table <- minet::validate(predicted.net.adj.matrix, true.net.adj.matrix)
  # AUC <- minet::auc.roc(table)
  #------------------------------------------------------------
  # End: Calculate Performance Metrics
  #------------------------------------------------------------
  
  Result[1, 1] <- TrPos
  Result[1, 2] <- TrNeg
  Result[1, 3] <- FlPos
  Result[1, 4] <- FlNeg
  Result[1, 5] <- TPR
  Result[1, 6] <- FPR
  Result[1, 7] <- FDR
  Result[1, 8] <- PPV
  Result[1, 9] <- ACC
  Result[1, 10] <- MCC
  Result[1, 11] <- F1
  # Result[1,8] <- AUC
  
  return (Result)
}
###########################################################################

###########################################################################
## Goal: Calculate performance metrics of the predicted unrolled net 
## adjacency matrices 'unrolled.net.adj.matrices' w.r.t. 
## the true unrolled net adjacency matrix 'true.net.adj.matrices'.
###########################################################################
CalcPerfDiNetUnrolled <- function(unrolled.net.adj.matrices, true.net.adj.matrices) {
  
  #------------------------------------------------------------
  ## Begin: Create the format for result
  #------------------------------------------------------------
  Result <- matrix(0, nrow = 1, ncol = 11)
  colnames(Result) <- list('TP', 'TN', 'FP', 'FN', 'TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F1')
  #------------------------------------------------------------
  ## End: Create the format for result
  #------------------------------------------------------------
  
  for (net.idx in 1:length(unrolled.net.adj.matrices)) {
    
    predicted.net.adj.matrix <- unrolled.net.adj.matrices[[net.idx]]
    
    ResultVsTrue <- CalcPerfDiNet(predicted.net.adj.matrix, true.net.adj.matrices[[net.idx]])
    Result <- rbind(Result, matrix(ResultVsTrue[1, ], nrow = 1, ncol = ncol(Result)))
    
    # rm(ResultVsTrue)
  }
  rm(net.idx)
  
  ## Print mean performance averaged over all time-varying networks
  ResultVsTrue <- colMeans(Result)
  ResultVsTrue <- matrix(colMeans(Result), nrow = 1, ncol = ncol(Result))
  colnames(ResultVsTrue) <- colnames(Result)
  # writeLines('Result EDISON vs True = \n')
  print(ResultVsTrue)
  rm(ResultVsTrue)
  
  return(Result)
}
###########################################################################