#!/usr/bin/env Rscript
## Goal of this program: Apply the TVDBN methods available in package 'EDISON' [1, 2].
## on the given datasets.
##
## Author: Saptarshi Pyne (saptarshipyne01@gmail.com)
## Last modified on: Mar 31, 2018
##
## The coding practices followed are documented at:
## https://github.com/sap01/coding_practices/blob/master/R_coding_practices.md
##
## Begin: References
## 1. Dondelinger, Frank, Sophie L?bre, and Dirk Husmeier. "Non-homogeneous 
## dynamic Bayesian networks with Bayesian regularization for inferring 
## gene regulatory networks with gradually time-varying structure." 
## Machine Learning 90.2 (2013): 191-230.
##
## 2. The 'EDISON' package in R: https://CRAN.R-project.org/package=EDISON
## End: References
##
## Usage:
## For Unix-alike OSes:
## Let us assume that this script is inside directory 
## '/home/saptarshi/R/R-3.3.2/projects/repoedisonr' and
## the Rscript file is inside directory '/home/saptarshi/R/R-3.3.2/bin'.
## Then, execute this script using the following commands (the '$' symbol
## represents the Bash command prompt):
## $ cd /home/saptarshi/R/R-3.3.2/projects/repoedisonr/  
## $ nohup time /home/saptarshi/R/R-3.3.2/bin/Rscript /home/saptarshi/R/R-3.3.2/projects/repoedisonr/EDISON.R input.json &
## where '/repoedisonr/asset/input.json' contains the user-defined parameters. A file 
## named 'nohup.out' will be generated inside 
## '/home/saptarshi/R/R-3.3.2/projects/repoedisonr/'.
##
## For Windows OSes:
## Let us assume that this script is inside directory 'D:\R\R-3.3.2\projects\repoedisonr' and
## the 'Rscript.exe' file is inside directory 'C:\Program Files\R\R-3.3.1\bin'.
## Then, execute this script using the following commands (the '>' symbol
## represents the DOS command prompt):
## >cd "D:\R\R-3.3.2\projects\repoedisonr"
## >"C:\Program Files\R\R-3.3.1\bin\Rscript.exe" EDISON.R input.json
## where '/repoedisonr/asset/input.json' contains the user-defined parameters.
##
## Input: A time series gene expression dataset with multiple time series.
## TODO (sap)
##
## Output: Time-varying Gene Regulatory Networks and a corresponding rolled up network.
##
## Remove all objects in the current workspace
rm(list = ls())

##------------------------------------------------------------
## Begin: Load the Required Packages
##------------------------------------------------------------
## For reading from and writing to '.json' files
library(rjson)

library(EDISON)
##------------------------------------------------------------
## End: Load the Required Packages
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Read User-defined input Params
##------------------------------------------------------------
input.args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1)
{
  stop("Exactly one input file must be supplied.", call.=FALSE)
}

input.params <- rjson::fromJSON(file = paste(getwd(), 'asset', input.args, sep = '/'))
rm(input.args)

## Input file for time-series gene expression data
input.data.filename <- input.params$input.data.filename
input.data.filename <- paste(getwd(), 'asset', input.data.filename, sep = '/')

## Number of time points (T)
num.timepts <- input.params$num.timepts

## Number of time series (S)
num.time.series <- input.params$num.time.series

## Please see the 'information.sharing' input param of the 
## EDISON::EDISON.run() function.
info.sharing.prior <- input.params$info.sharing.prior

## Num of iterations for MCMC simulation.
## Please see the 'num.iter' input param of the 
## EDISON::EDISON.run() function.
num.iters <- input.params$num.iters

rm(input.params)
##------------------------------------------------------------
## End: Read User-defined input Params
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Create the output directory
##------------------------------------------------------------
init.path <- getwd()

## Output directory name
output.dirname <- paste('output', format(Sys.time(), "%Y%m%d%H%M%S"), sep = '')

if(.Platform$OS.type == 'windows') {
  
  if(! output.dirname %in% shell("ls asset" , intern = TRUE)) {
    
    ## Output directory name for Windows OSes
    output.dirname <- paste('asset', output.dirname, sep = '/')
    output.dirname <- paste(init.path, output.dirname, sep = '/')
    
    ## Convert directory path to canonical form for the Windows OS.
    ## It raises the warning if the directory does not exist, which
    ## is expected. Therefore, please ignore the warning.
    output.dirname <- normalizePath(output.dirname, winslash = '\\', mustWork = NA)
    
    shell(paste('mkdir ', output.dirname, sep = ''), intern = TRUE, mustWork =NA)
  }
} else {
  ## .Platform$OS.type != 'windows'
  
  if(.Platform$OS.type == 'unix') {
    if(! output.dirname %in% system("ls asset" , intern = TRUE))
    {
      output.dirname <- paste('asset', output.dirname, sep = '/')
      output.dirname <- paste(init.path, output.dirname, sep = '/')
      
      system(paste('mkdir ', output.dirname, sep = ''))
    }
  } 
} 
##------------------------------------------------------------
## End: Create the output directory
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Load the Required External Functions
##------------------------------------------------------------
##
##------------------------------------------------------------
## End: Load the Required External Functions
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Main program
##------------------------------------------------------------

## Print the output dir name in 'nohup.out'
print('The output directory name is:')
print(output.dirname)
print('') ## to append a blank line

## Save console output in a file named 'output.txt' inside the output directory.
output.filename <- paste(output.dirname, 'output.txt', sep = '/')
output.file.conn <- file(output.filename, open = "wt")
sink(output.file.conn)

##------------------------------------------------------------
## Begin: Read input data file
##------------------------------------------------------------

## Begin: Find file extension of the input data file. Only '.tsv' and '.RData'
## are allowed.
## Split the string at every '.' and consider the last substring as the 
## file extension.
input.data.filename.ext <- unlist(strsplit(input.data.filename, '[.]'))
## End: Find file extension of the input data file. Only '.tsv' and '.RData'
## are allowed.

input.data <- NULL

if (input.data.filename.ext[length(input.data.filename.ext)] == 'tsv') {
  input.data <- read.table(input.data.filename, header = TRUE, sep="\t")
  
  ## Remove first col i.e. the time point names
  input.data <- input.data[, -1]
  
} else if (input.data.filename.ext[length(input.data.filename.ext)] == 'RData') {
  ## Loads an object named input.data
  load(input.data.filename)
}

num.nodes <- ncol(input.data)
##------------------------------------------------------------
## End: Read input data
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Learn the EDISON model
##------------------------------------------------------------

## start the timer
start.time <- proc.time()

node.names <- colnames(input.data)

## Input data is required to be a matrix where rows = nodes and
## cols = time pts.
## Please see input param 'input' of function EDISON::EDISON.run().
input.data <- t(input.data)

## Apply EDISON on one time series at a time
for (time.series.idx in 1:num.time.series) {
  
  ## Input data of the current time series
  input.data.curr.series <- input.data[, 1:num.timepts]
  
  ## Remaining input data
  input.data <- input.data[, -(1:num.timepts)]
  
  ## Run the user-defined TVDBN method on the
  ## input data of the current time series
  edison.result <- EDISON::EDISON.run(input.data.curr.series,
                                      information.sharing = info.sharing.prior,
                                      num.iter = num.iters)
  save(edison.result, 
       file = paste(output.dirname, '/edison.result.', time.series.idx, '.RData', sep = ''))
  rm(input.data.curr.series)
  
  ## Calculate posterior probabilities of changepoints
  change.pts <- EDISON::calculateCPProbabilities(edison.result)
  save(change.pts, 
       file = paste(output.dirname, '/change.pts.', time.series.idx, '.RData', sep = ''))
  rm(change.pts)
  
  ## Calculate marginal posterior probabilities of edges in the network
  unrolled.net.adj.matrices <- EDISON::calculateEdgeProbabilities(edison.result)
  save(unrolled.net.adj.matrices, 
       file = paste(output.dirname, 
                    '/unrolled.net.adj.matrices.', time.series.idx, '.RData', sep = ''))
  rm(edison.result, unrolled.net.adj.matrices)
  
  print(paste('Series', time.series.idx, 'is completed.', sep = ' '))

}
rm(time.series.idx)
##------------------------------------------------------------
## End: Learn the EDISON model
##------------------------------------------------------------

## Stop the timer
elapsed.time <- (proc.time() - start.time)
writeLines('elapsed.time = \n')
print(elapsed.time)

## Close output to the 'console_output.txt'
sink()
close(output.file.conn)

## Save R session info in a file named 'sessionInfo.txt' inside the output directory.
output.filename <- paste(output.dirname, 'sessionInfo.txt', sep = '/')
output.file.conn <- file(output.filename, open = "wt")
rm(output.filename)
sink(output.file.conn)
sessionInfo()
sink()
close(output.file.conn)
##------------------------------------------------------------
## End: Main Program
##------------------------------------------------------------
