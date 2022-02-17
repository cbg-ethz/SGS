##################
# transform the CPTs from Bestie to SGS


#' CPTs Conversion Bestie to SGS format
#'
#' Converts CPTs of Bestie format to SGS format
#'
#' @param DAGaugmented CPTs of type Bestie (created by Bestie:::DAGparameters())
#' @return CPTs of SGS format
#' @export
convertCPTsBestieToSGS <- function(DAGaugmented){
  # Convert CPTs of Bestie format to SGS format
  # CPTs of Bestie are created by Bestie:::DAGparameters()
  
  DAGBestie <- DAGaugmented$DAG
  CPTsBestie <- DAGaugmented$pmeans
  CPTsSGS <- list()
  
  # fill up CPTs node by node
  for(node in 1:dim(DAGBestie)[1]){
    
    # get parents of the node
    parentNodes <- get.allParents(DAGBestie, node)
    n_par <- length(parentNodes)
    
    # create empty table for CPT of SGS
    CPTsSGS[[node]] <- as.table(array(NA , dim = rep(2,n_par+1)))
    
    if(n_par==0){
      CPTsSGS[[node]][1:2] <- c(1-CPTsBestie[[node]],CPTsBestie[[node]])
    }else{
      # use permutations() to represent the possible states of the parents
      # for n_par>1, flip the columns to be in line with Bestie arrangement
      if(n_par==1){
        allVecPar <- permutations(2,n_par,repeats.allowed = TRUE)
      }else{
        allVecPar <- permutations(2,n_par,repeats.allowed = TRUE)[,n_par:1]
      }
      
      # some dirty coding to fill up the CPTs of SGS with the Bestie vectors
      for (i1 in 1:length(allVecPar[,1])){
        curVecPar <- allVecPar[i1,]
        mystring <- paste("CPTsSGS[[node]][", paste(as.character(curVecPar), collapse=", "),",] <- c(1-CPTsBestie[[node]][i1],CPTsBestie[[node]][i1])")
        eval(parse(text=mystring))
      }
    }
    # label the dimensions
    for(m1 in 1:(n_par+1)){
      dimnames(CPTsSGS[[node]])[[m1]] <- c(1:2)
    }
    names(dimnames(CPTsSGS[[node]])) <- as.character(c(parentNodes, node)) # dimnames = tempParents
  }
  return(CPTsSGS)
}

#' Bayesian network Conversion Bestie to SGS format
#'
#' Converts Bayesian network of Bestie format to SGS format
#' CPTs of Bestie are created by Bestie:::DAGparameters() 
#' and DAG of Bestie is created by BiDAG::iterativeMCMC()$DAG
#'
#' @param DAGaugmented Bayes net of type Bestie (created by Bestie:::DAGparameters())
#' @return Bayes net of SGS format
#' @export
convertBNBestieToSGS <- function(DAGaugmented){

  # get BN information from Bestie
  tempDAGSGS <- unname(DAGaugmented$DAG)
  tempNamesSGS <- as.character(1:dim(tempDAGSGS)[1])
  tempCPTsSGS <- convertCPTsBestieToSGS(DAGaugmented)
  N_nodes <- length(tempNamesSGS)
  
  # create new BN
  tempBN <- BN()
  name(tempBN) <- "converted Bayes net"
  num.nodes(tempBN) <- N_nodes
  variables(tempBN) <- tempNamesSGS
  discreteness(tempBN) <- rep(TRUE, N_nodes)
  node.sizes(tempBN) <- rep(2, N_nodes)
  cpts(tempBN) <- tempCPTsSGS
  dag(tempBN) <- tempDAGSGS
  wpdag(tempBN) <- matrix(0, N_nodes, N_nodes)
  
  return(tempBN)
}
