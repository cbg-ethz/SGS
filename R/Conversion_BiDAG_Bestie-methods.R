##################
# transform the CPTs from Bestie to SGS


#' CPTs Conversion Bestie to SGS format
#'
#' Converts CPTs of Bestie format to SGS format
#'
#' @param DAGaugmented CPTs of type Bestie (created by DAGparameters())
#' @return CPTs of SGS format
#' @export
convertCPTsBestieToSGS <- function(DAGaugmented){
  # Convert CPTs of Bestie format to SGS format
  # CPTs of Bestie are created by DAGparameters()
  
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
#' CPTs of Bestie are created by DAGparameters() 
#' and DAG of Bestie is created by BiDAG::iterativeMCMC()$DAG
#'
#' @param DAGaugmented Bayes net of type Bestie (created by DAGparameters())
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

#' Learn Bayesian Network
#'
#' Outputs the Bayesian network DAG and CPTs
#'
#' @param mydata input data for learning
#' @param bdepar hyper-parameters for structure learning
#' @return Bayesian network of type SubGroupSeparation
#' @export 
#' 
#' @importFrom BiDAG scoreparameters iterativeMCMC
#' @import Bestie
learnBN <- function(mydata, bdepar = list(chi = 0.5, edgepf = 2)){
  
  dataScore <- BiDAG::scoreparameters("bde", mydata, bdepar)
  
  learnedDAG <- BiDAG::iterativeMCMC(dataScore, mergetype = "skeleton", verbose = FALSE)
  
  DAGaugmented <- Bestie:::DAGparameters(learnedDAG$DAG, dataScore)
  
  BayesNetSGS <- convertBNBestieToSGS(DAGaugmented)
  
  return(BayesNetSGS)
}


# #' Learn CPTs (taken from Bestie package)
# #'
# #' Return the CPTs for a Bayesian network
# #'
# #' @param incidence DAG
# #' @param dataParams score parameters
# #' @return CPTs of the Bayes net
# #' @export
# #' @import Bestie
# DAGparameters <- function(incidence, dataParams){
#   if (dataParams$type != "bde") {
#     stop("The implementation is currently only for the BDe score.")
#   }
#   if (dataParams$DBN) {
#     n <- dataParams$n
#     nsmall <- dataParams$nsmall
#     bgn <- dataParams$bgn
#     slices <- dataParams$slices
#     dataParams_first <- dataParams$firstslice
#     if (bgn > 0) {
#       reorder <- c(1:bgn + nsmall, 1:nsmall)
#       dataParams_first$data <- dataParams_first$data[, 
#                                                      reorder]
#       dataParams_first$d1 <- dataParams_first$d1[, reorder]
#       dataParams_first$d0 <- dataParams_first$d0[, reorder]
#     }
#     params_first <- DAGparameters(incidence[1:n, 1:n], dataParams_first)
#     dataParams_other <- dataParams$otherslices
#     reorder <- c(1:n + nsmall, 1:nsmall)
#     dataParams_other$data <- dataParams_other$data[, reorder]
#     dataParams_other$d1 <- dataParams_other$d1[, reorder]
#     dataParams_other$d0 <- dataParams_other$d0[, reorder]
#     params <- DAGparameters(incidence, dataParams_other)
#     allalphas <- params$alphas
#     allalphas[1:n] <- params_first$alphas
#     allbetas <- params$betas
#     allbetas[1:n] <- params_first$betas
#     allpmeans <- params$pmeans
#     allpmeans[1:n] <- params_first$pmeans
#     if (slices > 2) {
#       nbig <- n + nsmall * (slices - 1)
#       incidence_unroll <- matrix(0, nbig, nbig)
#       incidence_unroll[1:nrow(incidence), 1:ncol(incidence)] <- incidence
#       allalphas_unroll <- vector("list", nbig)
#       allbetas_unroll <- vector("list", nbig)
#       allpmeans_unroll <- vector("list", nbig)
#       allalphas_unroll[1:ncol(incidence)] <- allalphas
#       allbetas_unroll[1:ncol(incidence)] <- allbetas
#       allpmeans_unroll[1:ncol(incidence)] <- allpmeans
#       for (ii in 1:(slices - 2)) {
#         block_rows <- n - nsmall + 1:(2 * nsmall)
#         block_cols <- n + 1:nsmall
#         incidence_unroll[block_rows + nsmall * ii, block_cols + 
#                            nsmall * ii] <- incidence[block_rows, block_cols]
#         allalphas_unroll[block_cols + nsmall * ii] <- allalphas[block_cols]
#         allbetas_unroll[block_cols + nsmall * ii] <- allbetas[block_cols]
#         allpmeans_unroll[block_cols + nsmall * ii] <- allpmeans[block_cols]
#         if (bgn > 0) {
#           block_rows <- 1:bgn
#           incidence_unroll[block_rows, block_cols + nsmall * 
#                              ii] <- incidence[block_rows, block_cols]
#         }
#       }
#       incidence <- incidence_unroll
#       allalphas <- allalphas_unroll
#       allbetas <- allbetas_unroll
#       allpmeans <- allpmeans_unroll
#     }
#   }
#   else {
#     n <- nrow(incidence)
#     allalphas <- vector("list", n)
#     allbetas <- vector("list", n)
#     allpmeans <- vector("list", n)
#     for (j in 1:n) {
#       parentNodes <- which(incidence[, j] == 1)
#       tempResult <- DAGparametersCore(j, parentNodes, dataParams)
#       allalphas[[j]] <- tempResult$alphas
#       allbetas[[j]] <- tempResult$betas
#       allpmeans[[j]] <- tempResult$pmeans
#     }
#   }
#   posteriorParams <- list()
#   posteriorParams$DAG <- incidence
#   posteriorParams$alphas <- allalphas
#   posteriorParams$betas <- allbetas
#   posteriorParams$pmeans <- allpmeans
#   return(posteriorParams)
# }

# DAGparametersCore <- function(j, parentNodes, param) {
#   # Taken from Bestie package
#   
#   # this function computes the parameters and their posterior distribution
#   # for a given node with parent set and for the data
#   
#   switch(param$type,
#          "bge" = {
#            coreParams <- list()
#            lp <- length(parentNodes) # number of parents
#            if (lp > 0) {# otherwise no regression coefficients
#              df <- param$awpN - param$n + lp + 1
#              R11 <- param$TN[parentNodes, parentNodes]
#              R12 <- param$TN[parentNodes, j]
#              
#              R11inv <- solve(R11) # could avoid inversions, but here for simplicity
#              mb <- R11inv %*% R12 # mean part
#              divisor <- param$TN[j, j] - R12 %*% mb
#              
#              coreParams$mus <- as.vector(mb)
#              coreParams$sigmas <- as.numeric(divisor/df) * R11inv
#              coreParams$dfs <- df
#            } else {
#              coreParams$mus <- NA
#              coreParams$sigmas <- NA
#              coreParams$dfs <- NA
#            }
#            return(coreParams)
#          },
#          "bde" = {
#            lp <- length(parentNodes) # number of parents
#            noParams <- 2^lp # number of binary states of the parents
#            chi <- param$chi
#            
#            alphas <- rep(NA, noParams)
#            betas <- rep(NA, noParams)
#            
#            switch(as.character(lp),
#                   "0"={ # no parents
#                     N1 <- sum(param$d1[, j])
#                     N0 <- sum(param$d0[, j])
#                     NT <- N0 + N1
#                     alphas <- (N1 + chi/(2*noParams))
#                     betas <- (N0 + chi/(2*noParams))
#                   },
#                   "1"={ # one parent
#                     summys <- param$data[, parentNodes]
#                     for (i in 1:noParams-1) {
#                       totest <- which(summys==i)
#                       N1 <- sum(param$d1[totest, j])
#                       N0 <- sum(param$d0[totest, j])
#                       NT <- N0 + N1
#                       alphas[i+1] <- (N1 + chi/(2*noParams))
#                       betas[i+1] <- (N0 + chi/(2*noParams))
#                     }
#                   },
#                   { # more parents
#                     summys <- colSums(2^(c(0:(lp-1)))*t(param$data[, parentNodes]))
#                     N1s <- collectC(summys, param$d1[, j], noParams)
#                     N0s <- collectC(summys, param$d0[, j], noParams)
#                     NTs <- N1s + N0s
#                     alphas <- (N1s + chi/(2*noParams))
#                     betas <- (N0s + chi/(2*noParams))
#                   }
#            )
#            
#            coreParams <- list()
#            coreParams$alphas <- alphas
#            coreParams$betas <- betas
#            coreParams$pmeans <- alphas/(alphas + betas)
#            
#            return(coreParams)
#          }
#   )
# }

# #' @import Bestie
# collectC <- function(xs, ys, n) {
#   # Taken from Bestie package
#   
#   .Call('_Bestie_collectC', PACKAGE = 'Bestie', xs, ys, n)
# }

