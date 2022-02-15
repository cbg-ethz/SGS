#' Sub Group Sampling (SGS) in Bayesian Networks
#'
#' Outputs the normalizing constant given some evidence of a Bayesian network
#'
#' @param BayesNet Bayesian network
#' @param obs list containing the evidence nodes and associated values
#' @return relevant DAG
#' @export
sample.subGroupSampling <- function(BayesNet, obs, N_samples = 100, relevantSubNet = FALSE, plot = TRUE)
{ 
  # Importance sampling in Bayesian networks given the evidence
  # Returned is the normalzing constant
  
  ### Compute the norm const and the updated Bayesian network of the subgroups using the junction-tree algorithm and LBP
  tempSGSresults <- subGroupSampling.propagation(BayesNet, obs)
  updatedNet <- tempSGSresults$BayesNet
  sampledNodes <- tempSGSresults$sampledNodes
  sampledEvidenceNodes <- tempSGSresults$sampledEvidenceNodes
  tempNormConst <- tempSGSresults$tempNormConst
  sortUpdatedNet <- TRUE

  # define range of sampled variables
  vars1 <- sampledNodes
  vars2 <- sampledEvidenceNodes
  
  ###  stochastic sampling
  net <- BayesNet
  lengthBN <- length(net@cpts)
  variablesBN <- variables(net)
  Nparents <- sapply(c(1:lengthBN),function(x) length(dim(net@cpts[[x]]))-1)
  
  observed.vars <- obs$observed.vars
  # if (class(obs$observed.vars) == "character") # hope that the user gives coherent input...
  #   observed.vars <- sapply(observed.vars, function(x) which(x == variablesBN))
  observed.vals <- obs$observed.vals
  
  # sort CPTs according to dimnames
  nodeNamesCPT <- sapply(c(1:lengthBN),function(x) match(names(dimnames(net@cpts[[x]])),variablesBN))
  for(g0 in 1:lengthBN){
    nodePosition <- match(g0, nodeNamesCPT[[g0]])
    tempPerm <- c(1:length(nodeNamesCPT[[g0]]))
    tempPerm <- c(tempPerm[-nodePosition],tempPerm[nodePosition])
    net@cpts[[g0]] <- aperm(net@cpts[[g0]], perm = tempPerm)
  }
  nodeParents <- sapply(c(1:lengthBN),function(x) match(head(names(dimnames(net@cpts[[x]])),-1),variablesBN))
  nodeNamesCPT <- sapply(c(1:lengthBN),function(x) match(names(dimnames(net@cpts[[x]])),variablesBN))
  
  # sort the CPTs of the updatedNet according to the CPTs of the original net
  if (sortUpdatedNet){
    # sort CPTs according to dimnames
    nodeNamesCPT2 <- sapply(c(1:lengthBN),function(x) match(names(dimnames(updatedNet@cpts[[x]])),variablesBN))
    nodesFree <- setdiff(1:lengthBN, obs$observed.vars)
    for(g1 in nodesFree){
      tempPerm <- match(nodeNamesCPT[[g1]], nodeNamesCPT2[[g1]])
      updatedNet@cpts[[g1]] <- aperm(updatedNet@cpts[[g1]], perm = tempPerm)
    }
    nodeParents2 <- sapply(c(1:lengthBN),function(x) match(head(names(dimnames(updatedNet@cpts[[x]])),-1),variablesBN))
    
  }
  
  # # faster version, but might mix up the ordering
  # nodeParents <- sapply(c(1:lengthBN),function(x) get.allParents(dag(net),x))
  
  # initialize vector of starting variables
  curVec <- rep(1, lengthBN)
  curVec[observed.vars] <- observed.vals
  
  # get topological order for sampling order
  topolOrder <- topo_sort(graph_from_adjacency_matrix(dag(net), "directed"), mode = c("out"))
  
  # get relevant subnetwork
  if (relevantSubNet==FALSE){
    consideredNodes <- c(1:lengthBN)
  }else{
    consideredNodes <- get.allAncestors(dag(net), observed.vars) 
  }
  
  # intialize partition function
  partitionFuncMean <- 0
  partitionFunc <- rep(0, N_samples)
  partitionFuncTempAll <- rep(0, N_samples)
  
  for (j1 in 1:N_samples){
    
    # reset partition function
    funcP1 <- 0
    funcP2 <- 0
    # funcP1 <- 1
    # funcP2 <- 1
    
    for (j2 in vars1){
      tempProb <- get.probNode(updatedNet, curVec, j2, nodeParents2)
      tempProb2 <- get.probNode(net, curVec, j2, nodeParents)
      
      # sample
      curVecNode <- sample(1:length(tempProb), size = 1, prob=tempProb)
      curVec[j2] <- curVecNode
      tempProbNode <- tempProb[curVecNode]
      
      # calculate importance function (use log scale to handle small numbers)
      funcP1 <- funcP1 + log(tempProb[curVecNode])
      # funcP1 <- funcP1*mpfr(tempProb[curVecNode], 50)
      # calculate probability function for weighting (use log scale to handle small numbers)
      funcP2 <- funcP2 + log(tempProb2[curVecNode])
      # funcP2 <- funcP2*mpfr(tempProb2[curVecNode], 50)
    }
    
    for (j3 in vars2){
      # get the prob of the evidence nodes
      parentsAndNode <- c(nodeParents[[j3]],j3)
      curVecPar <- rbind(curVec[parentsAndNode])
      tempProb <- net@cpts[[j3]][curVecPar]
      
      # no sampling for evidence nodes
      tempProbNode <- tempProb
      
      # calculate probability function for weighting (use log scale to handle small numbers)
      funcP2 <- funcP2 + log(tempProbNode)
      # funcP2 <- funcP2*mpfr(tempProbNode,50)
    }
    
    # calculate partition function (use log scale to handle small numbers)
    
    partitionFuncTemp <- exp(funcP2 - funcP1)*tempNormConst
    # partitionFuncTemp <- as.numeric(funcP2/funcP1)
    
    # caluclate mean of partition function
    partitionFuncMean <- partitionFuncMean*(j1-1)/j1 + partitionFuncTemp/j1
    
    partitionFunc[j1] <- partitionFuncMean
    partitionFuncTempAll[j1] <- partitionFuncTemp
  }
  
  results <- list("NormConst"=partitionFunc[N_samples],"partitionFunc"=partitionFunc,"partitionFuncTempAll"=partitionFuncTempAll)
  
  if (plot){
    plot.normConstCalc(results)
    #plot(1:N_samples, partitionFunc, type="l", col="blue", main="Parition Function", xlab="Sampling Step", ylab="Normalizing Constant")
  }
  
  return(results)
  
}

plot.DAG <- function(DAG){
  # plot DAG from matrix
  am.graph<-new("graphAM", adjMat=DAG, edgemode="directed")
  Rgraphviz::plot(am.graph)
}


plot.graph <- function(adjMatrix, edgemode="undirected"){
  # plot (undirected) graph from matrix
  am.graph<-new("graphAM", adjMat=adjMatrix, edgemode=edgemode)
  Rgraphviz::plot(am.graph)
}

