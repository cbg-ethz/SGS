#' Importance Sampling (IS) in Bayesian Networks
#'
#' Outputs the normalizing constant given some evidence of a Bayesian network
#'
#' @param BayesNet Bayesian network
#' @param obs List containing the evidence nodes and associated values
#' @param s_method Inference method
#' @param N_samples Number of samples
#' @param relevantSubNet Whether to consider the relevant subnet
#' @param plot If TRUE, plot results
#' @return Normalizing constant and vector of intermediate results
#' @export
#' @importFrom utils head
sample.normConst <- function(BayesNet, obs, s_method = "SGS", N_samples = 100, relevantSubNet = FALSE, plot = TRUE)
{ 
  # Importance sampling in Bayesian networks given the evidence
  # Returned is the normalzing constant
  
  # sub group sampling
  if(s_method=="SGS"){
    results <- sample.subGroupSampling(BayesNet, obs, N_samples = N_samples, plot = plot)
    return(results)
  }
  
  # Compute the updated Bayesian network using BP 
  sortUpdatedNet <- FALSE
  if (s_method=="BP-LW"){
    engineTemp <- InferenceEngine(BayesNet)
    engineTemp <- belief.propagation(engineTemp, obs)
    updatedNet <- engineTemp@updated.bn
    sortUpdatedNet <- TRUE
    updatedNetBP <- updatedNet
    # print(updatedNet)
  }
  # Compute the updated Bayesian network using LBP 
  if (s_method=="LBP-LW"){
    engineTemp <- InferenceEngine(BayesNet)
    engineTemp <- loopy_belief.propagation(engineTemp, obs)
    updatedNet <- engineTemp@updated.bn
    sortUpdatedNet <- TRUE
  }
  # Compute the updated Bayesian network using subgroup BP
  if (s_method=="SGSold"){
    updatedNet <- subGroup_belief.propagation(BayesNet, obs)
    sortUpdatedNet <- TRUE
    updatedNetSGS <- updatedNet
    # print(updatedNet)
  }
  # Compute the updated Bayesian network using subgroup BP
  if (s_method=="LBP"){
    updatedNet <- subGroup_loopyBelief.propagation(BayesNet, obs)
    sortUpdatedNet <- TRUE
  }
  
  net <- BayesNet
  lengthBN <- length(net@cpts)
  variablesBN <- variables(net)
  Nparents <- sapply(c(1:lengthBN),function(x) length(dim(net@cpts[[x]]))-1)
  
  observed.vars <- obs$observed.vars
  if (class(obs$observed.vars) == "character") # hope that the user gives coherent input...
    observed.vars <- sapply(observed.vars, function(x) which(x == variablesBN))

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
  minDimBN <- min(node.sizes(BayesNet)) # pick starting vector at random to avoid bias of Gibbs sampling
  curVec <- sample.int(minDimBN,lengthBN,replace = TRUE) #rep(1, lengthBN)
  curVec[observed.vars] <- observed.vals
  
  # get topological order for sampling order
  topolOrder <- topo_sort(graph_from_adjacency_matrix(dag(net), "directed"), mode = c("out"))
  
  # get relevant subnetwork
  if (relevantSubNet==FALSE){
    consideredNodes <- c(1:lengthBN)
  }else{
    consideredNodes <- get.allAncestors(dag(net), observed.vars) 
  }
  
  # define range of sampled variables
  vars1 <- intersect(topolOrder, setdiff(consideredNodes,observed.vars))
  vars2 <- intersect(topolOrder, observed.vars)
  
  # intialize partition function
  partitionFuncMean <- 0
  partitionFunc <- rep(0, N_samples)
  partitionFuncTempAll <- rep(0, N_samples)
  
  if(s_method=="GS"){
    nodeChildren <- sapply(c(1:lengthBN),function(x) get.allChildren(dag(net),x))
  }
  
  for (j1 in 1:N_samples){
    
    # reset partition function
    funcP1 <- 0
    funcP2 <- 0
    # funcP1 <- 1
    # funcP2 <- 1
    
    for (j2 in vars1){
      # get probabilities for the sampling
      if(s_method=="FS"){
        # get the prob. of node given the parents
        tempProb <- get.probNode(net, curVec, j2, nodeParents)
        tempProb2 <- tempProb
      }else if(s_method=="GS"){
        # get the prob. of node given the Markov blanket
        probsMB <- get.probNodeMB(net, curVec, j2, nodeParents, nodeChildren, consideredNodes)
        tempProb <- probsMB$probMB
        tempProb2 <- probsMB$tempProbNode
      }else if(s_method %in% c("BP-LW", "LBP-LW", "SGBP-LW","SGSold","LBP")){
        tempProb <- get.probNode(updatedNet, curVec, j2, nodeParents2)
        tempProb2 <- get.probNode(net, curVec, j2, nodeParents)
      }else{
        stop('s_method is not valid')
      }
      
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
    
    partitionFuncTemp <- exp(funcP2 - funcP1)
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


get.probNode <- function(net, curVec, node, nodeParents)
{
  # get probability of node given vector
  cpts <- cpts(net)
  parents <- nodeParents[[node]]
  if(length(parents)==0){
    tempProb <- cpts[[node]]
  }else{
    curVecPar <- curVec[parents]
    mystring <- paste("cpts[[node]][", paste(as.character(curVecPar), collapse=", "),",]")
    tempProb <- eval(parse(text=mystring))
  }
  
  return(tempProb)
}

get.probNodeMB <- function(net, curVec, node, nodeParents, nodeChildren, consideredNodes)
{
  # get probability of Markov blanket of node (used in Gibbs sampling)
  curVecTemp <- curVec
  
  tempProbNode <- get.probNode(net, curVec, node, nodeParents)
  
  tempChildren <- get.allChildren(dag(net), node)
  
  tempChildren <- intersect(tempChildren, consideredNodes)
  
  probMB <- rep(NA, length(tempProbNode))
  
  for (h2 in 1:length(tempProbNode)){
    curVecTemp[node] <- h2
    ProbChildAll <- 1
    for (h1 in tempChildren){
      tempProbChild <- get.probNode(net, curVecTemp, h1, nodeParents)
      ProbChildAll <- ProbChildAll*tempProbChild[curVec[h1]]
    }
    probMB[h2] <- ProbChildAll*tempProbNode[h2]
  }
  
  
  probMB <- probMB/sum(probMB)
  
  return(list("probMB"=probMB, "tempProbNode"=tempProbNode))
}


calcExactInference <- function(myBN, myObs, showTimes=FALSE){
  # calculate the exact normalizing constant of a Bayes net given observation
  
  start_time <- Sys.time() # measure computation time
  
  # consider only the ancestors of the observed variables
  relevantNodes <- get.allAncestors(myBN@dag, myObs$observed.vars)
  nodes_free <- setdiff(relevantNodes, myObs$observed.vars)
  
  # print approximate time
  if(showTimes==TRUE){
    print(paste("Approximated time:",3/(2^16)*2^length(nodes_free), "minutes"))
  }
  
  # create matrix containing all possible combinations in the net
  # this will then be used in the summation over possible values
  
  # create possible permuations of free nodes
  allVec <- permutations(2,length(nodes_free),c(1,2), repeats.allowed=TRUE)
  
  # create empty matrix
  newVec <- matrix(1,dim(allVec)[1],length(myBN@variables))
  
  # insert permutations into empty matrix  
  j <- 1
  for (i in nodes_free){
    newVec[,i] <- allVec[,j]
    j <- j+1
  }
  
  # insert observed values into empty matrix  
  j <- 1
  for (i in myObs$observed.vars){
    newVec[,i] <- myObs$observed.vals[j]
    j <- j+1
  }
  
  # calc the marginals
  nodeParents <- sapply(1:length(myBN@variables), function(x) get.allParents(myBN@dag, x))
  nodes <- relevantNodes

  N_length <- dim(newVec)[1]
  normConst <- 0
  for (i in 1:N_length){
    normConst <- normConst + get.probVector(myBN, newVec[i,], nodes, nodeParents)
  }
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  # print elapsed time
  if(showTimes==TRUE){
    print(paste("Elapsed time:",elapsed_time/60, "minutes"))
  }
  
  return(normConst)
  
}

get.probVector <- function(myBN, curVec, nodes, nodeParents, logSpace=FALSE){
  # get the product of the probability of multiple nodes given vector
  
  # work in log space to allow for small values (as they occur in high dims)
  tempVar <- 0
  for (iNode in nodes){
    curNodeVal <- curVec[iNode]
    tempVar <- tempVar + log(unname(
      get.probNode(myBN, curVec, iNode, nodeParents)[curNodeVal]))
  }
  
  # if necessary, return result in logspace
  if (logSpace==FALSE){
    tempVar <- exp(tempVar)
  }
  
  return(tempVar)
  
}




permutations <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{
  # function copied from the gtools v3.5.0 by Gregory R. Warnes
  
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
      0) 
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
      0) 
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE) 
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) 
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) 
    sub <- function(n, r, v) {
      if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                                                  ncol = ncol(inner), nrow = nrow(inner) * n, 
                                                  byrow = TRUE))
      }
    }
  else sub <- function(n, r, v) {
    if (r == 1) 
      matrix(v, n, 1)
    else if (n == 1) 
      matrix(v, 1, r)
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                                                        1, r - 1, v[-i])))
      X
    }
  }
  sub(n, r, v[1:n])
}




