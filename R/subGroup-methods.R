#' Get relevant DAG nodes (ancestors)
#'
#' Outputs the relevant subnetwork of the Bayesian network given the evidence
#'
#' @param DAG DAG as permutation matrix
#' @param evidenceNodes vector containing the evidence nodes (order: 1 instead of "A")
#' @return relevant DAG
#' @export
get.allAncestors <- function(dag,A)
{ # identify ancestors of node A in DAG
  # A as order and DAG as matrix
  ancest<-c();
  rv<-c(A);
  
  # search for immediate ancestors (parents)
  for(i0 in 1:length(A))
  {
    ancest<-c(ancest, get.allParents(dag, A[i0]))
  }
  if(length(ancest)>0)
  { rv<-c(A,get.allAncestors(dag,ancest));
  }
  return(unique(rv));
}

#' Get Conditionally Independent Subgroup (CIS)
#'
#' Outputs the corresponding CIS of a specific node in the Bayesian network given the evidence
#'
#' @param DAG Bayesian network of type bn.fit
#' @param evidenceNodes vector containing the evidence nodes
#' @param startNode node contained in the subgroup
#' @return all nodes of the corresponding CIS
#' @export
get.subGroup <- function(DAG, evidenceNodes, startNode)
{
  tempGroupNodes <- startNode
  relevantEvidenceNodes <- c()
  allTempPar <- c()
  allTempChi <- c()
  
  jjj <- 1
  
  repeat{
    currentNode <- tempGroupNodes[jjj]
    
    # get parents, children and parents of the children (Markov Blanket)
    tempPar <- get.allParents(DAG,currentNode)
    tempChi <- get.allChildren(DAG,currentNode)
    tempParChi <- c()
    for (i in tempChi){
      tempParChi <- c(tempParChi, get.allParents(DAG,i))
    }
    tempParChi <- setdiff(tempParChi,currentNode)

    # store them
    allTempPar <- c(allTempPar, tempPar)
    allTempChi <- c(allTempChi, tempChi)
    
    newGroupNodes <- c(tempPar,tempChi,tempParChi)
    
    relevantEvidenceNodes <- c(relevantEvidenceNodes,intersect(newGroupNodes,evidenceNodes))
    
    relevantEvidenceNodes <- unique(relevantEvidenceNodes)
    
    newGroupNodes <- setdiff(newGroupNodes,evidenceNodes)
    
    newGroupNodes <- setdiff(newGroupNodes,tempGroupNodes)
    
    tempGroupNodes <- c(tempGroupNodes,newGroupNodes)
    
    jjj <- jjj+1
    if (jjj>length(tempGroupNodes)){
      break
    }
  }
  
  parEvidence <- intersect(allTempPar, relevantEvidenceNodes)
  chiEvidence <- intersect(allTempChi, relevantEvidenceNodes)
  
  return(list("tempGroupNodes"=tempGroupNodes, "relevantEvidenceNodes"=relevantEvidenceNodes, "parEvidence"=parEvidence, "chiEvidence"=chiEvidence))
}



#' Get All Conditionally Independent Subgroups (CIS)
#'
#' Outputs all CISs of the Bayesian network given the evidence
#'
#' @param BayesNet Bayesian network of type bn.fit
#' @param evidenceNodes vector containing the evidence nodes
#' @return list of all CIS of the Bayesian network given the evidence
#' @export
get.allSubGroups <- function(DAG, evidenceNodes, visualize = F){
  
  lengthDAG <- dim(DAG)[1]
  
  # reduce to relevant Subnetwork
  relevantNodes <- get.allAncestors(DAG,evidenceNodes)
  subDAG <- as.matrix(as_adjacency_matrix(induced_subgraph(
    graph_from_adjacency_matrix(DAG,mode="directed"),relevantNodes)))
  nodesOrder <- intersect(c(1:lengthDAG), relevantNodes)
  evidenceNodesNew <- match(evidenceNodes, nodesOrder)
  
  # DAGgraph <- graph_from_adjacency_matrix(DAG,mode="directed")
  # DAGsubgraph <- induced_subgraph((DAGgraph,relevantNodes)
  # subDAG <- as_adjacency_matrix(DAGsubgraph)
  
  # check if subgroups are available
  if(length(setdiff(relevantNodes, evidenceNodes))==0){
    print("No sampling required as the evidence nodes have no parents.")
    
    allSubGroups <- list(numeric(0))
    allEvidenceSubGroups <- list(evidenceNodes)
    allParEvidence <- list(numeric(0))
    allChiEvidence <- list(evidenceNodes)
    
    # If as individual subgroup:
    # allSubGroups <- list()
    # allEvidenceSubGroups <- list()
    # allParEvidence <- list()
    # allChiEvidence <- list()
    # 
    # for (p0 in 1:length(evidenceNodes)){
    #   allSubGroups[[p0]] <- numeric(0)
    #   allEvidenceSubGroups[[p0]] <- evidenceNodes[p0]
    #   allParEvidence[[p0]] <- numeric(0)
    #   allChiEvidence[[p0]] <- evidenceNodes[p0]
    # }

    # summarize results in returned list
    return(list("allSubGroups"=allSubGroups, "allEvidenceSubGroups"=allEvidenceSubGroups, 
                "allParEvidence"=allParEvidence, "allChiEvidence"=allChiEvidence))
  }
  
  # create variables
  lengthsubDAG <- dim(subDAG)[1]
  remain <- setdiff(c(1:lengthsubDAG),evidenceNodesNew)
  
  allSubGroups <- list()
  allEvidenceSubGroups <- list()
  allParEvidence <- list()
  allChiEvidence <- list()
  allIncludedChiEvidence <- c()
  
  # get subgroups
  ii <- 0
  repeat{
    
    startNode <- remain[1]
    
    subGroupResult <- get.subGroup(subDAG, evidenceNodesNew, startNode)
    
    ii <- ii+1
    allSubGroups[[ii]] <- nodesOrder[subGroupResult$tempGroupNodes]
    allEvidenceSubGroups[[ii]] <- nodesOrder[subGroupResult$relevantEvidenceNodes]
    allParEvidence[[ii]] <- nodesOrder[subGroupResult$parEvidence]
    allChiEvidence[[ii]] <- nodesOrder[subGroupResult$chiEvidence]
    
    allIncludedChiEvidence <- c(allIncludedChiEvidence, subGroupResult$chiEvidence)
    
    remain <- setdiff(remain,subGroupResult$tempGroupNodes)
    
    if(length(remain)==0){
      break
    }
  }
  
  # include left evidence nodes
  leftEvidenceNodes <- setdiff(evidenceNodesNew, allIncludedChiEvidence)
  if (!is.integer0(leftEvidenceNodes)){
    leftEvidenceParents <- c()
    for (ll in leftEvidenceNodes){
      leftEvidenceParents <- c(leftEvidenceParents, get.allParents(subDAG, ll))
    }
    leftEvidenceParents <- setdiff(leftEvidenceParents,leftEvidenceNodes)
    
    allLeftEvidenceBNnodes <- c(leftEvidenceNodes, leftEvidenceParents)
    
    ii <- ii+1
    allSubGroups[[ii]] <- numeric(0)
    allEvidenceSubGroups[[ii]] <- nodesOrder[allLeftEvidenceBNnodes]
    allParEvidence[[ii]] <- nodesOrder[leftEvidenceParents]
    allChiEvidence[[ii]] <- nodesOrder[leftEvidenceNodes]
  }
  
  # visualize 
  if(visualize){
    dim2 <- ceiling(sqrt(length(allSubBNs)))
    dim1 <- ceiling(length(allSubBNs)/dim2)
    storeMar <- par("mar")
    par(mfrow=c(dim1,dim2),mar=c(1,1,1,1))
    for (jj in 1:length(allSubBNs)){
      # graphviz.plot(allSubBNs[[jj]])
      graphviz.plot(allSubBNs[[jj]], highlight = list(nodes = allEvidenceSubGroups[[jj]], col = "black", fill = "grey"))
    }
    par(mfrow = c(1,1))
    par(mar=storeMar)
  }
  
  # summarize results in returned list
  return(list("allSubGroups"=allSubGroups, "allEvidenceSubGroups"=allEvidenceSubGroups, 
              "allParEvidence"=allParEvidence, "allChiEvidence"=allChiEvidence))
  
}


get.subBN <- function(BayesNet,relevantNodes)
{
  # get the Bayes net of a subset of the nodes of the original one
  
  subBN <- BN()
  
  relNodes <- sort(relevantNodes)
  subBN@variables <- BayesNet@variables[relNodes]
  subBN@num.nodes <- length(subBN@variables)
  subBN@discreteness <- BayesNet@discreteness[relNodes]
  subBN@node.sizes <- BayesNet@node.sizes[relNodes]
  subBN@cpts <- BayesNet@cpts[relNodes]
  subDAG <- as.matrix(as_adjacency_matrix(induced_subgraph(
    graph_from_adjacency_matrix(BayesNet@dag,mode="directed"),relNodes)))
  subBN@cpts <- BayesNet@cpts[relNodes]
  subBN@dag <- subDAG
  subWPDAG <- as.matrix(as_adjacency_matrix(induced_subgraph(
    graph_from_adjacency_matrix(BayesNet@wpdag,mode="directed"),relNodes)))
  subBN@wpdag <- subWPDAG
  
  return(subBN)
}


sub_belief.propagation <- function(BayesNet, obs)
{
  # belief propagation by sub groups
  evidenceNodes <- obs[[1]]
  subGroups <- get.allSubGroups(dag(BayesNet), evidenceNodes)
  
  N_subGroups <- length(subGroups$allSubGroups)
  
  for (k0 in 1:N_subGroups)
  {
    if(!length(subGroups$allSubGroups[[k0]])==0){
      # create the sub BN with according observation
      tempSubGroups <- c(subGroups$allSubGroups[[k0]], subGroups$allEvidenceSubGroups[[k0]])
      tempSubGroups <- sort(tempSubGroups)
      tempSubBN <- get.subBN(BayesNet, tempSubGroups)
      tempObsVarsOrder <- subGroups$allEvidenceSubGroups[[k0]]
      orderObs <- match(tempObsVarsOrder, obs[[1]])
      tempObsVars <- match(tempObsVarsOrder,tempSubGroups)
      tempObsVals <- obs[[2]][orderObs]
      tempSubObs <- list("observed.vars" = tempObsVars, "observed.vals" = tempObsVals)
      
      # perform belief propagation on the sub BN
      tempSubEngine <- InferenceEngine(tempSubBN)
      tempSubEngine <- belief.propagation(tempSubEngine, tempSubObs)
      tempUpdatedNet <- tempSubEngine@updated.bn
      
      # copy the updated information (CPTs) back to the original Bayes Net
      cpts(BayesNet)[tempSubGroups] <- cpts(tempUpdatedNet)
    }# else{
      # print("To do: incorporate evidence in CPT (not necessary for sampling though).")
    # }
    
    ## TO DO: see above the comment. 
  }
  return(BayesNet)
}


subGroup_belief.propagation <- function(BayesNet, obs, group_limit = 15)
{
  # belief propagation by sub groups
  evidenceNodes <- obs[[1]]
  subGroups <- get.allSubGroups(dag(BayesNet), evidenceNodes)
  
  N_subGroups <- length(subGroups$allSubGroups)
  
  for (k0 in 1:N_subGroups)
  {
    if((!length(subGroups$allSubGroups[[k0]])==0) && (subGroups$allSubGroups[[k0]] <= group_limit)){
      # create the sub BN with according observation
      tempSubGroups <- c(subGroups$allSubGroups[[k0]], subGroups$allEvidenceSubGroups[[k0]])
      tempSubGroups <- sort(tempSubGroups)
      tempSubBN <- get.subBN(BayesNet, tempSubGroups)
      tempObsVarsOrder <- subGroups$allEvidenceSubGroups[[k0]]
      orderObs <- match(tempObsVarsOrder, obs[[1]])
      tempObsVars <- match(tempObsVarsOrder,tempSubGroups)
      tempObsVals <- obs[[2]][orderObs]
      tempSubObs <- list("observed.vars" = tempObsVars, "observed.vals" = tempObsVals)
      
      # perform belief propagation on the sub BN
      tempSubEngine <- InferenceEngine(tempSubBN)
      tempSubEngine <- belief.propagation(tempSubEngine, tempSubObs)
      tempUpdatedNet <- tempSubEngine@updated.bn
      
      # copy the updated information (CPTs) back to the original Bayes Net
      cpts(BayesNet)[tempSubGroups] <- cpts(tempUpdatedNet)
    }# else{
    # print("To do: incorporate evidence in CPT (not necessary for sampling though).")
    # }
    
    ## TO DO: see above the comment. 
  }
  return(BayesNet)
}

#' Create random Bayesian network
#'
#' Creates a random Bayesian networks where the DAG simulation is based on bnlearn
#' (Generates graphs whose node ordering is given by the order of the labels in the
#' nodes argument)
#'
#' @param N_nodes number of nodes of the network (30 by default)
#' @param N_neighbours expected number of neighbours per node (binary by default)
#' @param nodeDim number of possible assignments of each node (binary by default)
#' @param visualize If true, the network is plotted
#' @return random Bayesian network of type BN()
#' @export
randomBN <- function(N_nodes, N_neighbours = 2, nodeDim = 2, visualize = FALSE)
{
  # create random discrete Bayes net based on the randDAG function of "pcalg"
  
  # sample random DAG
  tempDAG <- unname(as(randDAG(N_nodes, N_neighbours, "regular", 
                               weighted = FALSE), "matrix"))
  
  # create probability tables
  AllBNprobs <- list()
  node.sizes <- rep(2,N_samples)
  for (m0 in 1:N_nodes){
    tempParents <- get.allParents(tempDAG, m0)
    tempParentsLength <- length(tempParents)+1
    
    # sample CPTs uniformly
    tempProb <- as.table(array(runif(nodeDim^tempParentsLength), 
                               dim = rep(nodeDim,tempParentsLength)))
    # #for homogeneous CPTs
    # tempProb <- as.table(array(1/nodeDim, dim = rep(nodeDim,tempParentsLength))) 
    
    # normalize CPTs
    pot <- tempProb
    dms <- c(tempParents, m0)
    if (length(dms) > 1){
      pot.bak <- pot
      dms.bak <- dms
      remaining <- (1:length(dms))[-which(dms == m0)]
      dms <- dms[remaining]
      pot <- apply(pot, remaining, sum)
      out <- divide(as.array(pot.bak),
                    c(1:length(dms.bak)),
                    as.array(pot),
                    c(1:(length(dms.bak)-1)), # out$vars,
                    node.sizes)
      tempProb <- out$potential
    }else{
      tempProb <- pot / sum(pot)
    }
    
    # label the dimensions
    for(m1 in 1:tempParentsLength){
      dimnames(tempProb)[[m1]] <- c(1:nodeDim)
    }
    names(dimnames(tempProb)) <- as.character(c(tempParents, m0)) # dimnames = tempParents
    
    # store in list
    AllBNprobs[[m0]] <- tempProb
  }
  names(AllBNprobs) <- as.character(c(1:N_nodes))

  # create new BN
  tempBN <- BN()
  name(tempBN) <- "random Bayes net"
  num.nodes(tempBN) <- N_nodes
  variables(tempBN) <- names(AllBNprobs)
  discreteness(tempBN) <- rep(TRUE, N_nodes)
  node.sizes(tempBN) <- rep(nodeDim, N_nodes)
  cpts(tempBN) <- AllBNprobs
  dag(tempBN) <- tempDAG
  wpdag(tempBN) <- matrix(0, N_nodes, N_nodes)
  
  # plot the DAG
  if (visualize==TRUE){
    plot(tempBN)
    plot.new()
  }
  
  return(tempBN)
}






# get.normConst <- function(results){
#   # estimate the normalizing constant
#   net <- results@bn
#   netNew <- results@updated.bn
#   
#   lengthBN <- length(net@variables)
#   
#   prob1 <- 0
#   prob2 <- 0
#   
#   vars1 <- setdiff(c(1:lengthBN),results@observed.vars)
#   vars2 <- results@observed.vars
#   
#   Nparents <- sapply(c(1:lengthBN),function(x) length(dim(net@cpts[[x]]))-1)
#   
#   for(i0 in vars1){
#     selectedElement <- 1
#     repeat{
#       if(!netNew@cpts[[i0]][selectedElement]==0) break
#       selectedElement <- selectedElement+1
#     }
#     prob1 <- prob1 + log(unname(net@cpts[[i0]][selectedElement]))
#     prob2 <- prob2 + log(unname(netNew@cpts[[i0]][selectedElement]))
#   }
#   
#   i2 <- 1
#   for(i1 in vars2){
#     selectedElement <- rbind(c(rep(1,Nparents[i1]), results@observed.vals[i2]))
#     prob1 <- prob1 + log(unname(net@cpts[[i1]][selectedElement]))
#     i2 <- i2+1
#   }
#   
#   normConst <- unname(exp(prob1-prob2))
#   
#   return(normConst)
# }





